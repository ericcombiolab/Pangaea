reads_type="10x"
# thread=150

root=`pwd`
conda_path="/home/comp/zmzhang/software/anaconda3"

#interleaved_reads=${root}/interleaved_link_reads.sorted.fastq
#metaspades_contigs=${root}/contigs_to_merge/metaspades_contigs.fasta
#athena_out=${root}/athena_out
interleaved_reads=`realpath $1`
metaspades_contigs=`realpath $2`
athena_out=$3
thread=$4
athena_thread=$thread
if [ $athena_thread -gt 128 ];then
    athena_thread=128
fi
echo $athena_thread


if [ ! -f $interleaved_reads ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` $interleaved_reads not exist"
    exit 1
fi
if [ ! -f $metaspades_contigs ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` $metaspades_contigs not exist"
    exit 1
fi

#Check if the conda env is xiaojin_pangaea
if [ $CONDA_DEFAULT_ENV != "xiaojin_pangaea" ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` Current conda env is not xiaojin_pangaea, activate it"
    source ${conda_path}/bin/activate xiaojin_pangaea
    if [ $CONDA_DEFAULT_ENV != "xiaojin_pangaea" ];then
        echo "`date "+%Y-%m-%d %H:%M:%S"` Activate xiaojin_pangaea failed"
        exit 1
    fi
fi


echo "`date "+%Y-%m-%d %H:%M:%S"` ============athena Pipeline start at `date`============"


echo "`date "+%Y-%m-%d %H:%M:%S"`  Branch 2 - step 3: run bwa index"
if [ -f ${metaspades_contigs}.bwt ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` already run bwa index"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` bwa index ${metaspades_contigs}"
    bwa index ${metaspades_contigs}
fi

echo "`date "+%Y-%m-%d %H:%M:%S"`  Branch 2 - step 4: run bwa mem"
if [ ! -d ${athena_out} ];then
    mkdir -p ${athena_out}
fi

athena_out=`realpath ${athena_out}`
if [ -f ${athena_out}/align-reads.contigs.bam ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` already run bwa mem"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` bwa mem -t $thread -C -p ${metaspades_contigs} ${interleaved_reads} | samtools sort -@ $thread -o ${athena_out}/align-reads.contigs.bam"
    bwa mem -t $thread -C -p ${metaspades_contigs} ${interleaved_reads} | samtools sort -@ $thread -o ${athena_out}/align-reads.contigs.bam
fi

echo "`date "+%Y-%m-%d %H:%M:%S"`  Branch 2 - step 5: run samtools index"
if [ -f ${athena_out}/align-reads.contigs.bam.bai ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` already run samtools index"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` samtools index -@ $thread ${athena_out}/align-reads.contigs.bam"
    samtools index -@ $thread ${athena_out}/align-reads.contigs.bam
fi


if [ ! -f ${athena_out}/config.json ];then
    echo "{
\"ctgfasta_path\": \"${metaspades_contigs}\",
\"reads_ctg_bam_path\": \"${athena_out}/align-reads.contigs.bam\",
\"input_fqs\": \"${interleaved_reads}\",
\"cluster_settings\": {
\"cluster_type\": \"multiprocessing\",
\"processes\": $athena_thread
}
}
" > ${athena_out}/config.json
fi


if [ ! -f ${athena_out}/results/olc/athena.asm.fa ]; then
    echo "`date "+%Y-%m-%d %H:%M:%S"`  Branch 2 - step 6: run athena"
    source ${conda_path}/bin/activate python2
    if [ $CONDA_DEFAULT_ENV == "python2" ];then
        echo "`date "+%Y-%m-%d %H:%M:%S"` python2 activated"
    else
        echo "`date "+%Y-%m-%d %H:%M:%S"` python2 not activated"
        exit
    fi
    cd ${athena_out}/
    echo "`date "+%Y-%m-%d %H:%M:%S"` athena-meta --config config.json"
    athena-meta --config config.json
    cd -
fi

# wait all steps to be done
wait



echo "`date "+%Y-%m-%d %H:%M:%S"` ============athena Pipeline end at `date`============"
