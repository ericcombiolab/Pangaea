reads_type="10x"
# thread=150

CURRENT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$CURRENT_DIR/../../

# Define the function to check and activate a Conda environment
activate_conda_env() {
    # The first argument is the target environment name
    local TARGET_ENV="$1"

    # Get the current active Conda environment
    local CURRENT_ENV="$CONDA_DEFAULT_ENV"

    # Check if the current environment is not the target environment
    if [[ "$CURRENT_ENV" != "$TARGET_ENV" ]]; then
        echo "The current Conda environment is not: $TARGET_ENV, it is: $CURRENT_ENV"
        echo "Activating the target environment: $TARGET_ENV"

        # Initialize Conda for this shell (required for using `conda activate` in scripts)
        eval "$(conda shell.bash hook)"
        
        # Activate the target Conda environment
        conda activate "$TARGET_ENV"

        # Verify if the activation was successful
        if [[ "$CONDA_DEFAULT_ENV" == "$TARGET_ENV" ]]; then
            echo "The target environment $TARGET_ENV has been successfully activated."
        else
            echo "Failed to activate the target environment $TARGET_ENV. Please check your Conda configuration."
            return 1
        fi
    else
        # If the current environment is already the target
        echo "The current Conda environment is already: $TARGET_ENV"
    fi
}


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

# print command
echo "$0 $*"
pangaea_env="pangaea-test"
athena_env="athena-meta"
# Activate the Conda environment
activate_conda_env "$pangaea_env"

# #Check if the conda env is xiaojin_pangaea
# if [ $CONDA_DEFAULT_ENV != "xiaojin_pangaea" ];then
#     echo "`date "+%Y-%m-%d %H:%M:%S"` Current conda env is not xiaojin_pangaea, activate it"
#     source ${conda_path}/bin/activate xiaojin_pangaea
#     if [ $CONDA_DEFAULT_ENV != "xiaojin_pangaea" ];then
#         echo "`date "+%Y-%m-%d %H:%M:%S"` Activate xiaojin_pangaea failed"
#         exit 1
#     fi
# fi


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
    # Activate the athena-meta Conda environment
    activate_conda_env "$athena_env"
    cd ${athena_out}/
    echo "`date "+%Y-%m-%d %H:%M:%S"` athena-meta --config config.json"
    athena-meta --config config.json
    cd -
fi

# wait all steps to be done
wait



echo "`date "+%Y-%m-%d %H:%M:%S"` ============athena Pipeline end at `date`============"
