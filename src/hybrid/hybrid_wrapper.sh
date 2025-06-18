#!/bin/bash
set -o pipefail
set -e
echo "`date` starting"
CURRENT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$CURRENT_DIR/../../
BINDIR=$root/src/bin
# the file path of this script
cd $root/src/cpptools/;cmake .;make ;cd -

while getopts "l:r:R:i:t:a:A:o:" opt; do
    case $opt in
        l) longreads="$OPTARG" ;;
        r) short_R1="$OPTARG" ;;
        R) short_R2="$OPTARG" ;;
        i) identity="$OPTARG" ;;
        t) type="$OPTARG" ;;
        a) athena_lc="$OPTARG" ;; # Optional Athena input file
        A) athena_out="$OPTARG" ;; # Optional Athena output file
        o) output_dir="$OPTARG" ;; # Optional output directory
        *) 
            echo "Usage: $0 -l <longreads> -r <short_R1> -R <short_R2> [-i <identity>] [-t <type>] [-a <athena_lc>] [-A <athena_out>] [-o <output_dir>]"
            exit 1
            ;;
    esac
done

# Check if the output directory exists, if not create it

if [ -z "$output_dir" ]; then
    echo "Warning: Not specify output directory, using default value"
    output_dir="$CURRENT_DIR/default_hybrid_out"
fi

if [ ! -d $output_dir ];then
    mkdir -p $output_dir
    echo "Output directory $output_dir created."
fi

# Check mandatory arguments
if [ -z "$longreads" ] || [ -z "$short_R1" ] || [ -z "$short_R2" ]; then
    echo "Error: Missing mandatory arguments. Usage: $0 -l <longreads> -r <short_R1> -R <short_R2> [-i <identity>] [-t <type>] [-a <athena_lc>] [-A <athena_out>] [-o <output_dir>]"
    exit 1
fi

# Set default values for optional arguments
if [ -z "$identity" ]; then
    echo "Warning: Not specify identity, use 60"
    identity=60
fi

if [ -z "$type" ]; then
    echo "Warning: Not specify type, use hybridspades"
    type="hybridspades"
fi

if [ -z "$athena_lc" ]; then
    echo "Warning: Not specify Athena input file, using default value"
    athena_lc="$output_dir/athena_out/results/olc/flye-input-contigs.fa"
fi
if [ -z "$athena_out" ]; then
    echo "Warning: Not specify Athena output file, using default value"
    athena_out="$output_dir/athena_out/results/olc/athena.asm.fa"
fi


barcode_list=$CURRENT_DIR/barcodelist.txt
threads=100

# Generate intermediate files directory
intermediate_files=$output_dir'/intermediate_files'
# athena_lc=$output_dir'/athena_out/results/olc/flye-input-contigs.fa'
# athena_out=$output_dir'/athena_out/results/olc/athena.asm.fa'
# if [ -f contigs/athena.asm.fa ] && [ -f contigs/flye-input-contigs.fa ];then
#     athena_lc=$output_dir'/contigs/flye-input-contigs.fa'
#     athena_out=$output_dir'/contigs/athena.asm.fa'
# fi
echo "athena_lc: $athena_lc"
echo "athena_out: $athena_out"

# print command
echo "$0 $*"
pangaea_env="pangaea"

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


if [ $# -lt 3 ];then
    echo "Usage: $0 longreads short_R1 short_R2 [identity] [type]"
    exit 1
fi

if [ -d $BINDIR ];then
    echo "bin directory $BINDIR found"
else
    echo "bin directory $BINDIR not found"
    exit 1
fi

if [ -z $identity ];then
    echo "Warning: Not specify identity, use 60"
    identity=60
fi

if [ -z $type ];then
    echo "Warning: Not specify type, use operams"
    type="operams"
fi

if [ $identity -lt 50 ];then
    echo "identity should be greater than 50"
    exit 1
fi

if [ ! -f $longreads ]; then
    echo "longreads $longreads not found"
    exit 1
fi

if [ ! -f $short_R1 ]; then
    echo "short_R1 $short_R1 not found"
    exit 1
fi

if [ ! -f $short_R2 ]; then
    echo "short_R2 $short_R2 not found"
    exit 1
fi

if [[ $longreads == *.gz ]];then
    echo "unzip $longreads"
    if [ -f ${longreads%.gz} ];then
        rm ${longreads%.gz}
    fi
    gunzip -k $longreads
    longreads=${longreads%.gz}
fi

if [ ! -f $barcode_list ];then
    echo "barcode_list $barcode_list not found"
    exit 1
fi

if [ ! -d ${intermediate_files} ];then
    mkdir ${intermediate_files}
fi

if [ -f ${intermediate_files}/barcode_maps.txt ];then
    echo "${intermediate_files}/barcode_maps.txt already done"
else
    
    # get the barcode map
    (
        set -x;
        awk '{                      
            if(NR % 4==1) { printf("%s\n",$1);}
            }' $longreads | tr "\t" "\n" > ${intermediate_files}/longreads.names.txt
        longreads_length=`cat ${intermediate_files}/longreads.names.txt | grep "^@"| wc -l`
        echo "longreads sequence number: $longreads_length"
        paste -d ' ' <(cat ${intermediate_files}/longreads.names.txt | sed 's/^@//') <(head -n $longreads_length $barcode_list | awk '{printf("BX:Z:%s\n",$0)}') > ${intermediate_files}/barcode_maps.txt
    )
fi


if [ -f $longreads.bwt ];then
    echo "bwa index already done"
else
   echo "bwa index"
   bwa index $longreads > /dev/null 2>&1
fi

if [ -f ${intermediate_files}/short2long.bam ];then
    echo "bwa mem already done"
else
    bwa mem -t $threads $longreads $short_R1 $short_R2 | samtools collate -@ $threads -o ${intermediate_files}/short2long.bam - > /dev/null 2>&1
fi

if [ -f ${intermediate_files}/short_reads_barcoded_map.txt ];then
    echo "add_barcode already done"
else
    $BINDIR/add_barcode -b ${intermediate_files}/short2long.bam -m ${intermediate_files}/barcode_maps.txt -l $identity -o ${intermediate_files}/short_reads_barcoded
fi

if [ -f ${intermediate_files}/interleaved_link_reads.fastq ] || [ -f interleaved_link_reads.sorted.fastq ];then
    echo "assign barcodes already done"
else
    echo "`date` assign barcodes start"
    $BINDIR/assign_barcodes --map ${intermediate_files}/short_reads_barcoded_map.txt --fastq1 $short_R1 --fastq2 $short_R2 --output ${intermediate_files}/interleaved_link_reads.fastq
    echo "`date` assign barcodes done"
fi

if [ -f ${intermediate_files}/interleaved_link_reads.sorted.fastq ];then
    echo "sort already done"
else
    echo "sort linked reads by barcode"
    awk '{
        if (NR % 8 == 0) {
            printf("%s\n", $0);
        }else if(NR % 8==1) {
            printf("%s BX:Z:XXXXXXXXXXXXXXX\t", $0);
        }else { printf("%s\t",$0);}
    }' ${intermediate_files}/interleaved_link_reads.fastq |  LANG=C sort -k 2.1,2.21 -k 1,1  |sed 's/ BX:Z:XXXXXXXXXXXXXXX//g '| tr "\t" "\n" > ${intermediate_files}/interleaved_link_reads.sorted.fastq
    echo "`date` sort done!"
fi

if [ ! -d logs ];then
    mkdir logs
fi

if [ ! -f $athena_lc ] || [ ! -f $athena_out ];then
    if [ -d athena_out ];then rm -r athena_out; fi
    bash $CURRENT_DIR/run_athena.sh ${intermediate_files}/interleaved_link_reads.sorted.fastq $root/example/hybrid_example/contigs/metaspades.fasta $output_dir/athena_out $threads > logs/athena.log 2>&1 &
fi

# # source $CONDA_PREFIX/bin/activate pangaea
# if [ $CONDA_DEFAULT_ENV == "pangaea" ];then
#     echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea activated"
# else
#     echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea not activated"
#     exit
# fi

activate_conda_env $pangaea_env
echo "running pangaea"
echo "python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -st 1,2,3 -o $output_dir/pangaea_out"
python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -st 1,2,3 -o $output_dir/pangaea_out &

wait

if [ ! -f $root/example/hybrid_example/contigs/$type.fasta ];then
    echo "contigs $root/example/hybrid_example/contigs/$type.fasta not found"
    # TODO 
    echo "Please run the hybrid assembly first"
    exit 1
fi

if [ ! -f $athena_lc ];then
    echo "contigs $athena_lc not found, check if you run athena first"
    exit 1
fi

if [ ! -f $athena_out ];then
    echo "contigs $athena_out not found, check if you run athena first"
    exit 1
fi

python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -sp $root/example/hybrid_example/contigs/$type.fasta -lc $athena_lc -at $athena_out -st 4 -o $output_dir/pangaea_out


echo "`date` all done!"
