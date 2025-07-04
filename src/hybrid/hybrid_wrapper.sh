#!/bin/bash
set -o pipefail
set -e
echo "`date` starting"
CURRENT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$CURRENT_DIR/../../
BINDIR=$root/src/bin
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

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

activate_conda_env $pangaea_env

# the file path of this script
cd src/cpptools
if [ -d "build" ]; then
  rm -rf build
fi
mkdir build; cd build;cmake ..;make
cd ../../../


while getopts "l:r:R:i:t:a:A:o:p" opt; do
    case $opt in
        l) longreads="$OPTARG" ;; # Optional. Long reads file, if not givenm we assume the short reads are linked reads with barcodes
        r) short_R1="$OPTARG" ;; # Mandatory. Short reads R1
        R) short_R2="$OPTARG" ;; # Mandatory. Short reads R2
        i) identity="$OPTARG" ;; # Optional. Identity for barcode assignment, default is 60
        t) type="$OPTARG" ;; # Optional. Hybrid assembly type, e.g., metaspades, hybridspades, metaplatanus, default is hybridspades
        a) athena_lc="$OPTARG" ;; # Optional. Local assembly result of athena-meta, default is $output_dir/athena_out/results/olc/flye-input-contigs.fa
        A) athena_out="$OPTARG" ;; # Optional. Result file of athena-meta, default is $output_dir/athena_out/results/olc/athena.asm.fa
        o) output_dir="$OPTARG" ;; # Optional. Output directory, default is current directory with 'default_hybrid_out'
        p) longreads_type="$OPTARG" ;; # Optional. Long reads type, e.g., pacbio or nanopore, default is pacbio
        *) 
            echo "Usage: $0 -l <longreads> -r <short_R1> -R <short_R2> [-i <identity>] [-t <hybrid_type>] [-a <athena_lc>] [-A <athena_out>] [-o <output_dir>] [-p <longreads_type>]"
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
if  [ -z "$short_R1" ] || [ -z "$short_R2" ]; then
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
    athena_lc="$output_dir/athena_out/results/olc/flye-input-contigs.fa"
    echo "Warning: Not specify Athena input file, using default path:$athena_lc"
fi
if [ -z "$athena_out" ]; then
    athena_out="$output_dir/athena_out/results/olc/athena.asm.fa"
    echo "Warning: Not specify Athena output file, using default path:$athena_out"
fi

if [ -z "$longreads_type" ]; then
    echo "Warning: Not specify long reads type, using default value"
    longreads_type="pacbio"
fi


if [ -z "$longreads" ]; then
    echo "Warning: Not specify long reads, assuming short reads are linked reads with barcodes"
else
    echo "Long reads file: $longreads with type $longreads_type"
fi

barcode_list=$CURRENT_DIR/barcodelist.txt
threads=50

# Generate intermediate files directory
intermediate_files=$output_dir'/intermediate_files'

echo "athena_lc: $athena_lc"
echo "athena_out: $athena_out"


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

if [ ! -f $short_R1 ]; then
    echo "short_R1 $short_R1 not found"
    exit 1
fi

if [ ! -f $short_R2 ]; then
    echo "short_R2 $short_R2 not found"
    exit 1
fi

if [ ! -z $longreads ] && [ -f $longreads ]; then

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

else
    #todo: transform 2 paired-end reads to interleaved linked reads
    echo "Warning: No long reads provided, assuming short reads are linked reads with barcodes"
    echo "Warning: If you want to use long reads, please provide the long reads file with -l option"
    if [ ! -d ${intermediate_files} ];then
        mkdir ${intermediate_files}
    fi
    if [ -f ${intermediate_files}/interleaved_link_reads.fastq ];then
        echo "interleaved_link_reads.fastq already exists, skipping interleaving"
    else
        echo "interleave short reads"
        zcat $short_R1 $short_R2 | perl -pe "s/\@(.+) (.+)/\@\$2 \$1/" > TMP.fq && fastq-sort --id TMP.fq | perl -pe "s/\@(.+) (.+)/\@\$2 \$1/" > ${intermediate_files}/interleaved_link_reads.fastq && rm -f TMP.fq
        echo "`date` interleave done!"
    fi
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


if [ ! -d $output_dir/logs ];then
    mkdir $output_dir/logs
fi

if [ -f $output_dir/metaspades_out/contigs.fasta ];then
    echo "metaspades_out/contigs.fasta already exists, skipping metaspades assembly"
else
    echo "running metaspades assembly"
    echo "metaspades.py --checkpoints last -1 $short_R1 -2 $short_R2 -t 50 -o $output_dir/metaspades_out"
    metaspades.py --checkpoints last -1 $short_R1 -2 $short_R2 -t 50 -o $output_dir/metaspades_out > $output_dir/logs/metaspades.log 2>&1 
fi

 fisrt_assembly_pid=0  # Set to 0 to indicate no assembly process is running
# if $type is specified, run hybrid assembly
if [ $type == "metaspades" ];then
    # check if the metaspades_out directory exists, if not exit
    if [ -d $output_dir/metaspades_out ];then
        echo "metaspades_out directory already exists, skipping metaspades assembly"
    else
        echo "error: metaspades_out directory does not exist, please run metaspades assembly first"
        exit 1
    fi
elif [ $type == "hybridspades" ];then
    if [ -f $output_dir/hybridspades_out/contigs.fasta ];then
        echo "hybridspades_out/contigs.fasta already exists, skipping hybridspades assembly"
    else
        echo "running hybridspades assembly"
        if [ $longreads_type == "nanopore" ];then
            echo "using nanopore long reads"
            echo "metaspades.py --checkpoints last -1 $short_R1 -2 $short_R2 --nanopore $longreads -t 50 -o $output_dir/hybridspades_out"
            metaspades.py --checkpoints last -1 $short_R1 -2 $short_R2 --nanopore $longreads -t 50 -o $output_dir/hybridspades_out > $output_dir/logs/hybridspades.log 2>&1 &
            fisrt_assembly_pid=$!  # Capture the PID of the hybridspades process
            echo "running hybrid assembly $type with PID $fisrt_assembly_pid"

        else
            echo "using pacbio long reads"
            metaspades.py --checkpoints last -1 $short_R1 -2 $short_R2 --pacbio $longreads -t 50 -o $output_dir/hybridspades_out > $output_dir/logs/hybridspades.log 2>&1 &
            fisrt_assembly_pid=$!  # Capture the PID of the hybridspades process
            echo "running hybrid assembly $type with PID $fisrt_assembly_pid"
        fi
    fi
elif [ $type == "metaplatanus" ];then
    echo "running metaplatanus assembly"
    mkdir -p hybrid_out_metaplatanus/metaplatanus_out/
    echo "/usr/bin/time -v metaplatanus -IP1 $short_R1 $short_R2 -p $longreads -t $threads -o $output_dir/metaplatanus_out/ -m 500 -t 64 "
    /usr/bin/time -v metaplatanus -IP1 $short_R1 $short_R2 -p $longreads -t $threads -o $output_dir/metaplatanus_out/ -m 500 -t 64 > $output_dir/logs/metaplatanus.log 2>&1 &
    fisrt_assembly_pid=$!  # Capture the PID of the metaplatanus process
else
    echo "Unknown type: $type. Please specify 'metaspades', 'hybridspades', or 'metaplatanus'. If you want to run operams, please install operams by your own and run final_merge.sh after this."
    exit 1
fi

athena_pid=0
if [ ! -f $athena_lc ] || [ ! -f $athena_out ];then
    if [ -d $output_dir/athena_out ];then rm -r $output_dir/athena_out; fi
    echo "bash $CURRENT_DIR/run_athena.sh ${intermediate_files}/interleaved_link_reads.sorted.fastq $output_dir/metaspades_out/contigs.fasta $output_dir/athena_out $threads"
    bash $CURRENT_DIR/run_athena.sh ${intermediate_files}/interleaved_link_reads.sorted.fastq $output_dir/metaspades_out/contigs.fasta $output_dir/athena_out $threads > $output_dir/logs/athena.log 2>&1 &
    athena_pid=$!  # Capture the PID of athena process
    echo "`date` running athena with PID $athena_pid"
fi



echo "running pangaea"
pangaea_pid=0
echo "python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -st 1,2,3 -o $output_dir/pangaea_out"
python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -st 1,2,3 -o $output_dir/pangaea_out &
pangaea_pid=$!  # Capture the PID of the pangaea process
echo "`date` running pangaea with PID $pangaea_pid"

if [ $athena_pid -eq 0 ]; then
    echo "Athena is not running, skipping Athena step"
else
    echo "`date` waiting for athena to finish"
    wait $athena_pid
    athena_status=$?
    if [ $athena_status -ne 0 ]; then
        echo "Athena failed with status $athena_status"
        pkill -P $$
        exit $athena_status
    else
        echo "`date` athena done"
    fi
fi

if [ $fisrt_assembly_pid -eq 0 ]; then
    echo "No assembly process is running, skipping assembly step"
else
    echo "`date` waiting for first assembly to finish"
    wait $fisrt_assembly_pid
    first_assembly_status=$?
    if [ $first_assembly_status -ne 0 ]; then
        echo "Assembly failed with status $first_assembly_status"
        pkill -P $$
        exit $first_assembly_status
    elif [ $type == "metaspades" ];then
        echo "`date` metaspades assembly done"
    elif [ $type == "hybridspades" ];then
        echo "`date` hybridspades assembly done"
    elif [ $type == "metaplatanus" ];then
        echo "`date` metaplatanus assembly done"
    fi
fi

if [ $pangaea_pid -eq 0 ]; then
    echo "No Pangaea process is running, skipping Pangaea step"
else
    echo "`date` waiting for Pangaea to finish"
    wait $pangaea_pid
    pangaea_status=$?
    if [ $pangaea_status -ne 0 ]; then
        echo "Pangaea failed with status $pangaea_status"
        pkill -P $$
        exit $pangaea_status
    else
        echo "`date` pangaea done"
    fi
fi

# if type is specified, check if the output contig file exists
if [ $type == "metaspades" ];then
    type_contigs="$output_dir/metaspades_out/contigs.fasta"
elif [ $type == "hybridspades" ];then
    type_contigs="$output_dir/hybridspades_out/contigs.fasta"
elif [ $type == "metaplatanus" ];then
    type_contigs="$output_dir/metaplatanus_out/contigs.fasta"
else
    echo "Unknown type: $type. Please specify 'metaspades', 'hybridspades', or 'metaplatanus'. If you want to run operams, please install operams by your own and run final_merge.sh after this."
    exit 1
fi
# verify if the file exists
if [ ! -f $type_contigs ];then
    echo "contigs $type_contigs not found, check if you run the $type assembly first"
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

python $root/src/pangaea.py -i ${intermediate_files}/interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -sp $type_contigs -lc $athena_lc -at $athena_out -st 4 -o $output_dir/pangaea_out


echo "`date` all done!"
