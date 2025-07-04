# ncbin output file
set -e
CURRENT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$CURRENT_DIR/../../
BINDIR=$root/src/bin
SCRIPTDIR=$root/src/scripts
# get the work directory by input parameters
if [ $# -lt 2 ];then
    echo "Usage: $0 <path_to_pangaea_root> <type>"
    echo "type: metaspades, hybridspades, metaplatanus"
    exit 1
fi
path=`realpath $1`
type=$2
if [ ! -d $path ];then
    echo "Path $path does not exist!"
    exit 1
else
    pangaea_out=$path/pangaea_out
    athena_local=$path/athena_out/results/olc/flye-input-contigs.fa
    athena=$path/athena_out/results/olc/athena.asm.fa
    cluster_dir=$path/pangaea_out/3.clustering
    assembly_dir=$path/pangaea_out/4.assembly
    pangaea_local=$path/pangaea_out/4.assembly/contigs.low_abd.binning.local.fa
    #assert all directories and files exist, use one line
    for dir in $pangaea_out $athena_local $athena $cluster_dir $assembly_dir; do
        if [ ! -d $dir ] && [ ! -f $dir ]; then
            echo "Directory or file $dir does not exist!"
            exit 1
        else
            echo "Directory or file $dir exists."
        fi
    done
fi


spades=$path/metaspades_out/contigs.fasta
hybridspades=$path/hybridspades_out/contigs.fasta
metaplatanus=$path/metaplatanus_out/_result/out_final.fa

pangaea_local=$assembly_dir/contigs.low_abd.binning.local.fa

if [ ! -d $final_asm ];then
    mkdir $final_asm
fi

if [ -f ${pangaea_local}.fai ];then
    rm ${pangaea_local}.fai
fi
    echo "`date` - Concatenating contigs"
    cat $assembly_dir/*.spades/contigs.fasta $cluster_dir/contigs.megahit.fa $athena_local > $pangaea_local
    $BINDIR/parse_header $pangaea_local contig_ > $assembly_dir/contigs.low_abd.binning.local.renamed.fa
    mv $assembly_dir/contigs.low_abd.binning.local.renamed.fa $pangaea_local
    echo "`date` - Concatenating done"

# types=( "metaplatanus" "metaspades" "hybridspades" )



if [ $type == "metaspades" ];then
    contig=$spades
elif [ $type == "hybridspades" ];then
    contig=$hybridspades
elif [ $type == "metaplatanus" ];then
    contig=$metaplatanus
else
    echo "Unknown type: $type. Please specify 'metaspades', 'hybridspades', or 'metaplatanus'. If you want to run operams, please install operams by your own"
    exit 1
fi
echo $type
echo $contig
if [ ! -f $contig ];then
    echo "contig:  $contig not exist!"
    continue
fi
if [ ! -f $assembly_dir/olc_$type/final.asm.fa ];then
    echo "`date` Performing olc based on $type"
    echo "$SCRIPTDIR/merge_olc.py $contig $pangaea_local $assembly_dir/olc_$type"
    $SCRIPTDIR/merge_olc.py $contig $pangaea_local $assembly_dir/olc_$type > /dev/null
    echo "`date` olc based on $type done"
else
    echo "`date` olc based on $type already done"
fi


    

if [ ! -f $assembly_dir/quickmerge_$type/merged_out.fasta ];then
    echo "`date` Performing final quickmerge based on $type"
    if [ ! -d $assembly_dir/quickmerge_$type ]; then
        mkdir $assembly_dir/quickmerge_$type
    fi
    athena=`realpath $athena`
    cd $assembly_dir/quickmerge_$type
    # quickmerge
    merge_wrapper.py $assembly_dir/olc_$type/final.asm.fa $athena > /dev/null
    $BINDIR/parse_header merged_out.fasta contig_ > merged_out.renamed.fasta
    mv merged_out.renamed.fasta merged_out.fasta
    echo "`date` quickmerge based on $type done"
else
    echo "`date` quickmerge based on $type already done"
fi


wait
echo "`date` all done"
