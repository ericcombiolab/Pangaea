# ncbin output file
set -e
BINDIR=`realpath "../../bin"`
root=$PWD
cluster_dir=$root/pangaea_out/3.clustering
assembly_dir=$root/pangaea_out/4.assembly
athena_local=$root/contigs/flye-input-contigs.fa
athena=$root/contigs/athena.asm.fa
spades=$root/contigs/metaspades.fasta
hybridspades=$root/contigs/hybridspades.fasta
operams=$root/contigs/operams.fasta

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

# types=( "operams" "metaspades" )
types=( $* )

for type in ${types[@]}
do
    if [ $type == "metaspades" ];then
        contig=$spades
    elif [ $type == "hybridspades" ];then
        contig=$hybridspades
    elif [ $type == "operams" ];then
        contig=$operams
    fi
    echo $type
    echo $contig
    if [ ! -f $contig ];then
        echo "contig:  $contig not exist!"
        continue
    fi
    if [ ! -f $assembly_dir/olc_$type/final.asm.fa ];then
        echo "`date` Performing olc based on $type"
        $BINDIR/merge_olc.py $contig $pangaea_local $assembly_dir/olc_$type > /dev/null
        echo "`date` olc based on $type done"
    else
        echo "`date` olc based on $type already done"
    fi
done

    
for type in ${types[@]}
do
    # run the block in parallel
    {
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
    } &
done

wait
echo "`date` all done"
