#!/bin/bash
set -o pipefail
set -e
echo "`date` starting"
cd ../../cpptools/;make ;cd -
longreads=$1
short_R1=$2
short_R2=$3
identity=$4
type=$5
barcode_list=barcodelist.txt
threads=100
BINDIR=`pwd`/../../bin
root=`pwd`/../../
intermediate_files='intermediate_files'
athena_lc=`pwd`'./athena_out/results/olc/flye-input-contigs.fa'
athena_out=`pwd`'./athena_out/results/olc/athena.asm.fa'
if [ -f contigs/athena.asm.fa ] && [ -f contigs/flye-input-contigs.fa ];then
    athena_lc=`pwd`'/contigs/flye-input-contigs.fa'
    athena_out=`pwd`'/contigs/athena.asm.fa'
fi
echo "athena_lc: $athena_lc"
echo "athena_out: $athena_out"

# print command
echo "$0 $*"

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

if [ -f interleaved_link_reads.sorted.fastq ];then
    echo "sort already done"
else
    echo "sort linked reads by barcode"
    awk '{                      
        if (NR % 8 == 0) {
            printf("%s\n", $0);
        }else if(NR % 8==1) {
            printf("%s BX:Z:XXXXXXXXXXXXXXX\t", $0);
        }else { printf("%s\t",$0);}
    }' ${intermediate_files}/interleaved_link_reads.fastq |  LANG=C sort -k 2.1,2.21 -k 1,1  |sed 's/ BX:Z:XXXXXXXXXXXXXXX//g '| tr "\t" "\n" > interleaved_link_reads.sorted.fastq
    echo "`date` sort done!"
fi

if [ ! -d logs ];then
    mkdir logs
fi

if [ ! -f $athena_lc ] || [ ! -f $athena_out ];then
    if [ -d athena_out ];then rm -r athena_out; fi
    ./run_athena.sh interleaved_link_reads.sorted.fastq metaspades_out/contigs.fasta athena_out $threads > logs/athena.log 2>&1 &
fi


python ../../pangaea.py -i interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -st 1,2,3 -o pangaea_out &

wait

if [ ! -f ./contigs/$type.fasta ];then
    echo "contigs ./contigs/$type.fasta not found"
    exit 1
fi

if [ ! -f $athena_lc ];then
    echo "contigs $athena_lc not found"
    exit 1
fi

if [ ! -f $athena_out ];then
    echo "contigs $athena_out not found"
    exit 1
fi

python $root/pangaea.py -i interleaved_link_reads.sorted.fastq  -md vae -tnf_k 4 -c 15 -lt 10,30 -sp ./contigs/$type.fasta -lc $athena_lc -at $athena_out -st 4 -o pangaea_out


echo "`date` all done!"
