echo "`date` operams start"
longreads=$1
short_R1=$2
short_R2=$3
threads=100

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

perl /home/comp/zmzhang/software/OPERA-MS-0.8.3/OPERA-MS.pl --long-read $longreads --short-read1 $short_R1 --short-read2 $short_R2 --num-processors $threads --out-dir operams_out > operams.log 2>&1 
echo "`date` operams end"