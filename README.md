# Pangaea
Pangaea is a linked-read assembler for linked-reads with high barcode specificity, using the variational autoencoder to bin linked-reads and multi-thresholding reassembly to assemble linked-reads.

## Installation
git clone git@github.com:ericcombiolab/Pangaea.git

### Dependencies
Pangaea depends on [numpy](https://numpy.org/install/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [sklearn](https://scikit-learn.org/stable/install.html), [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), [pysam](https://pysam.readthedocs.io/en/latest/installation.html), [torch](https://pytorch.org/get-started/locally/), [rph_kmeans](https://github.com/tinglabs/rph_kmeans), [pigz](https://zlib.net/pigz/), [bwa](https://github.com/lh3/bwa), [samtools](https://github.com/samtools/samtools), [seqtk](https://github.com/lh3/seqtk), [megahit](https://github.com/voutcn/megahit), [spades(>=v3.15.3)](https://github.com/ablab/spades), [flye](https://github.com/fenderglass/Flye), [quickmerge](https://github.com/mahulchak/quickmerge), and [jgi_summarize_bam_contig_depths](https://bitbucket.org/berkeleylab/metabat/src/master/).

### Installation using conda 
```
conda env create -f environment.yaml
conda activate pangaea

# Optional: Install Athena (based on python2), run  
conda env creat -f athena_enviroment.ymal
conda activate athena
```
Note that if you want to run Athena/Pangaea, you need to change conda enviroments accrodingly. 

### Compile cpp utils of Pangaea
```
cd Pangaea/cpptools && make && cd -
```

## Preprocessing of linked-reads
Run metaspades to obtain error-corrected reads:
```
metaspades.py --12 /path/to/reads -o /path/to/metaspades/out
```
The reads1 (input to Pangaea ```-1```), reads2 (input to Pangaea ```-2```) are in ```/path/to/metaspades/out/corrected/``` and spades contigs (input to Pangaea ```-sp```) are at ```/path/to/metaspades/out/contigs.fasta```.

Run athena to obtain local assembly contigs and athena contigs. The local assembly contigs (input to Pangaea ```-lc```) and athena contigs (input to Pangaea ```-at```) are ```/path/to/athena/out/results/olc/flye-input-contigs.fa``` and ```/path/to/athena/out/results/olc/athena.asm.fa```.
## Running Pangaea
```
usage: pangaea.py [-h] [-1 READS1] [-2 READS2] [-lreads LONG_READS]
                  [-lrtype {pacbio_raw,nanopore_raw,pacbio_corrected,nanopore_corrected}]
                  [-i INTERLEAVED_READS] -o OUTPUT [-l MIN_LENGTH] [-k KMER]
                  [-tnf_k TNF_KMER] [-s WINDOW_SIZE] [-v VECTOR_SIZE] [-r LR]
                  [-w WEIGHT_DECAY] [-e EPOCHS] [-b BATCH_SIZE] [-d DROPOUT]
                  [-p PATIENCE] [-wa WEIGHT_ALPHA] [-wk WEIGHT_KL]
                  [-ld LATENT_DIM] -c CLUSTERS [-t THREADS] [-g USE_CUDA]
                  [-n NUM_GPUS] [-sp SPADES] [-lc LOCAL_ASSEMBLY] [-at ATHENA]
                  [-lt LOW_ABD_CUT] [-md MODEL] [-wx WEIGHT_AUXILIARY]
                  [-ls LOSS_TYPE] [-cf CONFIDENCE] [-st STEPS]

optional arguments:
  -h, --help            show this help message and exit
  -1 READS1, --reads1 READS1
                        path to reads1 file (linked-reads)
  -2 READS2, --reads2 READS2
                        path to reads2 file (linked-reads)
  -lreads LONG_READS, --long_reads LONG_READS
                        path to reads file (long-reads)
  -lrtype {pacbio_raw,nanopore_raw,pacbio_corrected,nanopore_corrected}, --long_reads_type {pacbio_raw,nanopore_raw,pacbio_corrected,nanopore_corrected}
                        long reads type (default pacbio_raw)
  -i INTERLEAVED_READS, --interleaved_reads INTERLEAVED_READS
                        path to reads file (long-reads)
  -o OUTPUT, --output OUTPUT
                        output directory
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        min barcode length (default 2000)
  -k KMER, --kmer KMER  kmer for abundance (default 15)
  -tnf_k TNF_KMER, --tnf_kmer TNF_KMER
                        kmer for TNF (default 4, long reads should use 3)
  -s WINDOW_SIZE, --window_size WINDOW_SIZE
                        window size for abundance (default 10)
  -v VECTOR_SIZE, --vector_size VECTOR_SIZE
                        vector size for abundance (default 400)
  -r LR, --lr LR        learning rate (default 0.005)
  -w WEIGHT_DECAY, --weight_decay WEIGHT_DECAY
                        weight decay (default 0.0001)
  -e EPOCHS, --epochs EPOCHS
                        number of epochs (default 100)
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        batch size (defult 2048)
  -d DROPOUT, --dropout DROPOUT
                        dropout (default 0.2)
  -p PATIENCE, --patience PATIENCE
                        early stop patience (default 20)
  -wa WEIGHT_ALPHA, --weight_alpha WEIGHT_ALPHA
                        training weight for abundance and tnf (default 0.1)
  -wk WEIGHT_KL, --weight_kl WEIGHT_KL
                        training weight for KL (default 0.015)
  -ld LATENT_DIM, --latent_dim LATENT_DIM
                        latent dimension (default 32)
  -c CLUSTERS, --clusters CLUSTERS
                        number of clusters
  -t THREADS, --threads THREADS
                        number of threads (default 100)
  -g USE_CUDA, --use_cuda USE_CUDA
                        use cuda (default False)
  -n NUM_GPUS, --num_gpus NUM_GPUS
                        use gpu in parallel (if use cuda)
  -sp SPADES, --spades SPADES
                        path to original contigs
  -lc LOCAL_ASSEMBLY, --local_assembly LOCAL_ASSEMBLY
                        path to local assembly contigs
  -at ATHENA, --athena ATHENA
                        path to athena contigs
  -lt LOW_ABD_CUT, --low_abd_cut LOW_ABD_CUT
                        coverage for low abundance contigs
  -md MODEL, --model MODEL
                        model ( vae)
  -wx WEIGHT_AUXILIARY, --weight_auxiliary WEIGHT_AUXILIARY
                        training weight for auxiliary (default 0.1)
  -ls LOSS_TYPE, --loss_type LOSS_TYPE
                        reconstruction loss type (default ce)
  -cf CONFIDENCE, --confidence CONFIDENCE
                        clustering confidence
  -st STEPS, --steps STEPS
                        steps to run (default 1:feature extraction, 2:vae
                        trainning, 3:clutsering, 4:sub-assembly and final
                        assembly)
```

## Example of running Pangaea on linked reads
Under folder ```example/pangaea_out```, pangaea was run with
```
cd example
nohup python ../pangaea.py -1 reads1.fq.gz -2 reads2.fq.gz -sp contigs.fa -lc flye-input-contigs.fa -at athena.asm.fa -c 5 -o pangaea > pangaea.log 2>&1 &
```

## Example of running Pangaea on hybrid reads (short pairend reads and longreads)
```
cd example/hybrid
./hybrid_wrapper.sh atcc_longreads_small.fastq atcc_short_R1.fastq.gz atcc_short_R2.fastq.gz 60 operams > logs/log

```
###  Optional: Substituted the metaSPAdes in step 1 and Athena in step 2 with the corresponding hybrid assemblies (contigs generated from hybridSPAdes or OPERA-MS)
```
# type: operams, hybridspades
./final_merge.sh <type>
```
The new generated final assembly will be at ```pangaea_out/4.assembly/quickmerge_<type>/merged_out.fasta```
