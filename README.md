# Pangaea
Pangaea is designed to assemble short-reads with high specificity physical (linked-reads) or virtual barcodes (long-reads+short-reads). It includes (1) short-reads binning using variational autoencoder (2)multi-thresholding reassembly and (3) ensemble assembly.

## Installation
An all-in-one installation script. This may take about 20 minutes ~ 1 hours.
```
git clone https://github.com/ericcombiolab/Pangaea.git
cd Pangaea
./build.sh
```
### Dependencies
Pangaea depends on [numpy >= 1.23.5](https://numpy.org/install/), [pandas >= 1.5.3](https://pandas.pydata.org/docs/getting_started/install.html), [sklearn >= 1.2.2](https://scikit-learn.org/stable/install.html), [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), [pysam >= 0.20.0](https://pysam.readthedocs.io/en/latest/installation.html), [torch 1.10.0](https://pytorch.org/get-started/locally/), [rph_kmeans](https://github.com/tinglabs/rph_kmeans), [pigz 2.4](https://zlib.net/pigz/), [bwa >= 0.7.17](https://github.com/lh3/bwa), [samtools 1.9](https://github.com/samtools/samtools), [seqtk](https://github.com/lh3/seqtk), [megahit v1.2.9](https://github.com/voutcn/megahit), [spades(>=v3.15.3)](https://github.com/ablab/spades), [flye 2.8-b1674](https://github.com/fenderglass/Flye), [quickmerge](https://github.com/mahulchak/quickmerge), [Jellyfish 2.3.0](https://github.com/gmarcais/Jellyfish)and [jgi_summarize_bam_contig_depths](https://bitbucket.org/berkeleylab/metabat/src/master/).

The examples and experiments are mostly running on Linux servers with the following hardware specifications:
- Dell PowerEdge R6525
- CPU: Dual 64-core AMD EPYC 7742 2.25GHz 256MB L3 cache
- Memory: 1T


### Build conda environment（not needed if you have run ```./build.sh``` ）
```
conda env create -f environment.yaml
conda activate pangaea

# Optional: Install Athena (based on python2), run  
conda env creat -f athena_enviroment.ymal
conda activate athena
```
Note that if you want to run Athena/Pangaea, you need to change conda enviroments accrodingly. 

### Compile cpp utils of Pangaea（not needed if you have run ```./build.sh``` ）
```
cd Pangaea/cpptools && make && cd -
```

## Preprocessing of linked-reads
Run metaspades to obtain error-corrected reads and seed contigs for Athena and Pangaea:
```
metaspades.py --12 /path/to/reads -o /path/to/metaspades/out
```
The reads1 (input to Pangaea ```-1```), reads2 (input to Pangaea ```-2```) are in ```/path/to/metaspades/out/corrected/``` and spades contigs (input to Pangaea ```-sp```) are at ```/path/to/metaspades/out/contigs.fasta```.

Run athena to obtain local assembly contigs and athena contigs. The local assembly contigs (input to Pangaea ```-lc```) and athena contigs (input to Pangaea ```-at```) are ```/path/to/athena/out/results/olc/flye-input-contigs.fa``` and ```/path/to/athena/out/results/olc/athena.asm.fa```.

## Running Pangaea
```
usage: pangaea.py [-h] [-1 READS1] [-2 READS2] [-i INTERLEAVED_READS] -o
                  OUTPUT [-l MIN_LENGTH] [-k KMER] [-tnf_k TNF_KMER]
                  [-s WINDOW_SIZE] [-v VECTOR_SIZE] [-r LR] [-w WEIGHT_DECAY]
                  [-e EPOCHS] [-b BATCH_SIZE] [-d DROPOUT] [-p PATIENCE]
                  [-wa WEIGHT_ALPHA] [-wk WEIGHT_KL] [-ld LATENT_DIM] -c
                  CLUSTERS [-t THREADS] [-g USE_CUDA] [-n NUM_GPUS]
                  [-sp SPADES] [-lc LOCAL_ASSEMBLY] [-at ATHENA]
                  [-lt LOW_ABD_CUT] [-la LOW_ASSEMBLER] [-md MODEL]
                  [-ls LOSS_TYPE] [-st STEPS]

optional arguments:
  -h, --help            show this help message and exit
  -1 READS1, --reads1 READS1
                        path to reads1 file (linked-reads)
  -2 READS2, --reads2 READS2
                        path to reads2 file (linked-reads)
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
                        batch size (default 2048)
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
  -la LOW_ASSEMBLER, --low_assembler LOW_ASSEMBLER
                        local assembly method (spades or megahit)
  -md MODEL, --model MODEL
                        model ( vae)
  -ls LOSS_TYPE, --loss_type LOSS_TYPE
                        reconstruction loss type (default ce)
  -st STEPS, --steps STEPS
                        steps to run (default 1:feature extraction, 2:vae
                        training, 3:clutsering, 4:sub-assembly and final
                        assembly)
```

## Example of running Pangaea on short-reads with physical barcodes (linked reads)
```
conda activate pangaea
cd example
nohup python ../pangaea.py -1 reads1.fq.gz -2 reads2.fq.gz -sp contigs.fa -lc flye-input-contigs.fa -at athena.asm.fa -c 5 -o pangaea > pangaea.log 2>&1 &
```
This may take about 1~2 hours. 

## Example of running Pangaea on short-reads with virtual barcodes  (long-reads and short-reads)
```
cd example/hybrid
nohup ./hybrid_wrapper.sh atcc_longreads_small.fastq.gz atcc_short_R1.fastq.gz atcc_short_R2.fastq.gz 60 operams > hybrid.log 2>&1 &
```
This may take about 1~2 hours. 


###  Optional: Substitute the metaSPAdes in step 1 and Athena in step 2 with the corresponding hybrid assemblers (such as hybridSPAdes or OPERA-MS)
```
# type: operams, hybridspades
./final_merge.sh <type>
```
The new generated final assembly will be at ```pangaea_out/4.assembly/quickmerge_<type>/merged_out.fasta```
This may take about 10~30 minutes.