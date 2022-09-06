# Pangaea
Pangaea is a linked-read assembler for linked-reads with high specificity, using the variational autoencoder to binning linked-reads and multi-thresholding reassembly to assemble linked-reads.
 ## Installation
Pangaea depends on numpy, pandas, sklearn, snakemake, pysam, torch, rph_kmeans, bwa, samtools, seqtk, megahit, spades, flye, quickmerge, jgi_summarize_bam_contig_depths.

To install Pangaea, use the following script:
```
git clone git@github.com:ericcombiolab/Pangaea.git
cd Pangaea/cpptools && make
```

## Running
```
pangaea.py [-h] -1 READS1 -2 READS2 -o OUTPUT [-l MIN_LENGTH] [-k KMER]
                  [-s WINDOW_SIZE] [-v VECTOR_SIZE] [-r LR] [-w WEIGHT_DECAY]
                  [-e EPOCHS] [-b BATCH_SIZE] [-d DROPOUT] [-p PATIENCE]
                  [-wa WEIGHT_ALPHA] [-wk WEIGHT_KL] [-ld LATENT_DIM] -c
                  CLUSTERS [-t THREADS] [-g USE_CUDA] [-n NUM_GPUS] -sp SPADES
                  -lc LOCAL_ASSEMBLY -at ATHENA [-lt LOW_ABD_CUT]

optional arguments:
  -h, --help            show this help message and exit
  -1 READS1, --reads1 READS1
                        path to reads1 file (linked-reads)
  -2 READS2, --reads2 READS2
                        path to reads2 file (linked-reads)
  -o OUTPUT, --output OUTPUT
                        output directory
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        min barcode length (default 2000)
  -k KMER, --kmer KMER  kmer for abundance (default 15)
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
```