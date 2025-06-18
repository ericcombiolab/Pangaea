# Pangaea
Pangaea is designed to assemble short-reads with high specificity physical (linked-reads) or virtual barcodes (long-reads+short-reads). It includes (1) short-reads binning using variational autoencoder (2)multi-thresholding reassembly and (3) ensemble assembly.

## TODO
- [ ] Add more examples and wiki.
- [ ] Conda package through Bioconda channel.

## Docker
We have built a docker image for users to directly run Pangaea (without automaticly select cluster numbers).
```
docker pull jmelody/pangaea:std
git clone https://github.com/ericcombiolab/Pangaea.git
nohup docker run  -v $PWD/Pangaea/example/:/example -u $(id -u):$(id -g)  jmelody/pangaea:std /bin/bash /app/run_test.sh -1 /example/reads1.fq.gz -2 /example/reads2.fq.gz -sp /example/contigs.fa -lc /example/flye-input-contigs.fa -at /example/athena.asm.fa -o /example/pangaea -c 5 &
```
Please remember set ```-u $(id -u):$(id -g)``` to avoid using root user. Otherwise the Pangaea/example directory will also be root permission.
If you got error when running Docker. Please: 1. check the log file ```example/pangaea/log```. You can also open a issue and share the log file with us.  
Or, 2. you can run/follow the script ```build.sh``` which will create a conda envrionment and build the necessary packages step by step. 

## Installation

An all-in-one installation script. This may take about 20 minutes ~ 1 hours. 

Note! This script will download MetaPhlan4's latest database into this directory $PWD/metaplan4_DB, which requires about 20G storage and the process may take long. If you would like to specify a dedicated path for this, please run ```./build.sh -d [path]``` instead. 

Suggestion: MetaPhlan4 is used for clustering number choosing, if you tend to specify a clustering number by your own, you can just skip this installation step and run pangaea with -c [number] option. According to our experiments, we suggest using a number around 35 as the default c for complexed metagenomic data( in which their species number is above 200 and shannon diversity is around 3.5 ). For mock or simulation data, we suggest using a smaller or equal number to the real species number if the species number is lower than 35. 
```
git clone https://github.com/ericcombiolab/Pangaea.git
cd Pangaea
./build.sh
# set up metaphlan4
./build_db.sh (this will download the database in current directory named ./metaphlan4_DB)
# or 
./build_db.sh -d [your prefered directory(please make sure the space is larger than 25G)]
```
### Dependencies
Pangaea depends on [numpy >= 1.23.5](https://numpy.org/install/), [pandas >= 1.5.3](https://pandas.pydata.org/docs/getting_started/install.html), [sklearn >= 1.2.2](https://scikit-learn.org/stable/install.html), [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), [pysam >= 0.20.0](https://pysam.readthedocs.io/en/latest/installation.html), [torch 1.10.0](https://pytorch.org/get-started/locally/), [rph_kmeans](https://github.com/tinglabs/rph_kmeans), [pigz 2.4](https://zlib.net/pigz/), [bwa >= 0.7.17](https://github.com/lh3/bwa), [samtools 1.9](https://github.com/samtools/samtools), [seqtk](https://github.com/lh3/seqtk), [megahit v1.2.9](https://github.com/voutcn/megahit), [spades(>=v3.15.3)](https://github.com/ablab/spades), [flye 2.8-b1674](https://github.com/fenderglass/Flye), [quickmerge](https://github.com/mahulchak/quickmerge), [Jellyfish 2.3.0](https://github.com/gmarcais/Jellyfish)and [jgi_summarize_bam_contig_depths](https://bitbucket.org/berkeleylab/metabat/src/master/).

Examples and experiments run on a Linux server with the following specifications:

Hardware:
- Dell PowerEdge R6525
- CPU: Dual 64-core AMD EPYC 7742 2.25GHz 256MB L3 cache
- Memory: 1T

Software:
- Oracle Linux 8.7 (64-bit)
- gcc version 8.5.0
- conda 4.13.0
- python 3.8


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
                        steps to run (default 1:feature extraction, 2:vae trainning, 3:clutsering, 4:sub-assembly and final
                        assembly)
```

## Example of running Pangaea on short-reads with physical barcodes (linked reads)
```
conda activate pangaea
cd example
nohup python ../pangaea.py -1 reads1.fq.gz -2 reads2.fq.gz -sp contigs.fa -lc flye-input-contigs.fa -at athena.asm.fa -c 5 -o pangaea > pangaea.log 2>&1 &
```
The generated final assembly will be at ```pangaea/final.asm.fa```.

This may take about 1~2 hours. 

# Pangaea Hybrid
## Example of running Pangaea on short-reads with virtual barcodes  (long-reads and short-reads)
### Installation
- Install athena
```
 mamba env create -f athena_environment.yaml 
```
- Install pangea
```
# If you haven't install the pangaea environment yet
mamba env create -f environment.yaml
```

- compile and install necessary packages
```
bash build.sh
```

- run the hybrid example, see more details on how to understand each steps 
```
# bash src/hybrid/hybrid_wrapper.sh -l <longreads> -r <short_R1> -R <short_R2> [-i <identity>] [-t <type>] [-a <athena_lc>] [-A <athena_out>] [-o <output_dir>]
nohup bash src/hybrid/hybrid_wrapper.sh -l example/hybrid_example/atcc_longreads_small.fastq.gz -r example/hybrid_example/atcc_short_R1.fastq.gz -R example/hybrid_example/atcc_short_R2.fastq.gz -i 60 -t operams -o hybrid_out > hybrid.log 2>&1 &
```
If you encountered the Glibc problem
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found
```
You can install a compatible GCC version directly within Conda:

1. **Install GCC in Your Conda Environment**:
   ```bash
   conda install -c conda-forge gcc_linux-64 gxx_linux-64 cmake
   ```

2. **Verify the Installation**:
   ```bash
   strings $CONDA_PREFIX/lib/libstdc++.so.6 | grep GLIBCXX
   ```

   Ensure `GLIBCXX_3.4.26` is listed.

3. **Link the Updated Libraries**:
   Export the `LD_LIBRARY_PATH` to point to the Conda environment's libraries:
   ```bash
   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
   ```

---

### **Test the Fix**

After updating `libstdc++`, test your script again:

```bash
python /datahome/datasets/ericteam/csjxiao/Project/Pangaea_git/Pangaea_git_revise_250618/Pangaea/src/pangaea.py
```

If the issue is fixed, the script should run without errors.

---

### **Summary of Fixes**
1. Use `strings /lib64/libstdc++.so.6 | grep GLIBCXX` to check your `GLIBCXX` versions.
2. Update GCC and `libstdc++`:
   - System-wide (if you have root access).
   - Locally (if you don’t have root access).
   - Within a Conda environment (using `gcc_linux-64` from `conda-forge`).
3. Verify the new `GLIBCXX_3.4.26` is available and retry running your script.

This should resolve the `GLIBCXX_3.4.26` error. Let me know if you encounter further issues!

The generated final assembly will be at ```pangaea_out/final.asm.fa```.

This may take about 1~2 hours. 


###  Optional: Substitute the metaSPAdes in step 1 and Athena in step 2 with the corresponding hybrid assemblers (such as hybridSPAdes or OPERA-MS)
```
# type: operams, hybridspades
./final_merge.sh <type>
```
The new generated final assembly will be at ```pangaea_out/4.assembly/quickmerge_<type>/merged_out.fasta```.

This may take about 10~30 minutes.
