# Pangaea
Pangaea is designed for linked-read assembly or hybrid asembly using short- and long-reads. It includes (1) short-reads binning using variational autoencoder (2) multi-thresholding reassembly and (3) ensemble of different subassemblies.

## Installation

We provide an all-in-one installation script (`build.sh`) that sets up all dependencies, Conda environments, and builds the required tools. The process typically takes 20 minutes to 1 hour, depending on your system and internet speed.

### Quick Start
```bash
git clone https://github.com/ericcombiolab/Pangaea.git
cd Pangaea
./build.sh
```

### Step-by-step Installation Details

1. **Set up Conda environments**
    - The script uses `mamba` (a faster drop-in replacement for conda) to create environments.
    - Install the main `pangaea` environment:
      ```bash
      # if mamba is not installed, install mamba first
      conda install conda-forge::mamba
      # this will create conda env pangaea
      mamba env create -f environment.yaml
      ```
    - Install the `athena-meta` environment for Athena (required for hybrid assembly):
      ```bash
      # this will create conda env athena-meta
      mamba create -n athena-meta bioconda::athena_meta
      ```

2. **Install PyTorch**
    - Please follow [PyTorch's official instructions](https://pytorch.org/get-started/locally/) for your system and hardware. For CPU-only installation, you can use:
      ```bash
      pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
      ```

3. **Install rph_kmeans**
    - Install the `rph_kmeans` package from source:
      ```bash
      conda activate pangaea
      cd third_parties/rph_kmeans
      python setup.py install
      cd ../../
      ```

4. **Build C++ tools**
    - Compile the C++ binaries required by Pangaea:
      ```bash
      cd src/cpptools
      rm -rf build
      mkdir build && cd build
      cmake ..
      make
      cd ../../../
      ```

5. **Set up MetaPhlan4 database (optional)**
    - After running `build.sh`, you can choose to download the MetaPhlan4 database:
      ```bash
      # Download the database to ./metaphlan4_DB
      ./build_db.sh
      # Or specify a custom directory (ensure >25GB free space):
      ./build_db.sh -d [your_preferred_directory]
      ```
    - MetaPhlan4 is used for automatic cluster number selection. If you prefer to specify the cluster number manually (recommended, as results are not very sensitive to this parameter), you can skip the MetaPhlan4 database step and use the `-c [number]` option when running Pangaea.
    - The cluster number is a trade-off: a larger value produces more and smaller read bins (lower complexity per bin), while a smaller value keeps more reads from the same microbe together. We suggest `-c 30` for most datasets. For very complex datasets, try `-c 35` or `-c 40`. For mock or simulation data, we suggest using a value lower than the real species number or a small value (e.g., `-c 10`).

---


## Running Pangaea
```
Usage: ./src/run_pangaea [OPTIONS]
Required arguments:
  -s, --short_type <string>       Short reads type: short, stlfr, tellseq, 10x
                                    If 'short', hybrid assembly is performed and -l, -H, and -p are required.
                                    If 'stlfr', 'tellseq', or '10x', linked assembly is performed and -l, -H, and -p must NOT be set.
                                    For 'stlfr' or 'tellseq', original reads can be directly provided to Pangaea without preprocessing.
                                    For '10x', please ensure barcodes are in the BX:Z tag of the read headers (this can be done using 'longranger basic', followed by deinterleaving the reads).
  -r, --short_R1 <file>           Short reads R1 file
  -R, --short_R2 <file>           Short reads R2 file
  -I, --index <file>              Barcode index for Tell-Seq (required if -s is 'tellseq'; this file is provided with the reads)
  -o, --output_dir <dir>           Output directory (required)
Hybrid assembly (required if -s is 'short'):
  -l, --longreads <file>          Long reads file
  -H, --hybrid_asm <string>       Hybrid assembler: hybridspades, metaplatanus (default: hybridspades)
  -p, --longreads_type <string>   Long reads type: pacbio or nanopore
Optional arguments:
  -t, --threads <int>             Number of threads to use (default: 50; applied to all tools that support it)
  -h, --help                      Show this help message and exit
```

The assembled contigs will be in output_dir/final_asm.fa.

## Example of running linked-read assembly
```
cd example/linked_reads_example
run_pangaea -s stlfr -r reads1.fq.gz -R reads2.fq.gz -o pangaea
```

# Example of running hybrid assembly

```
cd example/hybrid_example
run_pangaea -s short -r atcc_short_R1.fastq.gz -R atcc_short_R2.fastq.gz -l atcc_longreads_small.fastq.gz -p pacbio -o pangaea
```
