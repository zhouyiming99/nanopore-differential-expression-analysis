# nanopore-differential-expression-analysis

This repository provides a comprehensive suite of scripts for processing Oxford Nanopore Technologies (ONT) Direct RNA sequencing data. The pipeline covers the entire workflow from raw signal processing to differential expression analysis.

## Pipeline Components

### 1. Preprocessing & Alignment
*   **`basecalling.sh`**: Basecalls raw FAST5 data using **Guppy**, merges FASTQ files, and sorts successful reads.
*   **`minimap2.sh`**: Performs splice-aware alignment to the reference genome/transcriptome using **Minimap2** and generates sorted BAM files.

### 2. RNA Modification Detection
*   **`nanopsu.sh`**: A full pipeline for **Pseudouridine ($\Psi$)** detection, including alignment cleaning, feature extraction, and machine learning prediction.
*   **`m6anet.sh`**: Detects **m6A** modifications by integrating **Nanopolish** (event alignment) and **m6anet** inference.
*   **`modkit.sh`**: Batch processes BAM files to extract **PolyA tail lengths** and base modification pileups using **ONT Modkit**.

### 3. Transcript Quantification & Discovery
*   **`isoquant.sh`**: Quantifies transcript and gene-level expression using **IsoQuant**, supporting both individual sample processing and batch merging.
*   **`bambu.R`**: An R script for novel transcript discovery and quantification using **Bambu**, specifically optimized for Nanopore DRS data.

### 4. Differential Expression Analysis
*   **`deseq2.R`**: Performs gene and transcript-level differential expression analysis (DEA) using **DESeq2**. Includes data cleaning, PCA plots, Volcano plots, and heatmaps.

## Prerequisites

Ensure the following tools are installed and available in your environment:
*   **Basecalling**: Guppy / Dorado
*   **Alignment/Bioinformatics**: Minimap2, Samtools
*   **Modification Tools**: Nanopolish, m6anet, Nanopsu, Modkit
*   **Quantification**: IsoQuant, Bambu (R package)
*   **Statistics**: R (with DESeq2, ggplot2, pheatmap)

## Usage

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```

2.  **Configure Paths**:
    Most scripts contain absolute paths for data and reference files. **You must edit the variables at the top of each script** (e.g., `REF`, `INPUT_DIR`, `THREADS`) to match your server environment.

3.  **Run the analysis**:
    ```bash
    # Set execution permissions
    chmod +x *.sh

    # 1. Basecalling and Alignment
    ./basecalling.sh
    ./minimap2.sh

    # 2. Modification Analysis
    ./nanopsu.sh
    ./m6anet.sh
    ./modkit.sh

    # 3. Quantification and DEA
    ./isoquant.sh
    Rscript bambu.R
    Rscript deseq2.R
    ```

## Key Outputs
*   `alignment/prediction.csv`: Pseudouridine ($\Psi$) site probabilities.
*   `m6anet_results/`: m6A modification probability scores.
*   `*_modifications.bed`: Base modification frequencies.
*   `*_polya_lengths.txt`: Estimated PolyA tail length per read.
*   `DESeq2_results/`: Differential expression tables and visualization (PCA, Volcano plots).

## License
This project is for research purposes. Please ensure you cite the respective tool authors when using these scripts for publication.
