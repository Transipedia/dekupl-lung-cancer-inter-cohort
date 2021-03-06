# DEkupl analysis of cancer datasets: a replicability study in lung cancer.

Exhaustive capture of biological variation in RNA-seq data through k-mer decomposition (article: https://doi.org/10.1186/s13059-017-1372-2, pre-print: http://biorxiv.org/content/early/2017/06/02/122937).

DE-kupl is a computational protocol that aims to capture all k-mer variation in an input set of RNA-seq libraries. To verify the replicability of DE-kupl, we developed this pipeline to compare the consistency of events between different cohorts. One lung cancer data is downloaded from the TCGA database (https://portal.gdc.cancer.gov/projects/TCGA-LUAD), which consists of 58 normal samples and 524 tumor samples. The other lung cancer data is downloaded from the SRA database (https://www.ncbi.nlm.nih.gov/sra?term=ERP001058), which consists of 77  paired normal and tumor samples.


## Dependencies

The Detector relies on the following python libraries and R packages: 

- **[numpy](https://www.numpy.org/)** NumPy is the fundamental package for scientific computing with Python. 
- **[pandas](https://www.pandas.org/)** Pandas is a python library that allows you to easily manipulate data to analyze. 
- **[limma](http://bioconductor.org/packages/release/bioc/html/limma.html)** Data analysis, linear models and differential expression for microarray data.
- **[HTSanalyzeR](https://www.bioconductor.org/packages/release/bioc/html/HTSanalyzeR.html)** This package provides classes and methods for gene and contig set enrichment. The over-representation analysis is performed based on hypergeometric tests.


- **Step 1: Run dekupl-run**. We first activate the conda environement where dekupl-run was installed, then we run the software. The description of parameters can be found from the repository of DEkupl (https://github.com/Transipedia/dekupl-run)
    ```
    conda install -n dekupl -c transipedia dekupl-run dekupl-annotation 
    source activate dekupl
    dekupl-run --configfile my-config.json  -jNB_THREADS --resources ram=MAX_MEMORY -p
    ``` 


- **Step 2: Run dekupl-annotation**. Then we ran DEkupl annotation on the output results from both two datasets. The reference files include the Genome sequence (GRCh38.p12) and annotation file (version 31). The main output files are the DiffContigInfo.tsv which include the annotation information of each contig.The description of parameters can be found from the repository of DEkupl (https://github.com/Transipedia/dekupl-annotation)
    ```
    source activate dekupl
    dkpl index -g toy/references/GRCh38-chr22.fa.gz -a toy/references/GRCh38-chr22.gff.gz -i test_index
    dkpl annot -i test_index toy/dkpl-run/merged-diff-counts.tsv.gz
    ```

- **Step 3: Extract shared events**. We ran this pipeline using the output results from two datasets generated by DEkupl run and DEkupl annotation as input.So we can compare the consistency between two datasets from both the gene's level and contig's level.
    ```
    python3 compare_contigs.py data/dkplanno_dataset1/DiffContigsInfos.tsv data/dkplanno_dataset2/DiffContigsInfos.tsv data/dkplrun_dataset1 data/dkplrun_dataset2/ data/genome.fa
    ```
## Input files

- Table `DiffContigsInfos.tsv`, summarizing for each contig, which is the DEkupl annotation output of dataset1/2.

- Path `dkplrun_dataset1`, the output directory name of DEkupl-run for dataset1

- Path `dkplrun_dataset2`, the output directory name of DEkupl-run for dataset2

- Fasta `genome.fa`, the fasta format file of the genome data for the downstream blast analysis.

## Output files
- Figure `enrichment.pdf`, the GSEA-like enrichment result using shared events.

- Figure `jaccardidx.pdf`, the table showing the comparison between two datasets using the Jaccard index.

- Table  `shared_contigs_dataset_DiffContigInfo.tsv`, the table containing all shared events and corresponding annotation data from each dataset. 

The criteria for comparing different categories of contigs from two datasets include:

- **SNV**: position of SNV 
- **LincRNA/intron**: position of center of contig +/- 30nt
- **splice/split**: positions of both splice sites +/-30nt
- **polyA**: position of 3'end of contig +/- 10nt
- **unmapped/repeat**: build k-mer contigs and annotate contigs based on sequence alignment.




