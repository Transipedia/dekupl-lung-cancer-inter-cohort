# Shared Events Detector

Exhaustive capture of biological variation in RNA-seq data through k-mer decomposition (article: https://doi.org/10.1186/s13059-017-1372-2, pre-print: http://biorxiv.org/content/early/2017/06/02/122937).

DE-kupl is a computational protocol that aims to capture all k-mer variation in an input set of RNA-seq libraries. To verify the reproducibility of DE-kupl, we developed this pipeline to compare the consistency of events between different cohorts.
The criteria for comparing different categories of contigs from two datasets include:
- **SNV**: position of SNV 
- **LincRNA/intron**: position of center of contig +/- 30nt
- **splice/split**: positions of both splice sites +/-30nt
- **polyA**: position of 3'end of contig +/- 10nt
- **unmapped/repeat**: build k-mer contigs and annotate contigs based on sequence alignment.

## Dependencies

The Detector relies on the following python libraries and R packages: 

- **[numpy](https://www.numpy.org/)** NumPy is the fundamental package for scientific computing with Python. 
- **[pandas](https://www.pandas.org/)** Pandas is a python library that allows you to easily manipulate data to analyze. 
- **[limma](http://bioconductor.org/packages/release/bioc/html/limma.html)** Data analysis, linear models and differential expression for microarray data.
- **[HTSanalyzeR](https://www.bioconductor.org/packages/release/bioc/html/HTSanalyzeR.html)** This package provides classes and methods for gene and contig set enrichment. The over-representation analysis is performed based on hypergeometric tests.



## Usage example
```
python3 compare_contigs.py DiffContigsInfos_dataset1.tsv DiffContigsInfos_dataset2.tsv dkplrundir_dataset1 dkplrundir_dataset2 genome.fa
```
## Input files

- Table `DiffContigsInfos_dataset1.tsv`, summarizing for each contig, which is the DEkupl annotation output of dataset1.

- Table `DiffContigsInfos_dataset2.tsv`, summarizing for each contig, which is the DEkupl annotation output of dataset2.

- Path `dkplrundir_dataset1`, the output directory name of DEkupl-run for dataset1

- Path `dkplrundir_dataset2`, the output directory name of DEkupl-run for dataset2

- Fasta `genome.fa`, the fasta format file of the genome data for the downstream blast analysis.

## Output files
- Figure `enrichment.pdf`, the GSEA-like enrichment result using shared events.

- Figure `jaccardidx.pdf`, the table showing the comparison between two datasets using jaccard index.

- Table  `shared_contigs_dataset_DiffContigInfo.tsv`, the table containing all shared events and corresponding annotation data from each dataset. 


