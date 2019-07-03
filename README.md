# T Cell–Dependent Affinity Maturation and Innate Immune Pathways Differentially Drive Autoreactive B Cell Responses in Rheumatoid Arthritis

Processed scRNAseq data, flow cytometry data, and analysis scripts

## Installation instructions
Clone this repo and in R run

```r
install.packages('RAData_1.4.tar.gz')
# or build the package from source
# devtools::install("RAData")
library(RAData)
# Deduplicated reads (UMIs) that were run through RSEM
data(radedup)

# Raw reads that were run through RSEM
data(radup)

```
These objects are lists of length 3, containing members `exprs` (expression array -- RSEM expected counts), `cData` (cellular metadata), `fData` (feature/gene metadat).

If you aren't running R, you can browse under
`RSEM Tables/` to find the pre- and post- deduplicated TPM and count tables.

## SRA

If want the .fastq, browse to [the SRA project page](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA487301).
