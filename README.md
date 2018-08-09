# Processed count data and scripts for "T Cell–Dependent Affinity Maturation and Innate Immune
Pathways Differentially Drive Autoreactive B Cell Responses in
Rheumatoid Arthritis"

## Installation instructions
Clone this repo and in R run

```r
devtools::install('RAData')
library(RAData)
# Deduplicated reads (UMIs) that were run through RSEM
data(radedup)

# Raw reads that were run through RSEM
data(radup)
```

If you aren't running R, you could browse under
RAData/inst/extdata/ to find the pre- and post- deduplicated TPM and count tables.