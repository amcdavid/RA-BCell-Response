
```r
library(data.table)
library(stringr)
dataDir <- '../inst/extdata' #count tables live here
idMaybe <- c('F', 'R1', 'R48', 'R26', 'R4', 'R59', 'R72')
idMaybeRE <- paste0('(', paste0(idMaybe, collapse='|'), ')')
treatmentMaybe <- c('AG', 'CIT', 'RF')
treatmentMaybeRE <- paste0('(', paste0(treatmentMaybe, collapse='|'), ')')

## Convert counts into counts per million
cpmfun <- function(ee){
    cpm <- ee/rowSums(ee, na.rm=TRUE)*1e6
    log2(cpm+1)
}

## Remove deduplicate files 
normalizeFiles <- function(files){
    ## detect _f_ deduplicates and remove
    filebase <- str_replace(files, '_f', '')
    filef <- str_detect(files, '_f')
    dtfile <- data.table(files, filebase, filef)
    dtfile[,N:=.N,keyby=filebase]
    ## is a singleton or is not _f_
    dtfile[N==1 | (N>1 & filef), files]
}

## Concatenate expression matrices, and cell data
## Intersect features
loadFilesIntoMast <- function(files, normalize=function(x) x, maxGeneId=48424){
    ncells <- ifelse(str_detect(files, 'Bulk'), 1000, 1)
    ## Manually added columns treatment and iid to this file...need to update this as new samples come in if their naming pattern doesn't match
    sampleSheet <- fread(file.path(dataDir, 'sampleSheet3.csv'), colClasses=c('subtype'='character', 'wellID'='character'))

    geneid <- Cdat <- EE <- vector('list', length=length(files))
    for( i in seq_along(files)){
        ff <- fread(files[[i]])
        setkey(ff, gene_id)
        ## genes present in file i
        geneid[[i]] <- as.vector(ff[,gene_id])
        ## take outer join
        ff <- ff[list(seq_len(maxGeneId)),,nomatch=NA]
        ## delete extra columns if present
        suppressWarnings(ff[, ':='(entrez_id=NULL, transcript_id=NULL, gene_symbol=NULL)])
        mat <- ff[,-1,with=FALSE]
        ee <- normalize(t(mat))
        sampid <- rownames(ee)
        this.ncells <- ncells[i]
        if(this.ncells==1){
            ## pull rawid (contains id and treatment) and wellID
            cdat <- data.frame(str_split_fixed(sampid, '_', 2), stringsAsFactors=FALSE)
            names(cdat) <- c('rawid', 'wellID')
            cdat$rawid <- toupper(cdat$rawid)
            idmaybe <- str_match(cdat$rawid, idMaybeRE)[,2]
            treatmaybe <- str_match(cdat$rawid, treatmentMaybeRE)[,2]
            if(all(!is.na(idmaybe))){
                cdat$id <- idmaybe
                if(all(!is.na(treatmaybe))){
                    cdat$treatment <- treatmaybe
                } else{
                    cdat$treatment <- 'AG'
                }
            }
            ## else leave null and fill in from sampleSheet
            ## cdat$treatment <- idtreat[,3]
            ## cdat$treatment[is.na(cdat$treatment)] <- 'AG'
        } else{
            cdat <- data.frame(str_split_fixed(sampid, '_', 3), stringsAsFactors=FALSE)
            names(cdat) <- c('rawid', 'subtype', 'wellID')
        }
            cdat$ncells <- this.ncells
            cdat$file <- basename(files[[i]])
            cdat$idx <- 1:nrow(cdat)
        cdat2 <- merge(cdat, sampleSheet, all.x=TRUE)
        setorder(cdat2, idx)
        cdat2$wellKey <- with(cdat2, paste(rawid, id, subtype, wellID, sep=':'))
        cdat2$idx <- NULL
        rownames(ee) <- cdat2$wellKey
        EE[[i]] <- ee
        Cdat[[i]] <- cdat2
    }
    geneidmin <- geneid[[1]]
    for(i in seq_along(files)){
        geneidmin <- intersect(geneidmin, geneid[[i]])
    }

    ##dt <- data.table(t(sapply(EE, dim)), files)
    for(i in seq_along(files)){
        EE[[i]] <- EE[[i]][,geneidmin]
    }
    EE <- do.call(rbind, EE)
    stopifnot(all(!is.na(EE)))
    Cdat <- do.call(rbind, Cdat)
    

    litGenes <- fread(file.path(dataDir, 'Subset genes.csv'))[, 1:4, with=FALSE]
    setnames(litGenes, c('naive', 'memory', 'asc', 'cd3'))
    litGenes <- melt.data.table(litGenes, measure.vars=names(litGenes))[!is.na(value),]
    setnames(litGenes, 'value', 'gene_id')
    litGenes <- litGenes[,list(litSubtype=paste(variable, collapse=':')),keyby=list(gene_id)]
    genetable <- fread(file.path(dataDir, 'RSEM Tables', 'Post-deduplication RSEM', 'Bulk', 'rsem_fdata.csv'))
    setkey(genetable, gene_id)
    geneid <- genetable[list(geneidmin),]
    geneid[,symbolOrIdx:=make.names(gene_symbol)]
    geneid[is.na(gene_symbol),symbolOrIdx:=paste0('p',as.character(gene_id))]
    geneid[,entrezOrIdx:=paste0('EN', entrez_id)]
    geneid[is.na(entrez_id),entrezOrIdx:=paste0('p',as.character(gene_id))]
    geneid <- merge(geneid, litGenes, by='gene_id', all.x=TRUE)

    btmGenes <- gdata::read.xls(file.path(dataDir, 'btm_annotation_table_InterestingModules 2016-05-27.xls'), as.is=TRUE)
    setDT(btmGenes)
    btmGenes[,Interest:=factor(Interest, levels=c('veryhigh', 'high', 'medium', 'low'))]
    btmGenesFlat <- btmGenes[,{
    genelist <- str_split(Module.member.genes, fixed(','))[[1]]
    ng <- length(genelist)
    .(Interest=rep(Interest, ng), ID=rep(ID, ng), Module.title=rep(Module.title, ng), gene_symbol=genelist)
    },by=.(seq_len(nrow(btmGenes)))]
    setorder(btmGenesFlat, gene_symbol, Interest)
    setkey(btmGenesFlat, gene_symbol)
    ordered <- merge(geneid, unique(btmGenesFlat[,.(gene_symbol, Interest, ID, Module.title)]),by.y='gene_symbol', by.x='symbolOrIdx', all.x=TRUE)
    geneid <- ordered[gene_id %in% geneidmin]
    setorder(geneid, gene_id)
    stopifnot(nrow(geneid)==ncol(EE))
    fVars <- geneid[, symbolOrIdx]
    geneset <- split(seq_along(fVars), fVars)
    fdOrder <- order(unlist(lapply(geneset, min)))
    geneset <- geneset[fdOrder]
    genefirst <- unlist(lapply(geneset, function(x) x[1]))
    skeleton <- EE[,genefirst]
    exprs.orig <- EE
    exprs.new <- lapply(geneset, function(colidx) {
        cols <- exprs.orig[, colidx, drop = FALSE]
        rowSums(cols)
    })
    ## This sets the colnames to names(exprs.new), hence symbolOrIdx
    EE.new <- do.call(cbind, exprs.new)
    geneid[,primerid:=symbolOrIdx]
    ## sort the fvars
    geneid <- geneid[genefirst,]

    stopifnot(all(abs(rowSums(EE.new)-rowSums(EE))<.01))
    
    Cdat$ngeneson <- rowSums(EE.new>0)
    Cdat$wellKey <- make.unique(Cdat$wellKey)
    Cdat$counts <- rowSums(EE.new)

    isotype_files <- list.files(file.path(dataDir, 'IGHC Counts'), pattern="*.csv", full.names=TRUE)
    itf <- str_match(list.files(file.path(dataDir, 'IGHC Counts'), pattern="*.csv"), '[A-Za-z0-9]+')
    names(isotype_files) <- itf
    #isotype_cw <- data.table(fname=itf, id=str_c('R', str_match(itf, '[0-9]+')), treatment=str_match(itf, 'a|r|c'))
    isotype_cw <- data.table(fname=c("1a", "235", "26a", "26c", "26r", 
                                     "48a", "48c", "4a", "4r", "59a", "59c", "59r", "72a", "72c", 
                                     "72r", "Fa"),
                             id = c("R1", "F", "R26", "R26", "R26", "R48", 
                                    "R48", "R4", "R4", "R59", "R59", "R59", "R72", "R72", "R72", 'F'),
                             treatment = c("AG", 'none', "AG", "CIT", "RF", "AG", "CIT", "AG", 
                                           "RF", "AG", "CIT", "RF", "AG", "CIT", "RF", "AG"))
    
    isotype <- rbindlist(lapply(isotype_files, fread), fill=TRUE, idcol='fname')
    isotype <- merge(isotype, isotype_cw, by='fname')
    isotype[,wellID:=str_match(cellID, '[0-9]+$')]
    Cdat <- merge(Cdat, isotype, c('id', 'treatment', 'wellID'), all.x=TRUE)
    Cdat <- Cdat[match(rownames(EE.new), Cdat$wellKey),]
    stopifnot(all(Cdat$wellKey==row.names(EE.new)))
    list(exprs=EE.new, cData=Cdat, fData=geneid)
}
```


```r
files <- normalizeFiles(list.files(file.path(dataDir, 'RSEM Tables', 'Post-deduplication RSEM'), pattern = 'rsem_count_matrix', full.names=TRUE, recursive=TRUE))
radedup <- loadFilesIntoMast(files)
cd <- as.data.table(radedup$cData)
cd[is.na(id) | is.na(treatment)]
```

```
## Empty data.table (0 rows) of 19 cols: id,treatment,wellID,rawid,subtype,ncells...
```

```r
cd[!str_detect(rawid, fixed(id))]
```

```
##      id treatment wellID rawid subtype ncells
##   1:  F      none    304     M       M      1
##   2:  F      none    305     M       M      1
##   3:  F      none    306     M       M      1
##   4:  F      none    307     M       M      1
##   5:  F      none    312     M       M      1
##                                     file   wellKey ngeneson    counts
##   1: rsem_count_matrix_f_235sc_RSEM2.csv M:F:M:304      915   4731.76
##   2: rsem_count_matrix_f_235sc_RSEM2.csv M:F:M:305      801  17054.09
##   3: rsem_count_matrix_f_235sc_RSEM2.csv M:F:M:306     2289 170615.63
##   4: rsem_count_matrix_f_235sc_RSEM2.csv M:F:M:307      501   2455.76
##   5: rsem_count_matrix_f_235sc_RSEM2.csv M:F:M:312     1899  97562.28
##      fname   cellID IgM IgD IgG IgE IgA Subset Isotype
##   1:   235 235sc304   1   0   3   0   0      M       G
##   2:   235 235sc305   0   0   0   0   4      M       A
##   3:   235 235sc306   0   0   0   0   0      M unknown
##   4:   235 235sc307   0   0   0   0   0      M unknown
##   5:   235 235sc312   0   0   1   0   4      M       A
##  [ reached getOption("max.print") -- omitted 6 rows ]
```

```r
cd[id=='F' & is.na(subtype)]
```

```
##      id treatment wellID rawid subtype ncells
##   1:  F        AG    208     F      NA      1
##   2:  F        AG    209     F      NA      1
##   3:  F        AG    210     F      NA      1
##   4:  F        AG    211     F      NA      1
##   5:  F        AG    212     F      NA      1
##                                     file    wellKey ngeneson     counts
##   1: rsem_count_matrix_f_FCHAg_RSEM2.csv F:F:NA:208      624   11473.76
##   2: rsem_count_matrix_f_FCHAg_RSEM2.csv F:F:NA:209      360    7620.00
##   3: rsem_count_matrix_f_FCHAg_RSEM2.csv F:F:NA:210      173    2166.38
##   4: rsem_count_matrix_f_FCHAg_RSEM2.csv F:F:NA:211      101    1307.60
##   5: rsem_count_matrix_f_FCHAg_RSEM2.csv F:F:NA:212      163    1945.64
##      fname cellID IgM IgD IgG IgE IgA Subset Isotype
##   1:    Fa  Fa208   0   0  65   0   0     NA      NA
##   2:    Fa  Fa209   0   0   0   0   0     NA      NA
##   3:    Fa  Fa210   0   0   0   0   0     NA      NA
##   4:    Fa  Fa211   0   0   0   0   0     NA      NA
##   5:    Fa  Fa212   0   0   1   0   0     NA      NA
##  [ reached getOption("max.print") -- omitted 6 rows ]
```
Inferred IDs (deduplicated)


```r
kable(cd[,c(.SD[1,], nwells=.N),keyby=list(file, rawid)])
```



|file                                            |rawid  |id  |treatment |wellID |subtype | ncells|wellKey           | ngeneson|     counts|fname |cellID   | IgM| IgD| IgG| IgE| IgA|Subset |Isotype | nwells|
|:-----------------------------------------------|:------|:---|:---------|:------|:-------|------:|:-----------------|--------:|----------:|:-----|:--------|---:|---:|---:|---:|---:|:------|:-------|------:|
|rsem_count_matrix_f_235b_RSEM2.csv              |235    |235 |none      |217    |M       |   1000|235:235:M:217     |    13217|  659931.39|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      4|
|rsem_count_matrix_f_235b_RSEM2.csv              |235-1  |235 |none      |265    |M       |   1000|235-1:235:M:265   |    13028|  793778.49|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      4|
|rsem_count_matrix_f_235sc_RSEM2.csv             |M      |F   |none      |304    |M       |      1|M:F:M:304         |      915|    4731.76|235   |235sc304 |   1|   0|   3|   0|   0|M      |G       |     48|
|rsem_count_matrix_f_235sc_RSEM2.csv             |N      |F   |none      |212    |N       |      1|N:F:N:212         |     3096|  140835.01|235   |235sc212 |   1|   0|   0|   0|   0|N      |M       |     46|
|rsem_count_matrix_f_235sc_RSEM2.csv             |P      |F   |none      |208    |P       |      1|P:F:P:208         |     3970|  385961.26|235   |235sc208 |   1|   0|   9|   0|  61|P      |A       |     48|
|rsem_count_matrix_f_235sc_RSEM2.csv             |T      |F   |none      |308    |T       |      1|T:F:T:308         |     2431|  103116.91|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |     45|
|rsem_count_matrix_f_832b_833b_withBTK_RSEM2.csv |832    |832 |none      |281    |M       |   1000|832:832:M:281     |    12274| 1326131.84|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      8|
|rsem_count_matrix_f_832b_833b_withBTK_RSEM2.csv |833    |833 |none      |280    |M       |   1000|833:833:M:280     |    10808|  652524.17|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      8|
|rsem_count_matrix_f_FCHAg_RSEM2.csv             |F      |F   |AG        |208    |NA      |      1|F:F:NA:208        |      624|   11473.76|Fa    |Fa208    |   0|   0|  65|   0|   0|NA     |NA      |    190|
|rsem_count_matrix_f_JG1Ag_RSEM2.csv             |R1     |R1  |AG        |208    |NA      |      1|R1:R1:NA:208      |     2435|   26941.42|1a    |1a208    |   0|   0|   0|   0|   0|NA     |NA      |    190|
|rsem_count_matrix_f_JG26A_RSEM2.csv             |R26AG  |R26 |AG        |208    |NA      |      1|R26AG:R26:NA:208  |     5280|  445207.09|26a   |26a208   |   0|   1|   0|   0|  46|NA     |NA      |    186|
|rsem_count_matrix_f_JG26cit_RSEM2.csv           |R26CIT |R26 |CIT       |208    |NA      |      1|R26CIT:R26:NA:208 |     4303|   49531.28|26c   |26c208   |   1|   0|   0| 109|  12|NA     |NA      |    186|
|rsem_count_matrix_f_JG26rf_RSEM2.csv            |R26RF  |R26 |RF        |208    |NA      |      1|R26RF:R26:NA:208  |     2466|  103553.29|26r   |26r208   |  55|   0|   1|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG48Ag_RSEM2.csv            |R48    |R48 |AG        |208    |NA      |      1|R48:R48:NA:208    |     4405|  101890.22|48a   |48a208   |   6|   0|   2|   0|   0|NA     |NA      |    189|
|rsem_count_matrix_f_JG48cit_RSEM2.csv           |R48CIT |R48 |CIT       |208    |NA      |      1|R48CIT:R48:NA:208 |     4144|   87164.67|48c   |48c208   |   1|   0|   0|   0|   7|NA     |NA      |    144|
|rsem_count_matrix_f_JG4Ag_RSEM2.csv             |R4     |R4  |AG        |208    |NA      |      1|R4:R4:NA:208      |     3976|  118497.07|4a    |4a208    | 185|   0|  42|   0|  22|NA     |NA      |    186|
|rsem_count_matrix_f_JG4rf_RSEM2.csv             |R4RF   |R4  |RF        |208    |NA      |      1|R4RF:R4:NA:208    |     3020|  291678.14|4r    |4r208    |  86|   0|   0|   0|   0|NA     |NA      |    188|
|rsem_count_matrix_f_JG59A_RSEM2.csv             |R59AG  |R59 |AG        |208    |NA      |      1|R59AG:R59:NA:208  |     3903|   76911.79|59a   |59a208   |  11|  64|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG59cit_RSEM2.csv           |R59CIT |R59 |CIT       |208    |NA      |      1|R59CIT:R59:NA:208 |     2634|   63895.59|59c   |59c208   |  45|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG59rf_RSEM2.csv            |R59RF  |R59 |RF        |208    |NA      |      1|R59RF:R59:NA:208  |     1692|   25516.88|59r   |59r208   |  15|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG72Ag_RSEM2.csv            |R72AG  |R72 |AG        |208    |NA      |      1|R72AG:R72:NA:208  |     3841|  219173.33|72a   |72a208   | 486|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG72cit_RSEM2.csv           |R72CIT |R72 |CIT       |208    |NA      |      1|R72CIT:R72:NA:208 |       74|     728.80|72c   |72c208   |   0|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG72rf_RSEM2.csv            |R72RF  |R72 |RF        |208    |NA      |      1|R72RF:R72:NA:208  |     2396|  227378.19|72r   |72r208   |   0|   0|   0|   0|   0|NA     |NA      |    186|

```r
kable(dcast.data.table(cd, id +subtype ~ treatment))
```

```
## Using 'Isotype' as value column. Use 'value.var' to override
```

```
## Aggregate function missing, defaulting to 'length'
```



|id  |subtype |  AG| BTK| CIT|  RF| none|
|:---|:-------|---:|---:|---:|---:|----:|
|235 |M       |   0|   0|   0|   0|    2|
|235 |N       |   0|   0|   0|   0|    2|
|235 |P       |   0|   0|   0|   0|    2|
|235 |T       |   0|   0|   0|   0|    2|
|832 |M       |   0|   1|   0|   0|    1|
|832 |N       |   0|   1|   0|   0|    1|
|832 |P       |   0|   1|   0|   0|    1|
|832 |T       |   0|   1|   0|   0|    1|
|833 |M       |   0|   1|   0|   0|    1|
|833 |N       |   0|   1|   0|   0|    1|
|833 |P       |   0|   1|   0|   0|    1|
|833 |T       |   0|   1|   0|   0|    1|
|F   |M       |   0|   0|   0|   0|   48|
|F   |N       |   0|   0|   0|   0|   46|
|F   |P       |   0|   0|   0|   0|   48|
|F   |T       |   0|   0|   0|   0|   45|
|F   |NA      | 190|   0|   0|   0|    0|
|R1  |NA      | 190|   0|   0|   0|    0|
|R26 |NA      | 186|   0| 186| 186|    0|
|R4  |NA      | 186|   0|   0| 188|    0|
|R48 |NA      | 189|   0| 144|   0|    0|
|R59 |NA      | 186|   0| 186| 186|    0|
|R72 |NA      | 186|   0| 186| 186|    0|
Number of wells per ID and subtype (deduplicated)


```r
with(subset(cd, is.na(IgM)), table(id, treatment))
```

```
##      treatment
## id    BTK CIT none
##   235   0   0    8
##   832   4   0    4
##   833   4   0    4
##   F     0   0   45
##   R48   0   1    0
```
Merge with isotype info




```r
files <- normalizeFiles(list.files(file.path(dataDir, 'RSEM Tables', 'Pre-deduplication RSEM'), pattern = 'rsem_count_matrix', full.names=TRUE, recursive=TRUE))
radup <- loadFilesIntoMast(files)
cd <- as.data.table(radup$cData)
cd[is.na(id) | is.na(treatment)]
```

```
## Empty data.table (0 rows) of 19 cols: id,treatment,wellID,rawid,subtype,ncells...
```

```r
cd[!str_detect(rawid, fixed(id))]
```

```
##      id treatment wellID rawid subtype ncells
##   1:  F      none    304     M       M      1
##   2:  F      none    305     M       M      1
##   3:  F      none    306     M       M      1
##   4:  F      none    307     M       M      1
##   5:  F      none    312     M       M      1
##                                     file   wellKey ngeneson    counts
##   1: rsem_count_matrix_f_235sc_RSEM1.csv M:F:M:304      906   5731.55
##   2: rsem_count_matrix_f_235sc_RSEM1.csv M:F:M:305      593   4340.53
##   3: rsem_count_matrix_f_235sc_RSEM1.csv M:F:M:306     2136 269676.14
##   4: rsem_count_matrix_f_235sc_RSEM1.csv M:F:M:307      373   1398.24
##   5: rsem_count_matrix_f_235sc_RSEM1.csv M:F:M:312     1781 142000.60
##      fname   cellID IgM IgD IgG IgE IgA Subset Isotype
##   1:   235 235sc304   1   0   3   0   0      M       G
##   2:   235 235sc305   0   0   0   0   4      M       A
##   3:   235 235sc306   0   0   0   0   0      M unknown
##   4:   235 235sc307   0   0   0   0   0      M unknown
##   5:   235 235sc312   0   0   1   0   4      M       A
##  [ reached getOption("max.print") -- omitted 6 rows ]
```

```r
cd[id=='F' & is.na(subtype)]
```

```
##      id treatment wellID rawid subtype ncells
##   1:  F        AG    208   FAG      NA      1
##   2:  F        AG    209   FAG      NA      1
##   3:  F        AG    210   FAG      NA      1
##   4:  F        AG    211   FAG      NA      1
##   5:  F        AG    212   FAG      NA      1
##                                   file      wellKey ngeneson     counts
##   1: rsem_count_matrix_FCHAg_RSEM1.csv FAG:F:NA:208      530   19999.99
##   2: rsem_count_matrix_FCHAg_RSEM1.csv FAG:F:NA:209      323   14534.34
##   3: rsem_count_matrix_FCHAg_RSEM1.csv FAG:F:NA:210      143    3638.14
##   4: rsem_count_matrix_FCHAg_RSEM1.csv FAG:F:NA:211       78    2308.88
##   5: rsem_count_matrix_FCHAg_RSEM1.csv FAG:F:NA:212      161    3528.01
##      fname cellID IgM IgD IgG IgE IgA Subset Isotype
##   1:    Fa  Fa208   0   0  65   0   0     NA      NA
##   2:    Fa  Fa209   0   0   0   0   0     NA      NA
##   3:    Fa  Fa210   0   0   0   0   0     NA      NA
##   4:    Fa  Fa211   0   0   0   0   0     NA      NA
##   5:    Fa  Fa212   0   0   1   0   0     NA      NA
##  [ reached getOption("max.print") -- omitted 6 rows ]
```
Inferred IDs (duplicates)


```r
kable(cd[,c(.SD[1,], nwells=.N),keyby=list(file, rawid)])
```



|file                                        |rawid  |id  |treatment |wellID |subtype | ncells|wellKey           | ngeneson|     counts|fname |cellID   | IgM| IgD| IgG| IgE| IgA|Subset |Isotype | nwells|
|:-------------------------------------------|:------|:---|:---------|:------|:-------|------:|:-----------------|--------:|----------:|:-----|:--------|---:|---:|---:|---:|---:|:------|:-------|------:|
|rsem_count_matrix_235B_RSEM1.csv            |235    |235 |none      |225    |P       |   1000|235:235:P:225     |    14855| 7049389.07|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      4|
|rsem_count_matrix_235B_RSEM1.csv            |235-1  |235 |none      |273    |P       |   1000|235-1:235:P:273   |    12160| 1105615.19|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      4|
|rsem_count_matrix_832bulk_withBTK_RSEM1.csv |832    |832 |none      |281    |M       |   1000|832:832:M:281     |    12023| 2829071.58|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      8|
|rsem_count_matrix_833bulk_withBTK_RSEM1.csv |833    |833 |none      |280    |M       |   1000|833:833:M:280     |    10439| 1297422.35|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |      7|
|rsem_count_matrix_FCHAg_RSEM1.csv           |FAG    |F   |AG        |208    |NA      |      1|FAG:F:NA:208      |      530|   19999.99|Fa    |Fa208    |   0|   0|  65|   0|   0|NA     |NA      |    188|
|rsem_count_matrix_JG26Ag_RSEM1.csv          |R26AG  |R26 |AG        |208    |NA      |      1|R26AG:R26:NA:208  |     4963| 1744702.25|26a   |26a208   |   0|   1|   0|   0|  46|NA     |NA      |    186|
|rsem_count_matrix_JG26rf_RSEM1.csv          |R26RF  |R26 |RF        |208    |NA      |      1|R26RF:R26:NA:208  |     2340|  255342.03|26r   |26r208   |  55|   0|   1|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_JG72Ag_RSEM1.csv          |R72AG  |R72 |AG        |208    |NA      |      1|R72AG:R72:NA:208  |     3635|  535366.77|72a   |72a208   | 486|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_JG72RF_RSEM1.csv          |R72RF  |R72 |RF        |208    |NA      |      1|R72RF:R72:NA:208  |     2410|  412312.05|72r   |72r208   |   0|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_R1Ag_RSEM1.csv            |R1AG   |R1  |AG        |208    |NA      |      1|R1AG:R1:NA:208    |     2263|   26041.37|1a    |1a208    |   0|   0|   0|   0|   0|NA     |NA      |    190|
|rsem_count_matrix_R48Ag_RSEM1.csv           |R48AG  |R48 |AG        |208    |NA      |      1|R48AG:R48:NA:208  |     4288|  181487.92|48a   |48a208   |   6|   0|   2|   0|   0|NA     |NA      |    189|
|rsem_count_matrix_R48cit_RSEM1.csv          |R48CIT |R48 |CIT       |208    |NA      |      1|R48CIT:R48:NA:208 |     3974|  216832.28|48c   |48c208   |   1|   0|   0|   0|   7|NA     |NA      |    144|
|rsem_count_matrix_R4Ag_RSEM1.csv            |R4AG   |R4  |AG        |208    |NA      |      1|R4AG:R4:NA:208    |     3869|  217133.82|4a    |4a208    | 185|   0|  42|   0|  22|NA     |NA      |    186|
|rsem_count_matrix_R4rf_RSEM1.csv            |R4RF   |R4  |RF        |208    |NA      |      1|R4RF:R4:NA:208    |     3008|  444974.79|4r    |4r208    |  86|   0|   0|   0|   0|NA     |NA      |    188|
|rsem_count_matrix_f_235sc_RSEM1.csv         |M      |F   |none      |304    |M       |      1|M:F:M:304         |      906|    5731.55|235   |235sc304 |   1|   0|   3|   0|   0|M      |G       |     48|
|rsem_count_matrix_f_235sc_RSEM1.csv         |N      |F   |none      |212    |N       |      1|N:F:N:212         |     3002|  206405.80|235   |235sc212 |   1|   0|   0|   0|   0|N      |M       |     46|
|rsem_count_matrix_f_235sc_RSEM1.csv         |P      |F   |none      |208    |P       |      1|P:F:P:208         |     3532|  275502.30|235   |235sc208 |   1|   0|   9|   0|  61|P      |A       |     48|
|rsem_count_matrix_f_235sc_RSEM1.csv         |T      |F   |none      |308    |T       |      1|T:F:T:308         |     2304|  164164.96|NA    |NA       |  NA|  NA|  NA|  NA|  NA|NA     |NA      |     45|
|rsem_count_matrix_f_JG59Ag_RSEM1.csv        |R59AG  |R59 |AG        |208    |NA      |      1|R59AG:R59:NA:208  |     3762|  179750.17|59a   |59a208   |  11|  64|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG59cit_RSEM1.csv       |R59CIT |R59 |CIT       |208    |NA      |      1|R59CIT:R59:NA:208 |     2573|  104045.24|59c   |59c208   |  45|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG59rf_RSEM1.csv        |R59RF  |R59 |RF        |208    |NA      |      1|R59RF:R59:NA:208  |     1657|   37330.49|59r   |59r208   |  15|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_JG72cit_RSEM1.csv       |R72CIT |R72 |CIT       |208    |NA      |      1|R72CIT:R72:NA:208 |       80|    1813.27|72c   |72c208   |   0|   0|   0|   0|   0|NA     |NA      |    186|
|rsem_count_matrix_f_R26cit_RSEM1.csv        |R26CIT |R26 |CIT       |208    |NA      |      1|R26CIT:R26:NA:208 |     4212|   92898.81|26c   |26c208   |   1|   0|   0| 109|  12|NA     |NA      |    186|

```r
kable(dcast.data.table(cd, id +subtype ~ treatment))
```

```
## Using 'Isotype' as value column. Use 'value.var' to override
```

```
## Aggregate function missing, defaulting to 'length'
```



|id  |subtype |  AG| BTK| CIT|  RF| none|
|:---|:-------|---:|---:|---:|---:|----:|
|235 |M       |   0|   0|   0|   0|    2|
|235 |N       |   0|   0|   0|   0|    2|
|235 |P       |   0|   0|   0|   0|    2|
|235 |T       |   0|   0|   0|   0|    2|
|832 |M       |   0|   1|   0|   0|    1|
|832 |N       |   0|   1|   0|   0|    1|
|832 |P       |   0|   1|   0|   0|    1|
|832 |T       |   0|   1|   0|   0|    1|
|833 |M       |   0|   1|   0|   0|    1|
|833 |N       |   0|   0|   0|   0|    1|
|833 |P       |   0|   1|   0|   0|    1|
|833 |T       |   0|   1|   0|   0|    1|
|F   |M       |   0|   0|   0|   0|   48|
|F   |N       |   0|   0|   0|   0|   46|
|F   |P       |   0|   0|   0|   0|   48|
|F   |T       |   0|   0|   0|   0|   45|
|F   |NA      | 188|   0|   0|   0|    0|
|R1  |NA      | 190|   0|   0|   0|    0|
|R26 |NA      | 186|   0| 186| 186|    0|
|R4  |NA      | 186|   0|   0| 188|    0|
|R48 |NA      | 189|   0| 144|   0|    0|
|R59 |NA      | 186|   0| 186| 186|    0|
|R72 |NA      | 186|   0| 186| 186|    0|
Number of wells per ID and subtype (duplicates)


```r
with(subset(cd, is.na(IgM)), table(id, treatment))
```

```
##      treatment
## id    BTK CIT none
##   235   0   0    8
##   832   4   0    4
##   833   3   0    4
##   F     0   0   45
##   R48   0   1    0
```
Merge with isotype info



```r
## Line up common samples to compare pre-and-post deduplicated
commonWell <- intersect(radup$cData$wellKey, radedup$cData$wellKey)
commonPid <- intersect(radup$fData$gene_id, radedup$fData$gene_id)
common
commonSca <- commonSca[commonWell,commonPid]
commonSca <- addlayer(commonSca, 'dup')
layer(commonSca) <- 'dup'
exprs(commonSca) <- exprs(radup)[commonWell,commonPid]

sca <- FromMatrix('SingleCellAssay', EE, cData=sampleSheet, fData=geneid)
sca <- subset(sca, !(rawid %in% c('empty', 'ercc')))
sca
```
