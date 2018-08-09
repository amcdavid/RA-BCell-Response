## ----library, error=FALSE------------------------------------------------
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
loadFilesIntoMast <- function(files, normalize=function(x) log2(x+1), maxGeneId=48424){
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

    browser()
    stopifnot(all(rowSums(EE.new)==rowSums(EE)))
    
    Cdat$ngeneson <- rowSums(EE.new>0)
    Cdat$wellKey <- make.unique(Cdat$wellKey)
    list(exprs=EE.new, cData=Cdat, fData=geneid)
}

## ----loadDedup, error=FALSE----------------------------------------------
files <- normalizeFiles(list.files(file.path(dataDir, 'RSEM Tables', 'Post-deduplication RSEM'), pattern = 'rsem_count_matrix', full.names=TRUE, recursive=TRUE))
radedup <- loadFilesIntoMast(files)
cd <- as.data.table(radedup$cData)
cd[is.na(id) | is.na(treatment)]
cd[!str_detect(rawid, fixed(id))]
cd[id=='F' & is.na(subtype)]

## ------------------------------------------------------------------------
kable(cd[,c(.SD[1,], nwells=.N),keyby=list(file, rawid)])
kable(dcast.data.table(cd, id +subtype ~ treatment))

## ----loadDup, error=FALSE------------------------------------------------
files <- normalizeFiles(list.files(file.path(dataDir, 'RSEM Tables', 'Pre-deduplication RSEM'), pattern = 'rsem_count_matrix', full.names=TRUE, recursive=TRUE))
radup <- loadFilesIntoMast(files)
cd <- as.data.table(radup$cData)
cd[is.na(id) | is.na(treatment)]
cd[!str_detect(rawid, fixed(id))]
cd[id=='F' & is.na(subtype)]

## ------------------------------------------------------------------------
kable(cd[,c(.SD[1,], nwells=.N),keyby=list(file, rawid)])
kable(dcast.data.table(cd, id +subtype ~ treatment))

## ---- eval=FALSE---------------------------------------------------------
## ## Line up common samples to compare pre-and-post deduplicated
## commonWell <- intersect(radup$cData$wellKey, radedup$cData$wellKey)
## commonPid <- intersect(radup$fData$gene_id, radedup$fData$gene_id)
## common
## commonSca <- commonSca[commonWell,commonPid]
## commonSca <- addlayer(commonSca, 'dup')
## layer(commonSca) <- 'dup'
## exprs(commonSca) <- exprs(radup)[commonWell,commonPid]
## 
## sca <- FromMatrix('SingleCellAssay', EE, cData=sampleSheet, fData=geneid)
## sca <- subset(sca, !(rawid %in% c('empty', 'ercc')))
## sca

