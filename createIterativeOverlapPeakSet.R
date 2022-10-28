#Originally written by Jeff Granja and Ryan Corces
#From: Corces & Granja et. al. Science 2018
################################################################################################################
#Load Required Packages
#--------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(BSgenome))

options(scipen = 999)
################################################################################################################
# Inputs -- 
#--------------
## uses optparse package to read in command line arguments
option_list <- list(
  make_option(c("--metadata"), action="store", type="character", default=NULL),
  make_option(c("--macs2dir"), action="store", type="character", default=NULL),
  make_option(c("--outdir"), action="store", type="character", default=NULL),
  make_option(c("--suffix"), action="store", type="character", default=NULL),
  make_option(c("--blacklist"), action="store", type="character", default=NULL),
  make_option(c("--genome"), action="store", type="character", default=NULL),
  make_option(c("--spm"), action="store", type="double", default=5),
  make_option(c("--rule"), action="store", type="character", default="2"),
  make_option(c("--extend"), action="store", type="integer", default=250)
)

opt = parse_args(OptionParser(option_list=option_list))

if (is.null(opt$metadata) || is.null(opt$macs2dir) || is.null(opt$outdir) || is.null(opt$suffix) || is.null(opt$blacklist) || is.null(opt$genome) || is.null(opt$spm) || is.null(opt$rule)) {
  # print usage information if required inputs are not properly specified
  cat("\n Usage: Rscript /path/to/ArchR_Create_Extended_PeakSet.R \n
  --metadata </path/to/metadata.txt> \n
  --macs2dir </path/to/macs2/dir/with/summit.bed> \n
  --outdir </path/to/output/dir> \n
  --suffix <suffix of MACS2 summit files> (typically \"_summits.bed\") \n
  --blacklist </path/to/blacklist.bed> \n
  --genome <genome name> (from BSGenome, i.e. \"hg38\") \n
  --spm <score per million cutoff> (recommended between 2 and 5 as a starting point), default 5 \n
  --rule <character string denoting the selection rule for peak filtering> \n
  --extend <integer of how many base pairs on either side of the summit to extend the peak>, default 250 \n\n")
  stop("ERROR --- check command line input")
}

write(paste0("\nCreating merged peak sets in output directory ",opt$outdir), stdout())

metadata <- opt$metadata
macs2_dir <- opt$macs2dir
out_dir <- opt$outdir
end <- opt$suffix
blacklist_bed <- opt$blacklist
genome <- opt$genome
spm <- opt$spm
selectionRules <- opt$rule
extend <- opt$extend

################################################################################################################
#FUNCTIONS
#--------------
extendedPeakSet <- function(df, BSgenome = NULL, blacklist = NULL, extend = 250, scorePerMillion = 2, selectionRules = "(n+1)/2"){

  #a lot faster than import.bed from rtracklayer
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  
  #Error-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  # stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #--------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #--------------
  
  chromSizes <- getChromSizes(BSgenome = BSgenome)
  groups <- unique(df$groups)

  groupSEList <- lapply(seq_along(groups), function(i){
    
    message(groups[i])

    df_group = df[which(df$groups==groups[i]),]

    #--- 1. read in summit files
    grList <- GenomicRanges::GRangesList(lapply(paste0(df_group$summits),readSummits))
    
    #--- 2. resize
    #--- 3. within chromsizes
    #--- 4. not in blacklist
    #--- 5. non-overlapping
    #--- 6. score-per-million
    grList <- lapply(seq_along(grList),function(x){
      extended_summits <- grList[[x]] %>%
        resize(., width = 2 * extend + 1, fix = "center") %>%     
        subsetByOverlaps(.,chromSizes,type="within") %>%
        subsetByOverlaps(.,blacklist,invert=TRUE) %>%
        convergeClusterGRanges(., by = "score", decreasing = T)
      mcols(extended_summits)$score <- edgeR::cpm(mcols(extended_summits)$score)
      return(extended_summits)
    })
    
    #--- 7. non-overlapping overlaps determined score-per-million
    grNonOverlapping <- GenomicRanges::GRangesList(grList) %>%
      unlist %>% 
      convergeClusterGRanges(., by = "score", decreasing = T) %>%
      sortSeqlevels %>% 
      sort
    
    #--- 8. select those meeting score-per-million threshold
    grNonOverlapping <- grNonOverlapping[which(mcols(grNonOverlapping)$score > scorePerMillion),]
    
    #--- 9. get score-per-million for each sample
    scores <- lapply(seq_along(grList), function(x){
      columnOverlaps(query = grNonOverlapping, 
                     subject = grList[[x]], 
                     colname = "score", 
                     decreasing = TRUE) 
    }) %>% data.frame()
    colnames(scores) <- df_group$samples
      
    #--- 10. create summarized experiment
    se <- SummarizedExperiment(assays = SimpleList(scores = as.matrix(scores)), rowRanges = grNonOverlapping)
    
    #free up some memory
    remove(grList,grNonOverlapping)
    
    #--- 11. apply selection rules
    n = ncol(se)
    minSamples <- eval(parse(text=selectionRules))
    if(n >= minSamples){
      message("Applying selection criteria: minimum of ", minSamples," samples with a score-per-million of ", scorePerMillion)
      keep <- which(rowSums(assay(se) > scorePerMillion) >= minSamples)
      message("Final ", groups[i], " Peak Set ", length(keep))
      se <- se[keep,]
    }else{
      warning("group with less than min samples required by selectionRules (most likely 1 summit for group)")
    }
    
    #--- 12. re-normalize score-per-million
    mcols(se)$score <- edgeR::cpm(mcols(se)$score) #Re-scale for final comparison #10^6*mcols(se)$score/sum(mcols(se)$score)
    mcols(se)$name <- paste0(groups[i],"_",seq_len(nrow(se))) #add short name
    return(se)
    
  }) #END groupSEList lapply
  
  names(groupSEList) <- groups
  
  #Create GRanges List for all Groups
  groupGRlist <- lapply(groupSEList, rowRanges)
  #Get Non Overlapping Peaks for all groups
  message("Now Getting Non Overlapping Peak Set Between Groups...")
  
  #there may be some R version-specific errors from this command
  finalGR <- Reduce("c",groupGRlist) %>%
    convergeClusterGRanges(., by = "score", decreasing = T) %>%
    sortSeqlevels %>%
    sort
  
  message("Determined ", length(finalGR), " Peaks")
  
  return(list(groupSE = groupSEList, finalGR = finalGR))
  
}

columnOverlaps <- function(query, subject, colname = "score", decreasing = TRUE){
  #First get overlaps
  o <- findOverlaps(query, subject) %>% data.frame()
  #Then append information
  o$col <- mcols(subject)[[colname]][o[,2]]
  #Order it by the factor to rank
  o <- o[order(o$col, decreasing = decreasing),]
  #Deduplicate
  o <- o[!duplicated(o$queryHits),]
  #Initialize
  val <- rep(0, length(query))
  #Fill Values
  val[o[,1]] <- o$col
  return(val)
}

convergeClusterGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
  stopifnot(by %in% colnames(mcols(gr)))
  i = 0
  gr_initial <- gr #initial gr i
  if(verbose){
    message("Converging", appendLF = FALSE)
  }
  while(length(gr_initial) > 0){
    if(verbose){
      message(".", appendLF = FALSE)
    }
    i = i + 1
    gr_clustered <- clusterGRanges(gr = gr_initial, filter = TRUE, by = by, decreasing = decreasing, verbose = verbose) #initial gr i
    gr_initial <- subsetByOverlaps(gr_initial ,gr_clustered, invert=TRUE) #blacklist called cluster
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_clustered
    }else{
      gr_all <- c(gr_all, gr_clustered)
    } 
  }
  return(gr_all)
}

clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE, verbose = TRUE){
  #reduce first
  gr <- sort(sortSeqlevels(gr))
  r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
  o <- findOverlaps(gr,r)
  mcols(gr)$cluster <- subjectHits(o)
  if(verbose){
    message(sprintf("found %s overlaps...",length(gr) - max(subjectHits(o))))
  }
  #filter by
  if(filter){
    if(any(toupper(colnames(mcols(gr))) %in% toupper(by))){
      if(verbose){
        message(sprintf("filtering overlaps by %s...",by))
      }
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      grn <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(grn))
    }else{
      if(verbose){
        message(sprintf("filtering by order..."))
      }
      gr <- gr[!duplicated(mcols(gr)$cluster),]
    }
    mcols(gr)$cluster <- NULL
  }

  return(gr)
}

getChromSizes <- function (BSgenome, filter = TRUE) {
  stopifnot(!is.null(BSgenome))
  if (inherits(BSgenome, "BSgenome")) {
    BSgenome <- BSgenome
  }
  else if (is.character(BSgenome)) {
    BSgenome <- BSgenome::getBSgenome(BSgenome, masked = FALSE)
  }
  else {
    stop("Error cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }
  grCS <- GRanges(seqlevels(BSgenome), IRanges(1, seqlengths(BSgenome)))
  if (filter) {
    grCS <- keepFilteredChromosomes(grCS)
  }
  seqlengths(grCS) <- end(grCS)
  return(grCS)
}

keepFilteredChromosomes <- function (x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode = "coarse") {
  if (standard) {
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
  }
  seq_names <- seqlevels(x)
  chr_remove <- c()
  if (underscore) {
    chr_remove <- c(chr_remove, which(grepl("_", seq_names)))
  }
  chr_remove <- c(chr_remove, which(seq_names %in% remove))
  if (length(chr_remove) > 0) {
    chr_keep <- seq_names[-chr_remove]
  }
  else {
    chr_keep <- seq_names
  }
  seqlevels(x, pruning.mode = pruning.mode) <- chr_keep
  return(x)
}

removeCliffedPeaks <- function(gr, BSgenome) {
  chromSizes <- getChromSizes(BSgenome = BSgenome)
  return(subsetByOverlaps(gr, chromSizes, type = "within"))
}

removeSpecificChromosome <- function(gr, remove = "chrY"){
  seq_names <- seqlevels(gr)
  if(!(remove %in% seq_names)) {
    print(paste0("Warning - ",remove," does not exist in provided GRange"))
    return(gr)
  }
  else {
    gr <- gr[which(!(as.vector(seqnames(gr)) == remove))]
    return(gr)
  }
}

removeNbasePeaks <- function(gr, BSgenome, percent = 0){
  BSgenome <- BSgenome::getBSgenome(BSgenome, masked = FALSE)
  seqs <- getSeq(BSgenome, gr)
  nuc_freqs <- as.data.frame(letterFrequency(seqs, c("A","C","G","T","N")))
  return(gr[which(nuc_freqs$N <= percent)])
}
################################################################################################################
#MAIN
#--------------

# Read Inputs
#--------------
# metadata_df <- readMetadata(metadata)
metadata_df <- read.table(metadata, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
df <- data.frame(samples = rownames(metadata_df), groups = metadata_df$Group,
                 summits = paste0(macs2_dir,rownames(metadata_df),end))
blacklist <- import.bed(blacklist_bed)

peakCalls <- extendedPeakSet(df = df, 
                             BSgenome = genome,
                             selectionRules = selectionRules,
                             blacklist = blacklist,
                             scorePerMillion = spm,
                             extend = extend)

#Save output
message("Filter peaks and saving output...")
#Save Bed formatted files with the following suffix
end_bed <- ".fwp.filter.non_overlapping.bed"
#Save Summarized Experiment as rds files with the following suffix
end_rds <- ".fwp.filter.non_overlapping.rds"

for(i in seq_along(peakCalls[[1]])){
  se_i <- peakCalls[[1]][[i]]
  gr_i <- rowRanges(se_i)
  group_i <- names(peakCalls[[1]])[i]
  dir.create(file.path(out_dir, group_i),recursive=TRUE,showWarnings=FALSE)
  bed_out <- file.path(out_dir, group_i, paste0(group_i, end_bed))
  rds_out <- file.path(out_dir, group_i, paste0(group_i, end_rds))
  
  message(paste0("Processing ",group_i,"..."))
  orig_length <- length(gr_i)
  #remove peaks that extend past chromosome ends
  gr_i <- removeCliffedPeaks(gr_i, BSgenome = genome)
  cliffed <- orig_length - length(gr_i)
  message(paste0("Cliffed - ",cliffed))
  
  #remove peaks that map to chrY
  gr_i <- removeSpecificChromosome(gr_i, remove = "chrY")
  chrY <- orig_length - cliffed - length(gr_i)
  message(paste0("chrY - ",chrY))
  
  #remove peaks that span a region containing any unknown "N" bases in the genome sequence
  gr_i <- removeNbasePeaks(gr_i, BSgenome = genome, percent = 0)
  Nbase <- orig_length - cliffed - chrY - length(gr_i)
  message(paste0("Nbase - ",Nbase))
  
  message(paste0("Total Peaks - ",length(gr_i),"\n"))
  
  export.bed(gr_i, bed_out)
  saveRDS(se_i, rds_out)
}

message("Processing merged peak set...")
gr_merge <- peakCalls[[2]]

orig_length <- length(gr_merge)
#remove peaks that extend past chromosome ends
gr_merge <- removeCliffedPeaks(gr_merge, BSgenome = genome)
cliffed <- orig_length - length(gr_merge)

#remove peaks that map to chrY
gr_merge <- removeSpecificChromosome(gr_merge, remove = "chrY")
chrY <- orig_length - cliffed - length(gr_merge)

#remove peaks that span a region containing any unknown "N" bases in the genome sequence
gr_merge <- removeNbasePeaks(gr_merge, BSgenome = genome, percent = 0)
Nbase <- orig_length - cliffed - chrY - length(gr_merge)

message(paste0("Cliffed - ",cliffed))
message(paste0("chrY - ",chrY))
message(paste0("Nbase - ",Nbase))
message(paste0("Total Peaks - ",length(gr_merge),"\n"))

# 3. Save bed and RDS containing Genomic Ranges
export.bed(gr_merge, file.path(out_dir, paste0("All_Samples", end_bed)))
saveRDS(gr_merge, file.path(out_dir, paste0("All_Samples", end_rds)))

