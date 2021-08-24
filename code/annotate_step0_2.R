#################################################
#################################################
#################################################
##### function from diffloop
# adaptation: do not collapse anything
loopsMake.mango.null <- function(beddir, samples, mergegap = 500, ext = "all", pad = 0, reduce = TRUE) {
    
    library(foreach)
    library(GenomicRanges)
    library(utils)
    library(dplyr)
    library(readr)
    ext="all"
    loops <- setClass("loops", slots = c(anchors = "GRanges",
    interactions = "nim", counts = "nim",
    colData = "data.frame", rowData = "data.frame"))
    
    
    ct <- list(col_character(), col_integer(), col_integer(),
    col_character(), col_integer(), col_integer(), col_integer(),
    col_number())
    ct <- list(col_character(), col_integer(), col_integer(),
    col_character(), col_integer(), col_integer(), col_integer(),
    col_number())
    
    # Iterate through files to set up anchors
    anchorsraw <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
        bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
        plyr::rbind.fill(bt[, 1:3], setNames(bt[, 4:6], names(bt[, 1:3])))
    }
    
    if (reduce) {
        anchors <- reduce(padGRanges(makeGRangesFromDataFrame(do.call(rbind,
        anchorsraw), ignore.strand = TRUE, seqnames.field = "X1",
        start.field = "X2", end.field = "X3"), pad = pad), min.gapwidth = mergegap)
    }
    if (!reduce) {
        anchors <- padGRanges(makeGRangesFromDataFrame(do.call(rbind,
        anchorsraw), ignore.strand = TRUE, seqnames.field = "X1",
        start.field = "X2", end.field = "X3"), pad = pad)
    }
    
    
    petlist <- foreach(sample = samples) %do% {
        fullfile <- file.path(beddir, paste(sample, "interactions", ext, "mango", sep = "."))
        bt <- read_delim(fullfile, "\t", col_types = ct, col_names = FALSE)
        
        anchors_left = makeGRangesFromDataFrame(bt[, 1:3], ignore.strand = TRUE, seqnames.field = "X1",
        start.field = "X2", end.field = "X3")
        anchors_right = makeGRangesFromDataFrame(bt[, 4:6], ignore.strand = TRUE, seqnames.field = "X4",
        start.field = "X5", end.field = "X6")
        counts = bt[, 7]
        
        leftanchor = 1:(length(anchors)/2)
        rightanchor = ((length(anchors)/2)+1):length(anchors)
        df <- data.frame(left = leftanchor, right = rightanchor)
        g <- as.data.frame(dplyr::group_by(df, left, right))
        d <- cbind(g, counts)
        colnames(d) <- c("left", "right", "counts")
        dag <- d
        #dag <- aggregate(counts ~ left + right, FUN = sum, data=d)
        colnames(dag) <- c("left", "right", sample)
        dag
    }
    
    .full_join <- function(a, b) {
        as.data.frame(dplyr::full_join(a, b, by = c("left", "right")))
    }
    
    # Map Counts
    pets <- Reduce(.full_join, petlist)
    iraw <- pets[, c("left", "right")]
    iraw <- t(apply(iraw, 1, function(x) {
        if (x[1] < x[2]) {
            x
        } else {
            x[c(2, 1)]
        }
    }))
    interactions <- iraw[order(iraw[, 1], iraw[, 2]), ]
    colnames(interactions) <- c("left", "right")
    counts <- as.matrix(pets[, -c(1:2)])[order(iraw[, 1], iraw[, 2]), ]
    counts[is.na(counts)] <- 0
    counts <- as.matrix(counts, ncol = length(samples))
    colnames(counts) <- samples
    
    # Initialize rowData slot (with loop widths)
    w <- (start(anchors[interactions[, 2]]) + end(anchors[interactions[, 2]]))/2 -
    (start(anchors[interactions[, 1]]) + end(anchors[interactions[, 1]]))/2
    w[w < 0] <- 0
    rowData <- as.data.frame(as.integer(w))
    colnames(rowData) <- c("loopWidth")
    
    # Remove 'chr' from anchors
    seqlevels(anchors) <- gsub("^chr(.*)$", "\\1", seqlevels(anchors))
    
    # Remove rownames from matrices
    row.names(interactions) <- NULL
    row.names(counts) <- NULL
    
    # Initialize colData slot
    groups <- rep("group1", length(samples))
    if(length(samples) == 1){
        sizeFactor <- 1
    } else {
        lc <- log2(counts)
        keep <- rowSums(counts > 0) == ncol(lc)
        lc <- lc[keep, ]
        target <- 2^rowMeans(lc)
        sizeFactor <- colMedians(sweep(2^lc, 1, target, FUN = "/"), na.rm = TRUE)
    }
    dfcd <- data.frame(sizeFactor, groups)
    rownames(dfcd) <- samples
    
    # Create loops object
    dlo <- loops()
    slot(dlo, "anchors", check = TRUE) <- anchors
    slot(dlo, "interactions", check = TRUE) <- interactions
    slot(dlo, "counts", check = TRUE) <- counts
    slot(dlo, "colData", check = TRUE) <- dfcd
    slot(dlo, "rowData", check = TRUE) <- rowData
    
    return(dlo)
}


annotateAnchors2 <- function(loops=all, features=eqtl_sig_gr, featureName="EQTL.TCGA", featureToAdd = "eQTLs", maxgap=0) {
    if (!(featureToAdd %in% names(mcols(features)))) stop("Do not have featureToAdd in data")
    lto = loops
    lto.df <- summary(lto)
    Lanchors <- GRanges(lto.df$chr_1, IRanges(lto.df$start_1, lto.df$end_1))
    Ranchors <- GRanges(lto.df$chr_2, IRanges(lto.df$start_2, lto.df$end_2))
    
    # Determine if right anchor overlaps GWAS region
    Rhits.p <- suppressWarnings(findOverlaps(features, Ranchors,
    maxgap = maxgap))
    Rvalues.p <- rep(FALSE, dim(lto.df)[1])
    Rvalues.p[unique(subjectHits(Rhits.p))] <- TRUE
    
    # Determine if left anchor overlaps GWAS region
    Lhits.p <- suppressWarnings(findOverlaps(features, Lanchors,
    maxgap = maxgap))
    Lvalues.p <- rep(FALSE, dim(lto.df)[1])
    Lvalues.p[unique(subjectHits(Lhits.p))] <- TRUE
    
    # Aggregate TSS
    Rtss <- data.frame(Rhits.p)
    Rtss <- cbind(Rtss, mcols(features[Rtss$queryHits]))
    Ltss <- data.frame(Lhits.p)
    Ltss <- cbind(Ltss, mcols(features[Ltss$queryHits]))
    # collapse
    r_ttss <- unique(Rtss[,c(2,which(names(Rtss) == featureToAdd))])
    r_tss <- aggregate(r_ttss[,featureToAdd]~subjectHits,paste,collapse=",",data=r_ttss)
    names(r_tss)[2] = featureToAdd
    
    l_ttss <- unique(Ltss[,c(2,which(names(Ltss) == featureToAdd))])
    l_tss <- aggregate(l_ttss[,featureToAdd]~subjectHits,paste,collapse=",",data=l_ttss)
    names(l_tss)[2] = featureToAdd
    
    anchor.1 <- rep("none", dim(lto)[2]) # record left is 1
    anchor.2 <- rep("none", dim(lto)[2]) # record right as 2
    
    anchor.1[l_tss$subjectHits] <- l_tss[,featureToAdd]
    anchor.2[r_tss$subjectHits] <- r_tss[,featureToAdd]
    
    lto@rowData$anchor.1 <- anchor.1
    lto@rowData$anchor.2 <- anchor.2
    
    colnames(lto@rowData)[which(colnames(lto@rowData) == "anchor.1")] = paste(featureName, "_1", sep="")
    colnames(lto@rowData)[which(colnames(lto@rowData) == "anchor.2")] = paste(featureName, "_2", sep="")
    
    return(lto)
}


annotateAnchorType2 <- function(
loops,
h3k27peakid = "h3k27peak",
gene = "HiChIP.Gene",
anchor=1) {
    
    df_anno <- summary(loops)
    loop.types <- 1:nrow(df_anno)
    description <- rep("none", length(loop.types))
    
    h3k27peakid = paste(h3k27peakid, anchor, sep="_")
    gene = paste(gene, anchor, sep="_")
    
    if ( !h3k27peakid %in% names(df_anno) ) stop("h3k27peakid column not in data")
    if ( !gene %in% names(df_anno) ) stop("gene column not in data")
    
    p1 = (
    df_anno[,h3k27peakid] != "none" &
    df_anno[,gene] !="none"
    )
    p2 = (
    df_anno[,h3k27peakid] == "none" &
    df_anno[,gene] !="none"
    )
    if ( length(intersect(which(p1),which(p2))) > 0) stop("Promoter definitions overlapping")
    p = p1 | p2
    #length(which(p))
    
    e = (
    df_anno[,h3k27peakid] != "none" &
    df_anno[,gene] =="none"
    )
    
    #if ( length(intersect(which(e1),which(e2)))  > 0) stop("Enhancer definitions overlapping")
    #e = e1 | e2
    #length(which(e))
    
    if ( length(intersect(which(p),which(e))) > 0) stop()
    
    other = (
    df_anno[,h3k27peakid] == "none" &
    df_anno[,gene] =="none"
    )
    
    if ( length(intersect(which(p),which(other))) > 0) stop("other definitions overlapping")
    if ( length(intersect(which(e),which(other))) > 0) stop("other definitions overlapping")
    
    
    # sum all: should be equal to number of loops
    description[p] <- "P"
    description[e] <- "E"
    description[other] <- "O"
    
    if ( "none" %in% names(table(description)) )
    print("Unclassified cases (none): look what these are!")
    loops@rowData$type <- description
    names(loops@rowData)[ncol(loops@rowData)] = paste("anchor.type_", anchor, sep="")
    #names(loops@rowData)[ncol(loops@rowData)] = paste("anchor.type.noamb_", anchor, sep="")
    return(loops)
    
}

# if summarizeby, return other columns grouped by loopid or anchorid
# if not summarizeby, return the original dataset with added true/false for matching
find_match = function(match1="HiChIP.Gene_1", match2="Gene.TWAS.anchor_1", df, idcol, split=FALSE, split_match = 2, summarizeby = TRUE, equal = TRUE) {
    if (!(match1 %in% colnames(df) )) stop("Column ", match1, " is missing")
    if (!(match2 %in% colnames(df) )) stop("Column ", match2, " is missing")
    if (!(idcol %in% colnames(df) )) stop("Column ", idcol, " is missing")
    if (any(grepl(",", df[,idcol]))) stop("idcol has already been grouped")
    
    
    t = df
    t = t[t[,match1] !="none" | t[,match2]!="none",] # to help with speed
    t = separate_rows(t, match1, sep = ",")
    t = separate_rows(t, match2, sep = ",")
    if (split & split_match == 1) {
        t[,match1] = as.character(lapply(strsplit(as.character(t[,match1]), ":", fixed=TRUE), "[", 2))
    }
    if (split & split_match == 2) {
        t[,match2] = as.character(lapply(strsplit(as.character(t[,match2]), ":", fixed=TRUE), "[", 2))
    }
    
    t = t[t[,match1] !="none",]
    t = t[t[,match2] !="none",]
    t = t[!is.na(t[,match1]),]
    t = t[!is.na(t[,match1]),]
    
    if (equal) {
        t =t[t[,match1] == t[,match2],]
    } else {
        t =t[t[,match1] != t[,match2],]
    }
    t = unique(t)
    if (!summarizeby) {
        id_interest = unique(t[,idcol])
        df$temp = ifelse(df[,idcol] %in% id_interest, TRUE, FALSE)
        names(df)[ncol(df)] = paste(match1, match2, sep="_")
        return(df)
    } else {
        # collapse by id
        # print(cbind(t[,match1], t[,match2]))
        library(plyr)
        # tss <- aggregate(.~t[,idcol], data = t, function(x) paste0(unique(x), collapse=","))
        # https://stackoverflow.com/questions/53341434/group-by-a-column-and-collapse-all-other-columns-without-na
        tss = data.frame (t %>%
        mutate_if(is.factor, as.character)  %>%
        mutate_all(~replace(., .=='NA', NA)) %>%
        group_by(t[,idcol]) %>%
        summarize_all(~paste(unique(na.omit(.)), collapse = ',')) )
        #summarise_all(funs(unique(paste(ifelse(is.na(.), "null", .), collapse = ",")))) )
        return(tss)
    }
}

#################################################
#################################################
#################################################

# Annotate enhancer, promoter, anchor type, loop type, and save the data for further annotations
# annotate_step0.R and annotate_step2.R

pad_anchors = 5000
P2PBckgr = "P2PBckgr_0"

# source("/u/project/pasaniuc/pasaniucdata/DATA/HiChIPprostate/HiCpro/scripts/annotate/myfunctions.R")
# source("myfunctions.R")
library(diffloop)
library(data.table)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(dplyr)


dir = "/u/project/pasaniuc/pasaniucdata/DATA/HiChIPprostate/"
setwd(dir)
outdir_main = "/u/project/pasaniuc/pasaniucdata/DATA/HiChIPprostate/manuscript1/AnnotatedData/"
outdir = paste(outdir_main, P2PBckgr, "/", sep="")
	if (!file.exists (outdir)) dir.create(outdir)

## 1. Import genes for promoters: Genes file does not change
genes_file= "EQTL/hg19/knownGeneRefSeqV19.gtf.gz"
genes = read.delim(genes_file, header=T, stringsAsFactors=F)
genes = genes[,c("chrom", "strand", "txStart", "txEnd", "name2"),]
chroms = c("chrX", "chrY", paste("chr", 1:22, sep=""))
genes = genes[genes$chrom %in% chroms,]
length(unique(genes$name2))
# For each gene, take longest transcript: should I *NOT* do this? Try to annotate without this, and name it "_allrefseqtranscripts"
# ** #
choose_longest = TRUE
if (choose_longest) {
genes$length = abs(genes$txStart - genes$txEnd)
genes = data.frame(genes %>% group_by(name2) %>% top_n(1, length))
genes = unique(genes)
genes = genes[!(duplicated(genes$name2) | duplicated(genes$name2, fromLast = TRUE)), ]
}
# ** #
names(genes)[which(names(genes)=="name2")] = "gene"

# define the TSS as the start position for each entry that is on the + strand and as the end position of every entry that is on the - strand
pad = 500
# or can use GenomicRanges function:
genes = data.frame(promoters(makeGRangesFromDataFrame(genes, ignore.strand = FALSE, seqnames.field = "chrom", start.field = "txStart", end.field = "txEnd", strand.field="strand", keep.extra.columns=T), upstream=pad, downstream=pad)) #, use.names=TRUE))
colnames(genes)[1:3] = c("chrom", "TSS_start", "TSS_end")

promoters = makeGRangesFromDataFrame(genes, seqnames.field="chrom", start.field="TSS_start", end.field="TSS_end", keep.extra.columns = T)
promoters = rmchr(promoters)
# Need positions of genes to overlap eQTL and TWAS data
genepos = data.frame(promoters)
genepos = genepos[,c(1,2,3,6)]
if (choose_longest) genepos_file = "EQTL/genepos.txt"
if (!choose_longest) genepos_file = "EQTL/genepos_allrefseqtranscripts.txt"
print(genepos_file)
# write.table(genepos, file = genepos_file, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")

## 2. Import chipseq for enhancers: this file will change depending on cell type
pad = 0 #1000 # 2500
# 2.1 import both broad and narrow peaks for enhancers
peaks_file = "other_peaks/H3K27Ac_LNCaP_regular_media/broad_narrow.bed"
h3k27_peaks = read.table(peaks_file, header=F, stringsAsFactors=F)
h3k27_peaks$h3k27peakid = paste(h3k27_peaks$V1, h3k27_peaks$V2, h3k27_peaks$V3, sep="_")
h3k27_peaks_gr = makeGRangesFromDataFrame(h3k27_peaks, seqnames.field="V1", start.field="V2", end.field="V3", keep.extra.columns = T)
h3k27_peaks_gr = rmchr(h3k27_peaks_gr)

enhancers = h3k27_peaks
if (pad !=0) {
enhancers$V2 = enhancers$V2 - pad
enhancers$V3 = enhancers$V3 + pad -1
}
enhancers = makeGRangesFromDataFrame(enhancers, seqnames.field="V1", start.field="V2", end.field="V3", keep.extra.columns = T)
enhancers = rmchr(enhancers)
# 2.2 Take promoters out of enhancers:
enhancers_lncap = enhancers[!(enhancers %over% promoters)]

# 3. Import data and annotate
############################################################
################## START ANNOTATION ########################
############################################################
# 1)

# if back1:
if (P2PBckgr == "P2PBckgr_1") {
filename = "FitHiChIP.interactions_FitHiC_Q0.01.bed"
# LNCaP
bed_dir_lncap= "replicates/combine_replicates_and_recall_loops/FitHiChIP/Rep1_Rep2_Rep3_Rep4_Rep5/FitHiChIP_Peak2ALL_b5000_L5000_U3000000/P2PBckgr_1/Coverage_Bias/FitHiC_BiasCorr/"
}

# if back0: 
if (P2PBckgr == "P2PBckgr_0") {
filename = "FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bed"

# LNCaP
bed_dir_lncap = "replicates/combine_replicates_and_recall_loops/FitHiChIP/Rep1_Rep2_Rep3_Rep4_Rep5/FitHiChIP_Peak2ALL_b5000_L5000_U3000000/P2PBckgr_0/Coverage_Bias/FitHiC_BiasCorr/Merge_Nearby_Interactions/"

}

message("Reading from bed dir ", bed_dir_lncap)
# either apply make_hichipper.R or this:
bedpe_files = paste(bed_dir_lncap, filename, ".gz", sep="")
out_file = paste(bed_dir_lncap, filename, ".interactions.all.mango", sep="")

if (!file.exists(out_file)) {
message("Creating hichipper input...")
        x = read.table(bedpe_files, header=T, stringsAsFactors = F)

        x_new = x[,c("chr1", "s1", "e1", "chr2", "s2", "e2", "cc", "Q.Value_Bias")]
        x_new$cc = as.integer(x_new$cc)
        x_new$Q.Value_Bias = as.numeric(x_new$Q.Value_Bias)


        message("File saved in: ", out_file)
        # write.table(x_new, file = out_file, row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")

}
all_lncap = loopsMake.mango.null(bed_dir_lncap, filename, mergegap=0, pad=pad_anchors, reduce=FALSE)
all_lncap <- rmchr(all_lncap)

all_lncap = annotateAnchors2(loops=all_lncap, features=h3k27_peaks_gr, featureName="h3k27peak", featureToAdd = "h3k27peakid", maxgap=0)
all_lncap = annotateAnchors2(loops=all_lncap, features=enhancers_lncap, featureName="Enhancerpeakid", featureToAdd = "h3k27peakid", maxgap=0)

all_lncap = annotateAnchors2(loops=all_lncap, features=promoters, featureName="HiChIP.Gene", featureToAdd = "gene", maxgap=0)

all_lncap = annotateAnchorType2(loops=all_lncap, anchor=1) # no ambiguous
all_lncap = annotateAnchorType2(loops=all_lncap, anchor=2) # ambiguous considered E

df_anno = summary(all_lncap)
df_anno$loopid = paste(df_anno$chr_1, df_anno$start_1, df_anno$end_1, df_anno$chr_2, df_anno$start_2, df_anno$end_2, sep="_")

df_anno$loop.type_details <- apply(df_anno[,c("anchor.type_1", "anchor.type_2")], 1, function(x) paste0(sort(x), collapse = "-"))
df_anno$loop.type_details_direction = paste(df_anno$anchor.type_1, df_anno$anchor.type_2, sep="-")

with_ambiguous = FALSE
if (with_ambiguous) {
	# if with_ambiguous = TRUE, loop.type_details is 9 categories, now collapse into 4 categories
	e.loops = df_anno$loop.type_details %in% c("A-E", "E-E", "E-O")
	p.loops = df_anno$loop.type_details %in% c("A-P", "O-P", "P-P")
	ep.loops = df_anno$loop.type_details %in% c("E-P")
	a.loops = df_anno$loop.type_details %in% c("A-A", "A-O")
	df_anno$loop.type = "none"
	df_anno$loop.type[e.loops] <- "E"
	df_anno$loop.type[p.loops] <- "P"
	df_anno$loop.type[ep.loops] <- "E-P"
	df_anno$loop.type[a.loops] <- "A"
}
if (!with_ambiguous) {
        e.loops = df_anno$loop.type_details %in% c("E-E", "E-O")
        p.loops = df_anno$loop.type_details %in% c("O-P", "P-P")
        ep.loops = df_anno$loop.type_details %in% c("E-P")
        df_anno$loop.type = "none"
        df_anno$loop.type[e.loops] <- "E"
        df_anno$loop.type[p.loops] <- "P"
        df_anno$loop.type[ep.loops] <- "E-P"
}
df_anno$loop.type_details = gsub("O-P", "P-O", df_anno$loop.type_details)

# rename PET counts
names(df_anno)[which(names(df_anno) %in% c("FitHiChIP.interactions_FitHiC_Q0.01.bed", "FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bed"))]="PET_Q0.01"
df_anno_lncap = df_anno
rm("df_anno")


# Save by loop
#outFile_loops = paste("AnnotatedData/", gsub("/", "_", bed_dir_lncap), "loops.txt", sep="")
outFile_loops = paste(outdir, "LNCAP", "_Q0.01_loops.txt", sep="")
if (pad_anchors != 0 ) outFile_loops = paste(outdir, "LNCAP", "_Q0.01_pad", pad_anchors, "_loops.txt", sep="")
if (!choose_longest) outFile_loops = gsub("_loops.txt", "_loops_allrefseqtx.txt", outFile_loops)

message("Loop data saved in ", outFile_loops)
# write.table(df_anno_lncap, file = outFile_loops, quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)

# Save bed file of unique anchors
anchors = rbind.data.frame(
		data.frame(chr=paste("chr", df_anno_lncap$chr_1, sep=""), start=df_anno_lncap$start_1, stop=df_anno_lncap$end_1), 
		data.frame(chr=paste("chr", df_anno_lncap$chr_2, sep=""), start=df_anno_lncap$start_2, stop=df_anno_lncap$end_2))
anchors = unique(anchors)

outFile_anchors = paste(outdir, "LNCaP_", "anchors.bed", sep="")
if (pad_anchors != 0 ) outFile_anchors = paste(outdir, "LNCaP", "_pad", pad_anchors, "_anchors.bed", sep="")
message("Anchors data saved in ", outFile_anchors)


# write.table(anchors, outFile_anchors, quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)
# Save reduced version too
# write.table(data.frame(reduce(makeGRangesFromDataFrame(anchors))), paste(outdir, "LNCaP_anchors_reduced.bed", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE)


#######################
message("################### SAVING ANNOTATED DATA ###################")
# Save by loop
#outFile_loops = paste("AnnotatedData/", gsub("/", "_", bed_dir), "loops.txt", sep="")
#message("Loop data saved in ", outFile_loops)
#write.table(df_anno, file = outFile_loops, quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)

d = all_lncap
d_df = df_anno_lncap
d_nm = "LNCaP"
############## ANNOTATE
outdir = "/u/project/pasaniuc/pasaniucdata/DATA/HiChIPprostate/manuscript1/AnnotatedData/P2PBckgr_0/"
outdir_annot = "/u/project/pasaniuc/pasaniucdata/DATA/HiChIPprostate/manuscript1/AnnotatedData/data_for_annotation/"

outFile_loops = paste(outdir, "AnnotatedData_", d_nm, "_loops.txt", sep="")
outFile_genes = paste(outdir, "AnnotatedData_", d_nm, "_genes.txt", sep="")
outFile_enhancers = paste(outdir, "AnnotatedData_", d_nm, "_enhancers.txt", sep="")

if (pad_anchors != 0 ) {
    outFile_loops = paste(outdir, "AnnotatedData_", d_nm, "_pad", pad_anchors, "_loops.txt", sep="")
    outFile_genes = paste(outdir, "AnnotatedData_", d_nm, "_pad", pad_anchors, "_genes.txt", sep="")
    outFile_enhancers = paste(outdir, "AnnotatedData_", d_nm, "_pad", pad_anchors, "_enhancers.txt", sep="")
}

if (acetylated) {
    outFile_loops = gsub("_loops.txt", "_loops_acetylated.txt", outFile_loops)
    outFile_genes = gsub("_genes.txt", "_genes_acetylated.txt", outFile_genes)
    outFile_enhancers = gsub("_enhancers.txt", "_enhancers_acetylated.txt", outFile_enhancers)
}

if (!choose_longest) {
    outFile_loops = gsub("_loops", "_loops_allrefseqtx", outFile_loops)
    outFile_genes = gsub("_genes", "_genes_allrefseqtx", outFile_genes)
    outFile_enhancers = gsub("_enhancers", "_enhancers_allrefseqtx", outFile_enhancers)
}

if (length(d)==0) {
    message("Dataset: does not exist in annotated data. Run step 0!")
    stop()
}
#######################
# 1) GWAS
# a) GWAS filtered by MAF and Flipped
#######################
options(scipen = 999)

gwas_file = paste(outdir_annot, "PrCa_Flipped_MAF0.001.txt.gz", sep="")
gwas= fread(paste("gunzip -c ", gwas_file, sep=""), header=T, stringsAsFactors=F, data.table=F)
gwas_sig = gwas[gwas$Pvalue<5*10^-8,]

# If get this error:
# no rows to aggregate
# remember to "rmchr"
gwas_sig.gr = rmchr(makeGRangesFromDataFrame(gwas_sig, seqnames.field="CHR", start.field="start", end.field="end", keep.extra.columns = T))

# Overlap with loops
t = annotateAnchors2(loops=d, features=gwas_sig.gr, featureName="GWAS_SNP", featureToAdd = "ch_pos", maxgap=pad)
t = annotateAnchors2(loops=t, features=gwas_sig.gr, featureName="GWAS_Pvalue", featureToAdd = "Pvalue", maxgap=pad)

#######################
# 1) GWAS Fine-mapping
# a) Paintor results
#######################
paintor_file = paste(outdir_annot, "paintor_1causals.txt.gz", sep="")
paintor = read.table(paintor_file, header=T, stringsAsFactors = F)
paintor_gr = rmchr(makeGRangesFromDataFrame(paintor, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))

t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_causal", featureToAdd = "ID", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_causal_rsid", featureToAdd = "rsid", maxgap=pad)

t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_Posterior_Prob", featureToAdd = "Posterior_Prob", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_region", featureToAdd = "region", maxgap=pad)

#t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_index_chrpos", featureToAdd = "index_chrpos", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor_gr, featureName="Paintor_index_rsid", featureToAdd = "index_rsid", maxgap=pad)

# Add numb SNPs in the causal 95% set
#t = annotateAnchors2(loops=d, features=paintor_gr, featureName="Paintor_Ncredset", featureToAdd = "n", maxgap=pad)

# Add paintor causals passing PP>=0.1 and PP>=0.2
paintor0.2_gr = rmchr(makeGRangesFromDataFrame(paintor[paintor$Posterior_Prob >= 0.2,], seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))

t = annotateAnchors2(loops=t, features=paintor0.2_gr, featureName="Paintor_causal0.2", featureToAdd = "ID", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor0.2_gr, featureName="Paintor_causal0.2_rsid", featureToAdd = "rsid", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor0.2_gr, featureName="Paintor_Posterior_Prob0.2", featureToAdd = "Posterior_Prob", maxgap=pad)

paintor0.1_gr = rmchr(makeGRangesFromDataFrame(paintor[paintor$Posterior_Prob >= 0.1,], seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))

t = annotateAnchors2(loops=t, features=paintor0.1_gr, featureName="Paintor_causal0.1", featureToAdd = "ID", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor0.1_gr, featureName="Paintor_causal0.1_rsid", featureToAdd = "rsid", maxgap=pad)
t = annotateAnchors2(loops=t, features=paintor0.1_gr, featureName="Paintor_Posterior_Prob0.1", featureToAdd = "Posterior_Prob", maxgap=pad)

#######################
# 1) GWAS
# a) Dadaev results
#######################

dad_file = paste(outdir_annot, "dadaev.txt.gz", sep="")
dad = read.table(dad_file, header=T, stringsAsFactors = F)
dad_gr = rmchr(makeGRangesFromDataFrame(dad, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))
# Overlap with loops
t = annotateAnchors2(loops=t, features=dad_gr, featureName="Dadaev_causal", featureToAdd = "ch_pos", maxgap=0)
#PP_best_tag
t = annotateAnchors2(loops=t, features=dad_gr, featureName="Dadaev_PP_best_tag", featureToAdd = "PP_best_tag_rsid", maxgap=0)

#######################
# 2) Thib # computed, taking the top SNP (best pval) per gene
#######################
computed =TRUE
if (computed) {
    #eqtl_sig_file = paste(outdir_annot, "Thib_eqtl_sig_fdr0.05.txt.gz", sep="")
    #eqtl_sig_file = paste(outdir_annot, "Thib_eqtl_sig_bonf.txt.gz", sep="")
    # This is the top SNP per gene
    eqtl_sig_file = paste(outdir_annot, "Thib_eqtl_sig_bonf_bestSNP.txt.gz", sep="")
    #bed_file = paste(outdir_bed, "Thib_eGenes", sep="")
    
    eqtl_sig = read.table(eqtl_sig_file, header=T, stringsAsFactors = F)
    eqtl_sig$CHR = gsub("chr", "", eqtl_sig$CHR)
    eqtl_sig$chr_pos = gsub("chr", "", eqtl_sig$chr_pos)
    eqtl_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))
    egenes_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig[,c("CHR_gene", "start_gene", "stop_gene", "egenes")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))
    
    t = annotateAnchors2(loops=t, features=eqtl_sig_gr, featureName="EQTL.Thibo.PRAD.anchor", featureToAdd = "EQTL.eGENE", maxgap=pad)
    
    # Annotate if there is an eGene TSS overlapping an anchor
    # Use promoters position
    t = annotateAnchors2(loops=t, features=egenes_sig_gr, featureName="eGene.Thibo.PRAD.anchor", featureToAdd = "egenes", maxgap=pad)
    
    message("File: ", eqtl_sig_file, "; nSNPS:",  nrow(eqtl_sig), " ; nGenes:", length(unique(eqtl_sig$egenes)))
    rm("eqtl_sig_file")
    #rm("bed_file")
}
#######################
# 3) TCGA:
# a) reported: difference in duplicate removal
# b) computed myself: # computed, taking the top SNP (best pval) per gene
#######################
reported = FALSE
if (reported) {
    eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_reported.txt.gz", sep="")
    #bed_file = paste(outdir_bed, "TCGA_eGenes_reported", sep="")
} else {
    #eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_sig_fdr0.05.txt.gz", sep="")
    #eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_sig_bonf.txt.gz", sep="")
    eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_sig_bonf_bestSNP.txt.gz", sep="")
}
eqtl_sig = read.table(eqtl_sig_file, header=T, stringsAsFactors = F)
eqtl_sig$CHR = gsub("chr", "", eqtl_sig$CHR)
eqtl_sig$chr_pos = gsub("chr", "", eqtl_sig$chr_pos)
eqtl_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))
egenes_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig[,c("CHR_gene", "start_gene", "stop_gene", "egenes")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))

# This is the top SNP per gene
t = annotateAnchors2(loops=t, features=eqtl_sig_gr, featureName="EQTL.TCGA.PRAD.anchor", featureToAdd = "EQTL.eGENE", maxgap=pad)
# Annotate if there is an eGene TSS overlapping an anchor
# Use promoters position
t = annotateAnchors2(loops=t, features=egenes_sig_gr, featureName="eGene.TCGA.PRAD.anchor", featureToAdd = "egenes", maxgap=pad)

message("File: ", eqtl_sig_file, "; nSNPS:",  nrow(eqtl_sig), " ; nGenes:", length(unique(eqtl_sig$egenes)))
rm("eqtl_sig_file")

#######################
# 4) Reproted
# a) reported GTEx
#######################
eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_reported.txt.gz", sep="")
eqtl_sig = read.table(eqtl_sig_file, header=T, stringsAsFactors = F)
eqtl_sig$CHR = gsub("chr", "", eqtl_sig$CHR)
eqtl_sig$chr_pos = gsub("chr", "", eqtl_sig$chr_pos)
eqtl_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))
egenes_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig[,c("CHR_gene", "start_gene", "stop_gene", "egenes")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))

# This is the top SNP per gene
t = annotateAnchors2(loops=t, features=eqtl_sig_gr, featureName="EQTL.TCGA.PRAD.reported.anchor", featureToAdd = "EQTL.eGENE", maxgap=pad)
# Annotate if there is an eGene TSS overlapping an anchor
# Use promoters position
t = annotateAnchors2(loops=t, features=egenes_sig_gr, featureName="eGene.TCGA.PRAD.reported.anchor", featureToAdd = "egenes", maxgap=pad)
message("File: ", eqtl_sig_file, "; nSNPS:",  nrow(eqtl_sig), " ; nGenes:", length(unique(eqtl_sig$egenes)))
rm("eqtl_sig_file")

#eqtl_sig_file = paste(outdir_annot, "TCGA_eqtl_sig_bonf_bestSNP.txt.gz", sep="")
eqtl_sig_file = paste(outdir_annot, "GTEx_eqtl_reported.txt.gz", sep="")
eqtl_sig = read.table(eqtl_sig_file, header=T, stringsAsFactors = F)
eqtl_sig$CHR = gsub("chr", "", eqtl_sig$CHR)
eqtl_sig$chr_pos = gsub("chr", "", eqtl_sig$chr_pos)

eqtl_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig, seqnames.field="CHR", start.field="start", end.field="stop", keep.extra.columns = T))
egenes_sig_gr = rmchr(makeGRangesFromDataFrame(eqtl_sig[,c("CHR_gene", "start_gene", "stop_gene", "egenes")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))
# This is the top SNP per gene
t = annotateAnchors2(loops=t, features=eqtl_sig_gr, featureName="EQTL.GTEx.Prost.anchor", featureToAdd = "EQTL.eGENE", maxgap=pad)
# Annotate if there is an eGene TSS overlapping an anchor
# Use promoters position
t = annotateAnchors2(loops=t, features=egenes_sig_gr, featureName="eGene.GTEx.Prost.anchor", featureToAdd = "egenes", maxgap=pad)

message("File: ", eqtl_sig_file, "; nSNPS:",  nrow(eqtl_sig), " ; nGenes:", length(unique(eqtl_sig$egenes)))
rm("eqtl_sig_file")


#######################
# 5) TWAS
#######################
twas_file = paste(outdir_annot, "TWAS_sig.txt.gz", sep="")

twas = read.table(twas_file, header=T, stringsAsFactors = F)

twasgenes_gr = rmchr(makeGRangesFromDataFrame(twas[,c("CHR_gene", "start_gene", "stop_gene", "Gene", "id")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))

t = annotateAnchors2(loops=t, features=twasgenes_gr, featureName="Gene.TWAS.anchor", featureToAdd = "Gene", maxgap=0)

#######################
# 6) COLOC file
#######################
coloc_file_tcga = paste(outdir_annot, "COLOC_TCGA.txt.gz", sep="") # COLOC was run only on the 147 regions
coloc_147_tcga = read.table(coloc_file_tcga, header=T, stringsAsFactors=F)
coloc_147_tcga = coloc_147_tcga[coloc_147_tcga$pp4 >= 0.75,]
coloc_147_tcga_gr = rmchr(makeGRangesFromDataFrame(coloc_147_tcga[,c("gwas_gene", "gwas_region_chr", "gwas_region_start", "gwas_region_stop", "geneid",  "pp4")], seqnames.field="gwas_region_chr", start.field="gwas_region_start", end.field="gwas_region_stop",keep.extra.columns = T))
coloc_147_tcga_genes.gr = rmchr(makeGRangesFromDataFrame(coloc_147_tcga[,c("CHR_gene", "start_gene",  "stop_gene", "geneid", "gwas_region_chr", "gwas_region_start", "gwas_region_stop", "pp4")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))

coloc_file_thib = paste(outdir_annot, "COLOC_Thib.txt.gz", sep="")
coloc_147_thib = read.table(coloc_file_thib, header=T, stringsAsFactors=F)
coloc_147_thib = coloc_147_thib[coloc_147_thib$pp4 >= 0.75,]
coloc_147_thib_gr = rmchr(makeGRangesFromDataFrame(coloc_147_thib[,c("gwas_gene", "gwas_region_chr", "gwas_region_start", "gwas_region_stop", "geneid",  "pp4")], seqnames.field="gwas_region_chr", start.field="gwas_region_start", end.field="gwas_region_stop",keep.extra.columns = T))
coloc_147_thib_genes.gr = rmchr(makeGRangesFromDataFrame(coloc_147_thib[,c("CHR_gene", "start_gene",  "stop_gene", "geneid", "gwas_region_chr", "gwas_region_start", "gwas_region_stop", "pp4")], seqnames.field="CHR_gene", start.field="start_gene", end.field="stop_gene", keep.extra.columns = T))

# annotate anchor with GWAS region used
t = annotateAnchors2(loops=t, features=coloc_147_tcga_genes.gr, featureName="Gene.COLOC_tcga.anchor", featureToAdd = "geneid", maxgap=0)
t = annotateAnchors2(loops=t, features=coloc_147_thib_genes.gr, featureName="Gene.COLOC_thib.anchor", featureToAdd = "geneid", maxgap=0)


##################

df_anno = summary(t)
dupl_cols = colnames(df_anno)[duplicated(colnames(df_anno))]
if (length(dupl_cols) > 0 ) stop("Data has been annotated twice: restart from step 0.")

message("Done annotating dataset, creating final columns for matches...")
#######################
#######################

df_anno$loopid = paste(df_anno$chr_1, df_anno$start_1, df_anno$end_1, df_anno$chr_2, df_anno$start_2, df_anno$end_2, sep="_")

#1.a. is the hichip gene a target gene for an eqtl overlapping the *opposite* anchor?
t1 = find_match(match1="HiChIP.Gene_1", match2="EQTL.Thibo.PRAD.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="EQTL.Thibo.PRAD.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]] ))
#id_interest = unique(c(t1$anchorid, t2$anchorid))
df_anno$EQTLThib_and_HiChIP_opposite = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)

# 1.c same for TCGA
t1 = find_match(match1="HiChIP.Gene_1", match2="EQTL.TCGA.PRAD.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="EQTL.TCGA.PRAD.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]] ))
#id_interest = unique(c(t1$anchorid, t2$anchorid))
df_anno$EQTLTCGA_and_HiChIP_opposite = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)

# 3.a is the hichip gene a colocalized gene in general?
t1 = find_match(match1="HiChIP.Gene_1", match2="Gene.COLOC_tcga.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="Gene.COLOC_tcga.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t3 = find_match(match1="HiChIP.Gene_1", match2="Gene.COLOC_tcga.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t4 = find_match(match1="HiChIP.Gene_2", match2="Gene.COLOC_tcga.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
df_anno$COLOCTCGA_and_HiChIP = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)

t1 = find_match(match1="HiChIP.Gene_1", match2="Gene.COLOC_thib.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="Gene.COLOC_thib.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t3 = find_match(match1="HiChIP.Gene_1", match2="Gene.COLOC_thib.anchor_1", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
t4 = find_match(match1="HiChIP.Gene_2", match2="Gene.COLOC_thib.anchor_2", df_anno, idcol="loopid", split=TRUE, split_match = 2, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
df_anno$COLOCThib_and_HiChIP  = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)


# 4. is the hichip gene a TWAS gene overlapping anchor 1? # *any* anchor
t1 = find_match(match1="HiChIP.Gene_1", match2="Gene.TWAS.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="Gene.TWAS.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t3 = find_match(match1="HiChIP.Gene_1", match2="Gene.TWAS.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t4 = find_match(match1="HiChIP.Gene_2", match2="Gene.TWAS.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
df_anno$TWAS_and_HiChIP = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)


# 6. Is it an eGene even if it is not overlapping an eQTL?
# 6.a
t1 = find_match(match1="HiChIP.Gene_1", match2="eGene.TCGA.PRAD.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="eGene.TCGA.PRAD.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t3 = find_match(match1="HiChIP.Gene_1", match2="eGene.TCGA.PRAD.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t4 = find_match(match1="HiChIP.Gene_2", match2="eGene.TCGA.PRAD.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
df_anno$egeneTCGA_and_HiChIP = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)
# 6.b
t1 = find_match(match1="HiChIP.Gene_1", match2="eGene.Thibo.PRAD.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t2 = find_match(match1="HiChIP.Gene_2", match2="eGene.Thibo.PRAD.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t3 = find_match(match1="HiChIP.Gene_1", match2="eGene.Thibo.PRAD.anchor_2", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
t4 = find_match(match1="HiChIP.Gene_2", match2="eGene.Thibo.PRAD.anchor_1", df_anno, idcol="loopid", split=FALSE, summarizeby = FALSE)
id_interest = unique(c(t1$loopid[t1[,ncol(t1)]], t2$loopid[t2[,ncol(t2)]], t3$loopid[t3[,ncol(t3)]],t4$loopid[t4[,ncol(t4)]]))
df_anno$egeneThibo_and_HiChIP = ifelse(df_anno$loopid %in% id_interest, TRUE, FALSE)
# 7. reported
# "eGene.TCGA.PRAD.anchor_1"
# "eGene.Thibo.PRAD.anchor_1"
# "eGene.GTEx.Prost.anchor_1"
# "eGene.TCGA.PRAD.reported.anchor_1"


# rename PET counts
names(df_anno)[which(names(df_anno) %in% c("FitHiChIP.interactions_FitHiC_Q0.01.bed", "FitHiChIP.interactions_FitHiC_Q0.01_MergeNearContacts.bed"))]="PET_Q0.01"
# merge with other data
s = df_anno[,!(names(df_anno) %in% names(d_df))]
s = cbind.data.frame(df_anno$loopid, s)
names(s)[1] = "loopid"

if (nrow(d_df)!=nrow(df_anno)) stop("Data mismatch")
d_df = merge(d_df, s, by="loopid")

if (acetylated) {
    d_df = subset(d_df, d_df$h3k27peak_1 != "none" & d_df$h3k27peak_2 != "none")
}

message("Writing annotated loop file in: ", outFile_loops)
write.table(d_df, file=outFile_loops, quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)


###########################################
# SAVE FILES BY PROMOTER AND BY ENHANCER
# from script summarize_by_connectivity.R
###########################################

################################
#### GENE CONNECTIVITY #####
################################

message("Writing annotated loop file in: ", outFile_genes)

p_df = d_df[d_df$anchor.type_1=="P" | d_df$anchor.type_2=="P",]
gene = summarize_by2cols(p_df, col_nms=c("HiChIP.Gene_1", "HiChIP.Gene_2"), add_cols = c("PET_Q0.01", "loop.type", "Enhancerpeakid_1", "Enhancerpeakid_2"), groupby = "loopid", means_cols = "PET_Q0.01")
#To get sums or median or means, you can separate that column by doing:
# x$PET_SUM = sapply(strsplit(as.character(x$PET_counts), "[,]"),
# function(x) sum(as.numeric(x), na.rm=TRUE))
gene$nPET = sapply(strsplit(as.character(gene$PET_Q0.01), "[,]"),
function(x) sum(as.numeric(x), na.rm=TRUE))

# Add whether gene connects or overlaps a Fine-mapped SNP
x=d_df[d_df$Paintor_causal_rsid_1!="none" | d_df$Paintor_causal_rsid_2!="none",]
x1 = separate_rows(x, "HiChIP.Gene_1", sep = ",")
x2 = separate_rows(x, "HiChIP.Gene_2", sep = ",")
g = unique(c(x1$HiChIP.Gene_1, x2$HiChIP.Gene_2))
gene$GWAS = ifelse(gene$HiChIP.Gene %in% g, TRUE, FALSE)

write.table(gene, file=outFile_genes, quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)


################################
#### ENHANCER CONNECTIVITY #####
################################

message("Writing annotated loop file in: ", outFile_enhancers)
e_df = d_df[d_df$anchor.type_1=="E" | d_df$anchor.type_2=="E",]
enhancer = summarize_by2cols(e_df, col_nms=c("Enhancerpeakid_1", "Enhancerpeakid_2"), add_cols = c("PET_Q0.01", "loop.type", "HiChIP.Gene_1", "HiChIP.Gene_2"), groupby = "loopid", means_cols = "PET_Q0.01")

enhancer$nPET = sapply(strsplit(as.character(enhancer$PET_Q0.01), "[,]"),
function(x) sum(as.numeric(x), na.rm=TRUE))

# Add whether enhancer connects or overlaps a Fine-mapped SNP
d_df=d_df[d_df$Paintor_causal_rsid_1!="none" | d_df$Paintor_causal_rsid_2!="none",]
x1 = separate_rows(x, "Enhancerpeakid_1", sep = ",")
x2 = separate_rows(x, "Enhancerpeakid_2", sep = ",")
e = unique(c(x1$Enhancerpeakid_1, x2$Enhancerpeakid_2))
enhancer$GWAS = ifelse(enhancer$Enhancerpeakid %in% e, TRUE, FALSE)

write.table(enhancer, file=outFile_enhancers, quote = FALSE, sep = "\t", row.names = FALSE, col.names=TRUE)


}

# create E-P loops with GWAS at enhancer:
# ARCHIVED/annotate_step3_GWAS.R


