#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Plot TAD boundaries for all 60 Hi-C samples for given chromosome")
parser$add_argument("-c", help="chromosome (ex. chr1)", type="character", dest="c", required=TRUE)
args <- parser$parse_args()

# Get 60 Hi-C samples
hic_names <- read.table("hic_samples.txt", header=TRUE, sep="\t")
prefix <- "__hg19__genome__C-40000-iced__"
suffix <- "--is520000--nt0--ids320000--ss0--immean--mbs0--mts0.tads.bed"


get_conservation <- function(i, df) {
	# Return conservation for start 
	# and end boundaries in row i of df

	# R1 : cyan
	# R2 : dodgerblue
	# P: steelblue1
	# R1,R2: mediumorchid2
	# R1,P: purple
	# R2,P: magenta2
	# R1,R2,P: red

	# Get block of replicates br
	cellname <- df[i,"cell"]
	br <- df[df$cell == cellname,]
	# Get conservaion of start boundary
	boundary_start <- df[i, "TADstart"]
	conserved <- br[br$TADstart == boundary_start | br$TADend == boundary_start,]
	if (nrow(conserved) > 3) {
		print(conserved)
		print('ERROR: more than 3 replicate types')
		quit()
	}
	cons_rep_start <- as.character(conserved$biorep)
	cons_rep_start <- paste(cons_rep_start[order(cons_rep_start)], collapse="_")
	# Get conservaion of end boundary
	boundary_end <- df[i, "TADend"]
	conserved <- br[br$TADstart == boundary_end | br$TADend == boundary_end,]
	if (nrow(conserved) > 3) {
		print(conserved)
		print('ERROR: more than 3 replicate types')
		quit()
	}	
	cons_rep_end <- as.character(conserved$biorep)
	cons_rep_end <- paste(cons_rep_end[order(cons_rep_end)], collapse="_")

	return(c(cons_rep_start, cons_rep_end))
}
# Build Hi-C TAD dataframe
for (i in 1:nrow(hic_names)) {
	if (!exists("hic_df")) {
		hic_df <- read.table(paste(hic_names$sample_name[i], prefix, args$c, suffix, sep=""), 
			header=FALSE, skip=1, sep="\t")
		sample <- rep(hic_names$sample_name[i], nrow(hic_df))
		cell <- rep(hic_names$cell_type[i], nrow(hic_df))
		biorep <- rep(hic_names$replicate[i], nrow(hic_df))
		hic_df <- cbind(hic_df, sample, cell, biorep)
		colnames(hic_df) <- c("chrom", "TADstart", "TADend", "description", "score",
		 "sample", "cell", "biorep")
	}
	else {
		tmp_df <- read.table(paste(hic_names$sample_name[i], prefix, args$c, suffix, sep=""), 
			header=FALSE, skip=1, sep="\t")
		sample <- rep(hic_names$sample_name[i], nrow(tmp_df))
		cell <- rep(hic_names$cell_type[i], nrow(tmp_df))
		biorep <- rep(hic_names$replicate[i], nrow(tmp_df))
		tmp_df <- cbind(tmp_df, sample, cell, biorep)
		colnames(tmp_df) <- c("chrom", "TADstart", "TADend", "description", "score",
		 "sample", "cell", "biorep")
		hic_df <- rbind(hic_df, tmp_df)
	}	
}

# Subtract 1 bp from TADstart for matching conservation 
# between boundaries at start and end of TAD
hic_df$TADstart <- hic_df$TADstart - 1


# Find conservation
cons_start <- c();
cons_end <- c();
for (i in 1:nrow(hic_df)){
	if (i%%1000 == 0) {
		print(i)
	}
	clist <- get_conservation(i, hic_df)
	cons_start <- c(cons_start, clist[1])
	cons_end <- c(cons_end, clist[2])
}

hic_df <- cbind(hic_df, cons_start, cons_end)
head(hic_df,2000)