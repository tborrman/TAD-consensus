#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Plot TAD boundaries for all 60 Hi-C samples for given chromosome and deep seq K562")
parser$add_argument("-c", help="chromosome (ex. chr1)", type="character", dest="c", required=TRUE)
args <- parser$parse_args()

# Get 60 Hi-C samples
hic_names <- read.table("hic_samples.txt", header=TRUE, sep="\t")
prefix <- "__hg19__genome__C-40000-iced__"
suffix <- "--is520000--nt0--ids320000--ss0--immean--mbs0--mts0.tads.bed"


get_conservation <- function(i, df) {
	# Return conservation for start 
	# and end boundaries in row i of df

	# Get block of replicates br
	cellname <- df[i,"cell"]
	br <- df[df$cell == cellname,]
	# Get conservation of start boundary
	boundary_start <- df[i, "TADstart"]
	conserved <- br[br$TADstart == boundary_start | br$TADend == boundary_start,]
	# if (nrow(conserved) > 3) {
	# 	print(conserved)
	# 	print('ERROR: more than 3 replicate types')
	# }
	cons_rep_start <- levels(as.factor(as.character(conserved$biorep)))
	cons_rep_start <- paste(cons_rep_start[order(cons_rep_start)], collapse="_")

	# Get conservation of end boundary
	boundary_end <- df[i, "TADend"]
	conserved <- br[br$TADstart == boundary_end | br$TADend == boundary_end,]
	cons_rep_end <- levels(as.factor(as.character(conserved$biorep)))
	cons_rep_end <- paste(cons_rep_end[order(cons_rep_end)], collapse="_")

	return(c(cons_rep_start, cons_rep_end))
}

add_colors <- function(df) {
	# Add conservation colors to dataframe
	# P: steelblue1
	# R1 : cyan
	# R2 : dodgerblue
	# R1_R2: mediumorchid2
	# P_R1: purple
	# P_R2: magenta2
	# P_R1_R2: red
	
	color_map <- c("steelblue1", "cyan", "dodgerblue",
		"mediumorchid2", "purple", "magenta2", "red")  
	names(color_map) <- c("P", "R1", "R2", "R1_R2", "P_R1",
		"P_R2", "P_R1_R2")

	constart <- df$cons_start
	consend <- df$cons_end
	# Order errors when trying to use simple indexing or sapply
	# Using stupid loops instead, I hate R
	cons_start_color <- c()
	cons_end_color <- c()
	for (cst in constart) {
		cons_start_color <- c(cons_start_color, color_map[cst])
	}
	for (csd in consend) {
		cons_end_color <- c(cons_end_color, color_map[csd])
	}
	cdf <- cbind(df, cons_start_color, cons_end_color)
	return(cdf)
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
	if (i%%250 == 0) {
		print(i)
	}
	clist <- get_conservation(i, hic_df)
	cons_start <- c(cons_start, clist[1])
	cons_end <- c(cons_end, clist[2])
}

hic_df <- cbind(hic_df, cons_start, cons_end)
hic_df <- add_colors(hic_df)

# K562
K562_df <- read.table(paste("Houda_Ctrl_DpnII_K562.40000_", args$c,".txt", suffix, sep=""), 
                     header=FALSE, skip=1, sep="\t")
colnames(K562_df) <- c("chrom", "TADstart", "TADend", "description", "score")
K562_df$TADstart <- K562_df$TADstart - 1

# write.table(hic_df, paste("plot_boundaries_table_", args$c, ".txt", sep=""), 
# 	sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE )

# hic_df <- read.table(paste("plot_boundaries_table_", args$c, ".txt", sep=""), 
# 	sep="\t", header=TRUE, comment.char="")

# Plot
options(scipen=100)
Mb <- 1000000
png(paste("plot_boundaries_k562/plot_boundaries_", args$c, ".png", sep=""), height=1700, width=5000, res= 300)
par(mar=c(5, 6, 4, 2) + 0.1)
cells <- levels(hic_df$cell)
cell_df <- hic_df[hic_df$cell == cells[1],]
plot(cell_df$TADstart / Mb, rep(1, nrow(cell_df)), col= as.character(cell_df$cons_start_color), pch= "|", 
	xlab= "bp (Mb)" , ylim= c(0,21), main= args$c, yaxt= "n", ylab="")
axis(2, at=1:21, labels=c(cells, "K562"), las=2) 
points(cell_df$TADend / Mb, rep(1, nrow(cell_df)), col= as.character(cell_df$cons_end_color), pch= "|")

for (i in 2:length(cells)) {
	cell_df <- hic_df[hic_df$cell == cells[i],]
	points(cell_df$TADstart / Mb, rep(i, nrow(cell_df)), col= as.character(cell_df$cons_start_color), pch= "|")
	points(cell_df$TADend / Mb, rep(i, nrow(cell_df)), col= as.character(cell_df$cons_end_color), pch= "|")
}
points(K562_df$TADstart / Mb, rep(length(cells) + 1, nrow(K562_df)),col="black", pch="|")
points(K562_df$TADend / Mb, rep(length(cells) + 1, nrow(K562_df)),col="black", pch="|")
dev.off()
