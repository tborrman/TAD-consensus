k562_ins <- read.table("Houda_Ctrl_DpnII_K562.40000_chr10.txt--is520000--nt0--ids320000--ss0--immean.insulation.bedGraph",
                       sep="\t", header=FALSE, skip=1)
colnames(k562_ins) <- c("chr","start", "end", "insulation")

k562_b <- read.table("Houda_Ctrl_DpnII_K562.40000_chr10.txt--is520000--nt0--ids320000--ss0--immean.insulation.boundaries.bed",
                     sep="\t", header=FALSE, skip=1)
colnames(k562_b) <- c("chr", "start", "end", "description", "score")

k562_ins_win <- k562_ins[k562_ins$start >=70000000 & k562_ins$end <=80000000,]
k562_b_win <- k562_b[k562_b$start >=70000000 & k562_b$end <=80000000,]

pdf("k562_insulation_chr5_50Mb-80Mb.pdf", height = 4, width=20)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(k562_ins_win$start/1000000 , k562_ins_win$insulation, type="o", col="blue",
     pch=21, bg="blue", xlab = "bp (Mb)", ylab = "Insulation score", cex.lab=1.5,
     ylim=c(-1.5, 1.0))
abline(v=k562_b_win$start/1000000, col="red", lwd =1.5)
dev.off()

A549_ins <- read.table("ENCODE3-A549C-HindIII-R1__hg19__genome__C-40000-iced__chr10--is520000--nt0--ids320000--ss0--immean.insulation.bedGraph",
                       sep="\t", header=FALSE, skip=1)
colnames(A549_ins) <- c("chr","start", "end", "insulation")

A549_b <- read.table("ENCODE3-A549C-HindIII-R1__hg19__genome__C-40000-iced__chr10--is520000--nt0--ids320000--ss0--immean.insulation.boundaries.bed",
                     sep="\t", header=FALSE, skip=1)
colnames(A549_b) <- c("chr", "start", "end", "description", "score")

A549_ins_win <- A549_ins[A549_ins$start >=70000000 & A549_ins$end <=80000000,]
A549_b_win <- A549_b[A549_b$start >=70000000 & A549_b$end <=80000000,]

pdf("A549_insulation_chr5_50Mb-80Mb.pdf", height = 4, width=20)
par(mar=c(5, 6, 4, 2) + 0.1)
plot(A549_ins_win$start/1000000 , A549_ins_win$insulation, type="o", col="blue",
     pch=21, bg="blue", xlab = "bp (Mb)", ylab = "Insulation score", cex.lab=1.5,
     ylim=c(-1.5, 1.0))
abline(v=A549_b_win$start/1000000, col="red", lwd =1.5)
dev.off()
