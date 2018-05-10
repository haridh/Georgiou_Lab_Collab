map <- read.csv("file_map.csv", stringsAsFactor = FALSE, header = TRUE)
map$File <- paste(map$File, ".sort.dup.bam", sep = "")
map <- data.frame(map, row.names = 2)
norm <- read.table("Normalized_counts.txt", sep = "\t", stringsAsFactor = FALSE, header = TRUE)
m <- norm[grepl("FCGR2+", rownames(norm)), ]
s <- colSums(m)
s <- data.frame(s)
names(s)[1] <- paste("FCGR2")
s <- s/15.545
write.table(s, "FCGR2_norm_counts.txt", sep = "\t", quote = FALSE, col.names  = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCGR2 <- smap$FCGR2 + 1
smap$log_FCGR2 <- log(smap$FCGR2, 2)
ma <- ifelse(max(smap$log_FCGR2) > 10, max(smap$log_FCGR2), 10)
library(ggplot2)
pdf("FCGR2_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCGR2, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()


m <- norm[grepl("FCGR3+", rownames(norm)), ]
s <- colSums(m)
s <- data.frame(s)
names(s)[1] <- paste("FCGR3")
s <- s/6.031
write.table(s, "FCGR3_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCGR3 <- smap$FCGR3 + 1
smap$log_FCGR3 <- log(smap$FCGR3, 2)
ma <- ifelse(max(smap$log_FCGR3) > 10, max(smap$log_FCGR3), 10)
pdf("FCGR3_exp.pdf")
library(ggplot2)
p <- ggplot(smap, aes(Sample, smap$log_FCGR3, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL1", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL1")
s <- s/3.935
write.table(s, "FCRL1_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL1 <- smap$FCRL1 + 1
smap$log_FCRL1 <- log(smap$FCRL1, 2)
ma <- ifelse(max(smap$log_FCRL1) > 10, max(smap$log_FCRL1), 10)
pdf("FCRL1_exp.pdf")
library(ggplot2)
p <- ggplot(smap, aes(Sample, smap$log_FCRL1, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL2", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL2")
s <- s/5.058
write.table(s, "FCRL2_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL2 <- smap$FCRL2 + 1
smap$log_FCRL2 <- log(smap$FCRL2, 2)
ma <- ifelse(max(smap$log_FCRL2) > 10, max(smap$log_FCRL2), 10)
pdf("FCRL2_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCRL2, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL3", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL3")
s <- s/6.538
write.table(s, "FCRL3_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL3 <- smap$FCRL3 + 1
smap$log_FCRL3 <- log(smap$FCRL3, 2)
ma <- ifelse(max(smap$log_FCRL3) > 10, max(smap$log_FCRL3), 10)
pdf("FCRL3_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCRL3, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL4", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL4")
s <- s/3.875
write.table(s, "FCRL4_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL4 <- smap$FCRL4 + 1
smap$log_FCRL4 <- log(smap$FCRL4, 2)
ma <- ifelse(max(smap$log_FCRL4) > 10, max(smap$log_FCRL4), 10)
pdf("FCRL4_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCRL4, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL5", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL5")
s <- s/12.419
write.table(s, "FCRL5_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL5 <- smap$FCRL5 + 1
smap$log_FCRL5 <- log(smap$FCRL5, 2)
ma <- ifelse(max(smap$log_FCRL5) > 10, max(smap$log_FCRL5), 10)
pdf("FCRL5_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCRL5, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

m <- norm[grepl("FCRL6", rownames(norm)), ]
s <- t(m)
s <- data.frame(s)
names(s)[1] <- paste("FCRL6")
s <- s/2.671
write.table(s, "FCRL6_norm_counts.txt", sep = "\t", quote = FALSE, col.names = TRUE)
smap <- merge(s, map, by = "row.names")
smap$FCRL6 <- smap$FCRL6 + 1
smap$log_FCRL6 <- log(smap$FCRL6, 2)
ma <- ifelse(max(smap$log_FCRL6) > 10, max(smap$log_FCRL6), 10)
pdf("FCRL6_exp.pdf")
p <- ggplot(smap, aes(Sample, smap$log_FCRL6, fill = Sample))
p + geom_boxplot(outlier.colour = "white") + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, binwidth = 0.4) + theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid.major = element_blank()) + scale_y_continuous(breaks = seq(0,10,1)) + coord_cartesian(ylim = c(0,ma))
dev.off()

