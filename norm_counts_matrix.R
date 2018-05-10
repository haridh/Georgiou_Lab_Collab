file.names <- dir("./", pattern = ".txt")
for(i in 1:length(file.names)){
p <- read.table(file.names[i],header=TRUE,sep = "\t", stringsAsFactor = FALSE)
q <- p[, c(1,7)]
if(i==1){
out <- q} else{
out <- merge(out, q, by = "Geneid")
}}
write.table(out, "Raw_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cts <- data.frame(out, row.names = 1)
cts <- round(cts)
suppressMessages(library(DESeq2))
sampleTable <- data.frame(sampleName = file.names, condition = file.names)
suppressMessages(dds <- DESeqDataSetFromMatrix(countData = cts, colData = sampleTable, design = ~ condition))
suppressMessages(dds <- DESeq(dds))
norm <- counts(dds, normalized = TRUE)
write.table(norm, "Normalized_counts.txt", sep = "\t", quote = FALSE)

