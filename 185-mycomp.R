library("DESeq2")
library("tximport")

##### List the files and set up metadata #####
# Note, you should change this to use the files in your home directory
files <- c("C:/Users/megan/Downloads/GSE74202_RAW/GSM1914489_C1_48hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914490_C2_48hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914491_S1_48hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914492_S2_48hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914493_C1_72hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914494_C2_72hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914495_S1_72hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914496_S2_72hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914497_C1_24hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914498_C2_24hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914499_S1_24hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914500_S2_24hpf.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914501_Sample_para_mut1.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914502_Sample_para_mut2.sorted.bam.count.ft.txt",           
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914503_Sample_wt_para1.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914504_Sample_wt_para2.sorted.bam.count.ft.txt",
           #"C:/Users/megan/Downloads/GSE74202_RAW/GSM1914505_Sample_wt_para3.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914506_Sample_ache_mut1.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914507_Sample_ache_mut2.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914508_Sample_wt_ache1.sorted.bam.count.ft.txt",
           "C:/Users/megan/Downloads/GSE74202_RAW/GSM1914509_Sample_wt_ache2.sorted.bam.count.ft.txt")

pca <- read.delim(files[20],header=FALSE)
pca <- temp[1:33737,]
for (i in c(19,18,17,16,15,14,13,12)) {
    temp <- read.delim(files[i],header=FALSE)
    temp <- temp[1:33737,]
    pca <- cbind(pca[1],temp[2],pca[2])
    print(i)
} #transcripts as rows, samples as columns

pcaS <- read.delim(files[11],header=FALSE)
for (i in c(10,9,8,7,6,5,4,3,2,1)) {
    temp <- read.delim(files[i],header=FALSE)
    pcaS <- cbind(pcaS[1],temp[2],pcaS[2])
    print(i)
} #transcripts as rows, samples as columns

for (j in 1:33737) {
    if (as.character(pca[j,1]) != as.character(pcaS[j,1])){
        pcaS <- rbind(pcaS[1:j-1,],c(pca[j,1],rep(NA,11)),pcaS[1:j,])
        print(j)
    }
}


cbind(t, z[, "symbol"][match(rownames(t), rownames(z))])

pca <- rlog(pca)

pca <- prcomp(t(MyReadCountMatrix))
summary(pca)

#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
pca.proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100

conditions <- c("unc45bWT_48hpf", "unc45bWT_48hpf",
                "unc45b-/-_48hpf", "unc45b-/-_48hpf",
                "unc45bWT_72hpf", "unc45bWT_72hpf",
                "unc45b-/-_72hpf", "unc45b-/-_72hpf",
                "unc45bWT_24hpf", "unc45bWT_24hpf",
                "unc45b-/-_24hpf", "unc45b-/-_24hpf",
                "hsp90a-/-", "hsp90a-/-",
                "hsp90aWT", "hsp90aWT",
                "ache-/-", "ache-/-",
                "acheWT", "acheWT")
samples <- data.frame("run"=c("unc45b-48-wt1", "unc45b-48-wt2",
                "unc45b-48-mut1", "unc45b-48-mut2",
                "unc45b-72-wt1", "unc45b-72-wt2",
                "unc45b-72-mut1", "unc45b-72-mut2",
                "unc45b-24-wt1", "unc45b-24-wt2",
                "unc45b-24-mut1", "unc45b-24-mut2",
                "hsp90a-72-mut1", "hsp90a-72-mut2",
                "hsp90a-72-wt1", "hsp90a-72-wt2",
                "ache-72-mut1", "ache-72-mut2",
                "ache-72-wt1", "ache-72-wt2"),
                    "condition"=conditions)
#names(files) = samples$run

##### Use "tximport" to convert RSEM results to the format needed by DESeq2 #####
#txi <- tximport(files, type="rsem")
#txi$length[txi$length == 0] <- 1 # add pseudocount of 1 to lengths to fix error with 0-length transcripts
#ddsTxi <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)

# Filter things with very low counts to we don't waste time on those
#keep <- rowSums(counts(ddsTxi)) >= 10
#ddsTxi <- ddsTxi[keep,]

##### Perform deseq2 #####
#dds <- DESeq(ddsTxi)
#res <- results(dds)

##### Write results to chow_vs_hfd_deseq2.csv #####
#write.csv(as.data.frame(res), file="~/week4/chow_vs_hfd_deseq2.csv")

#data <- read.csv("~/week4/chow_vs_hfd_deseq2.csv")
#plot(data$log2FoldChange, -log10(data$pvalue), xlim=c(-10,10), ylim=c(0,40),
#     xlab="log2 fold change", ylab="-log10 p-value", main="chow vs. hfd differential expression results",
#     type="n")
# then add the points
#sel <- which(data$padj<0.05)
#points(data[sel,"log2FoldChange"], -log10(data[sel,"pvalue"]),col="green",pch=20,cex = .6)
#sel <- which(data$padj>=0.05)
#points(data[sel,"log2FoldChange"], -log10(data[sel,"pvalue"]),col="red",pch=20,cex=.6)

plot(pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(pca.proportionvariances[2], 2), "%"))
points(pca$x, col="black", pch=16, cex=1)

