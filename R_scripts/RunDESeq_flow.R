suppressPackageStartupMessages({library(optparse)})
suppressPackageStartupMessages({library(stringr)})

option_list = list(
make_option(c("-i", "--infile"), type="character", default="./", help="Input file", metavar="character"),
make_option(c("-d", "--design_table"), type="character", default="./", help="Design sheet", metavar="character"),
make_option(c("-o", "--outpath"), type="character", default="./", help="Path to output files [default= %default]", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

if (is.null(options$infile)) {
	print ("Error: No input file specified.")
	stop()
} else {
	cat ("Input coverage file: ", options$infile, "\n")
	in_file <- options$infile
}


if (is.null(options$design_sheet)) {
	print ("Error: No design_table specified.")
	stop()
} else {
	cat ("Design_table: ", options$design_table, "\n")
	d_table <- options$design_table
}


if (is.null(options$outpath)) {
	print ("Warning: No path to output files specified. Set as default: current directory")
} else {
	cat ("Path to output files: ", options$outpath, "\n")
	out_path <- options$outpath
}



suppressPackageStartupMessages({library(DESeq2)})

count_data <- read.table(in_file, sep="\t", row.names = 1, header=T)
design_table <- read.table(d_table, sep="\t", row.names = 1, header=T)
conditions <- as.factor(design_table[,1])
exp_conditions <- data.frame(row.names=colnames(count_data), condition = conditions)

Deseq_dataset <- DESeqDataSetFromMatrix(countData = count_data, colData = exp_conditions, design = ~condition)
Deseq_dds <- DESeq(Deseq_dataset)
Deseq_res <- results(Deseq_dds)
norm.count <- counts(Deseq_dds, normalized=TRUE)
print ("Write normalized count table.")
write.table(as.data.frame(norm.count), file=paste(out_path,"Deseq_normcount.txt",sep=""),sep="\t")
print ("Done.")

write.csv(as.data.frame(Deseq_res), file=paste(out_path,"Deseq2_results.csv",sep=""))

pdf(paste(out_path, "Dispersion.pdf", sep=""))
plotDispEsts(Deseq_dds)
dev.off()

rlog_tf <- rlog(Deseq_dds)
distances <- dist(t(assay(rlog_tf)))
distance_matrix <- as.matrix(distances)
rownames(distance_matrix) <- colnames(distance_matrix)

suppressPackageStartupMessages({library(RColorBrewer)})
suppressPackageStartupMessages({library(gplots)})

palette <-colorRampPalette(brewer.pal(9,"GnBu"))(100)
hierarchy <- hclust(distances)
suppressPackageStartupMessages({library(pheatmap)})
print ("Print Dendrogram_and_heatmap.")
pdf(paste(out_path, "Dendrogram_and_heatmap.pdf", sep=""))
heatmap.2(distance_matrix,Rowv=as.dendrogram(hierarchy), symm=TRUE, trace="none", col=rev(palette))
dev.off()
print ("Done.")
Norm_counts <- norm.count[order(rowSums(norm.count),decreasing = T),] #rowsums 내림차순으로 정렬
Norm_counts <- Norm_counts[rowSums(Norm_counts)>0,]
#Norm_counts <- sapply(count_data, function(x){(x-min(x))/(max(x)-min(x))})

print ("Print each gene heatmap")
pdf(paste(out_path, "Geneheatmap.pdf", sep=""))
pheatmap(Norm_counts[1:100,],palette,show_rownames = F,border_color = NA) 
dev.off()
print ("Done.")

unique_cond <- unique(conditions)
comparison_condi <- paste(unique_cond[1],"vs",unique_cond[2],sep="")
suppressPackageStartupMessages({library(EnhancedVolcano)})
print ("Print volcano plot.")
pdf(paste(out_path, "volcanoplot.pdf", sep=""))
EnhancedVolcano(Deseq_res, lab= rownames(Deseq_res),x="log2FoldChange",y="pvalue", pCutoff =0.05, FCcutoff = 1,legendLabSize=10,title =paste(comparison_condi," Volcano plot",sep=""))
dev.off()
print ("Done.")


print ("Print PCA plot.")
pdf(paste(out_path, "PCA.pdf", sep=""))
plotPCA(rlog_tf, intgroup=c("condition"))
dev.off()
print ("Done.")

