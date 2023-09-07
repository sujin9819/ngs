suppressPackageStartupMessages({library(optparse)})
suppressPackageStartupMessages({library(stringr)})

option_list <- list(
    make_option(c("-r", "--ref"), type="character", default=NULL, help="input keyword", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

if (is.null(options$ref)) {
	print ("Error: No select file.")
	stop()
} else {
	cat ("select input keyword: ", options$ref, "\n")
	in_file <- options$ref
}

dir <- getwd()
file_list <- list.files(dir,pattern="cov") 
group <-in_file
select <- grep(group ,file_list)
select_file <- file_list[select]
temp <- read.delim(paste(dir,select_file[1],sep='/'), header=F, sep="\t")
data <- data.frame(temp[,c(1,9)])
colnames(data)<-c("contig","locus")
n=3
for (file in select_file){
	temp <- read.delim(paste(dir,file,sep='/'), header=F, sep="\t")
	temp <- temp[,c(9,10)]
	colnames(temp) <- c("locus",file)
	data <- merge(data, temp,by='locus',all.y=TRUE)
	n <- n+1
}
colnames(data)<-gsub(".cov","",colnames(data))
df_sum <- data.frame(colSums(data[-c(1,2)]))
annot <- data[,c(1:2)]
annot$n <- str_split_fixed(annot$locus,"_",2)[,2]
annot$accession <- paste(annot$contig,"_",annot$n,sep="")
annot <- annot[,-c(2,3)]
contig_count<-data[,-1]
contig_count[, -1] = sapply(contig_count[, -1], FUN = "as.numeric")
contig_count <- aggregate(contig_count[,-1], by=list(contig = contig_count$contig), FUN=sum)
count <- merge(data, annot,by='locus',all.y=TRUE)
write.table(data,paste(group,"data.count.txt",sep=''),row.names=F,sep='\t',quote=F)
write.table(count,paste(group,"data.count",sep=''),row.names=F,sep="\t",quote=F)
write.table(df_sum,paste(group,"total.count",sep=''),row.names=T,sep='\t',quote=F)
write.table(contig_count,paste(group,"contig.count",sep=''),row.names=F,sep='\t',quote=F)
