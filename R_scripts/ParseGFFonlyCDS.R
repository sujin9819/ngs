suppressPackageStartupMessages({library(optparse)})
suppressPackageStartupMessages({library(stringr)})

option_list <- list(
    make_option(c("-i", "--input"), type="character", default="./", help="input gff file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="./", help="output file name", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

if (is.null(options$input)) {
	print ("Error: No select file.")
	stop()
} else {
	cat ("select input file name: ", options$input, "\n")
	in_file <- options$input
}
if (is.null(options$output)) {
	print ("Error: No output file name.")
	stop()
} else {
	cat ("out file name: ", options$output, "\n")
	output <- options$output
}

dir <- getwd()
gff <- read.delim(paste(dir,in_file,sep='/'), header=F, sep="\t",comment.char="#")
gff_splited <- str_split_fixed(gff$V9, ";",2)[,1] 
gff_splited <- str_split_fixed(gff_splited, "=",2)[,2]
gff[,9] <- gff_splited
write.table(gff,output,row.names=F,col.names=F,sep="\t",quote=F)
