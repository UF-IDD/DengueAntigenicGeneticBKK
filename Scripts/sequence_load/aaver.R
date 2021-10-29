library(seqinr)

    # global variable "aaver" is required to be present

aaFiles = list.files(file.path("02-processedData/",aaver)
    , pattern="\\.fasta$"
    , full.names = T
)

aa = lapply( aaFiles, read.fasta)
names(aa) = gsub("\\.fasta$","",basename(aaFiles))
rm(aaFiles)

