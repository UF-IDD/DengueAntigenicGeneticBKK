
library(argparse)
library(dplyr)
library(seqinr)
    

parser = ArgumentParser(description = "Concatenate sequences of the same taxa")
parser$add_argument('outfile')
parser$add_argument('infiles', nargs='+')

# test
inputArg = parser$parse_args(c(
    "02-processedData/AA_E_62NS2A/concat.fasta"
    , "02-processedData/AA_Envelope/1.fasta"
    , "02-processedData/AA_NS2A_beyondE/subset.fasta"
))

inputArg = parser$parse_args()


Seqs = lapply(inputArg$infiles, read.fasta)
seqnames = Reduce(lapply(Seqs, names), f = union)
if(length(Reduce(lapply(Seqs, names), f = intersect)) != length(seqnames)){
    stop("Non-total overlap of taxa between input files!")
}

Seqs = lapply(seqnames, function(x){
    lapply(Seqs, function(Seq){ Seq[[x]]}) %>%
    do.call(what = c)
})


dir.create(dirname(inputArg$outfile), recursive = T)
write.fasta(
    Seqs
    , seqnames
    , inputArg$outfile
)