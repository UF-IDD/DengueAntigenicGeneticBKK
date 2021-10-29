#!/usr/bin/env Rscript --vanilla


    #   Generate AA fasta file by
    #   subsampling sites within specified protein(s) to meet the specified length.
    #   Defaults to sampling from the whole polyprotein.
    #   : this sensitivity analysis would only make sense for proteins that have
    #   ; sequence longer than the specified length as redundant sites would emerge
    #   : and would be removed on colinearity reduction step

library(argparse)
library(dplyr)
library(readr)

aaver = "AA"
source("Scripts/sequence_load/aaver.R") # AA sequence

parser = ArgumentParser(description = "sample input protein(s) to specified length")
parser$add_argument('outdir')
parser$add_argument('--outLength', required = T, type = 'integer')
parser$add_argument('--nSamples', required = T, type = 'integer')
parser$add_argument('-ip', type = 'integer', nargs='+'
    , default = c(
        1   # anchored capsid
        , 5 # prM
        , 3 # envelope
        , 6 # NS1
        , 7 # NS2A
        , 8 # NS2B
        , 9 # NS3
        , 10 # NS4A
        , 12 # 2K
        , 11 # NS4B
        , 14 # NS5
    )
)
parser$add_argument('--sortsite', action = 'store_true', default = F
    , help = "sort the sampled sites by their source positions"
)

inputArg = parser$parse_args()
attach(inputArg)
dir.create(outdir, recursive = T)

# concatenate protein sequences of the same entry
aa_proteins = lapply(seq_along(aa[[1]]), function(i){
    lapply(aa[ip], function(aap) aap[[i]]) %>%
    unlist(use.names = F)
})

# document source of positions
pLength = sapply(aa[ip], function(aap){
    aap[[1]] %>% length
})
data.frame(
    protein_name = names(pLength)
    , posStart = cumsum(pLength) - pLength + 1
    , posEnd = cumsum(pLength)
) %>%
write_csv(file.path(outdir, "position_source.csv"))


# perform sampling
lapply(1:nSamples, function(iSample){

    if(length(aa_proteins[[1]]) < outLength){
        # return nothing if the protein is not longer than the specified length
        return(NULL)
    }
    
    set.seed(outLength * iSample)
    # generate positions of residues to include
    i = sample.int(n = length(aa_proteins[[1]]), size = outLength, replace = F)
    if(sortsite){ i = sort(i) }
    
    # generate the sampled sequence
    aa.sampled = lapply(aa_proteins, function(aa_protein){
        aa_protein[i]
    })
    
    # write to FASTA
    write.fasta(aa.sampled
        , names = names(aa[[1]])
        , file.out = paste0(
            file.path(outdir, iSample)
            , ".fasta"
        )
    )
    
    # log sampled positions
    writeLines(as.character(i), paste0(file.path(outdir, iSample),".txt"))
})

