#!/usr/bin/env Rscript --vanilla


    #   Generate AA fasta file by
    #   randomly sampling from residues of all proteins to meet the specified length

library(tidyverse)
    
inputArgs = commandArgs(trailing = T)

outdir = inputArgs[1]
outLength = as.integer(inputArgs[2]) # 901 for NS5
nSamples = as.integer(inputArgs[3])  # number of samples to generate
aaver = "AA"
source("Scripts/sequence_load/aaver.R") # AA sequence


# concatenate protein sequences of the same entry
# to recover the polyprotein
iPoly = c(
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
aa.concat = lapply(seq_along(aa[[1]]), function(i){
    lapply(aa[iPoly], function(aap) aap[[i]]) %>%
    unlist(use.names = F)
})


dir.create(outdir, recursive = T)

lapply(1:nSamples, function(iSample){

    set.seed(outLength * iSample)

    # generate positions of residues to include
    i = sample.int(n = length(aa.concat[[1]]), size = outLength, replace = F)

    lapply(aa.concat, function(aap) aap[i]) %>%
        write.fasta(
            names = names(aa[[1]])
            , file.out = paste0(
                file.path(outdir, iSample)
                , ".fasta"
            )
        )
})
