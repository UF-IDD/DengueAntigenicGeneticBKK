#!/usr/bin/env Rscript --vanilla


    #   Generate AA fasta file by
    #   subsetting specified sites within specified protein(s).
    #   Defaults to sampling from the whole polyprotein.

library(argparse)
library(dplyr)
library(readr)


parser = ArgumentParser(description = "subset input protein(s) to specified sites")
parser$add_argument('outdir')
parser$add_argument('--positions', nargs='+', type = 'integer')
parser$add_argument('--positionFile')
parser$add_argument('--exclude', action = "store_true", default = F
    , help = "Exclude specified sites instead of including them."
)
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
parser$add_argument('--aaver', default = "AA")
parser$add_argument('--subname', default = "subset")
parser$add_argument('--untrack', action = "store_true", default = F
    , help = "Do not track source of AA and positions subsetted."
)

inputArg = parser$parse_args()


# combine sites specified from file and command line input
inputArg$positions = c(
        inputArg$positions
        , readLines(inputArg$positionFile)
    ) %>%
    unique %>%
    as.integer %>%
    sort

attach(inputArg)
dir.create(outdir, recursive = T)


# import AA sequences
source("Scripts/sequence_load/aaver.R")

# concatenate protein sequences of the same entry
aa_proteins = lapply(seq_along(aa[[1]]), function(i){
    lapply(aa[ip], function(aap) aap[[i]]) %>%
    unlist(use.names = F)
})

# if specified sites should be exluded instead of included
if(exclude){
    positions = setdiff(seq_along(aa_proteins[[1]]), positions)
}

# document source of positions
if(!untrack){
    pLength = sapply(aa[ip], function(aap){
        aap[[1]] %>% length
    })
    data.frame(
        protein_name = names(pLength)
        , posStart = cumsum(pLength) - pLength + 1
        , posEnd = cumsum(pLength)
    ) %>%
    write_csv(file.path(outdir, "position_source.csv"))
}


# perform subsetting
aa.sampled = lapply(aa_proteins, function(aa_protein){
    aa_protein[positions]
})

# write to FASTA
write.fasta(aa.sampled
    , names = names(aa[[1]])
    , file.out = file.path(outdir, paste0(subname, ".fasta"))
)

# log sampled positions
if(!untrack){
    writeLines( as.character(positions)
        , file.path(outdir, paste0(subname, ".txt"))
    )
}