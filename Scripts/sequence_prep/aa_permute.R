

    #   Generate AA fasta file by
    #   permuting residues across viruses within each site.


library(argparse)
library(dplyr)
library(readr)
    

parser = ArgumentParser(description = "permute residues of specified protein(s) across viruses")
parser$add_argument('outdir')
parser$add_argument('--aaver', default = "AA")
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
parser$add_argument('-n', type = 'integer'
    , required = T
    , help = "number of permutations to perform"
)
parser$add_argument('--seed', type = 'integer'
    , help = "set seed for reproducibility"
)


inputArg = parser$parse_args()

attach(inputArg)
dir.create(outdir, recursive = T)


# import AA sequences
source("Scripts/sequence_load/aaver.R")

# concatenate protein sequences of the same entry
aa_proteins = lapply(seq_along(aa[[1]]), function(i){
    lapply(aa[ip], function(aap) aap[[i]]) %>%
    unlist(use.names = F)
})


if(is.null(seed)){ seed = runif(1) * 997620 }
set.seed(seed)

lapply(1:n, function(iperm){
    # permute
    aa.perm = aa_proteins %>%
        do.call(what = rbind) %>%
        apply(2, sample, size = length(aa_proteins), replace = F) %>%
        t %>%
        as.data.frame(stringsAsFactors = F) %>%
        as.list

    # write to FASTA
    write.fasta(aa.perm
        , names = names(aa[[1]])
        , file.out = file.path(outdir, paste0(iperm, ".fasta"))
    )
})

