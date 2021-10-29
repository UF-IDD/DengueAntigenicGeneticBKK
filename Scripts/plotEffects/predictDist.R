#!/usr/bin/env Rscript --vanilla
# rm(list=ls()); Sys.setenv(SLURM_ARRAY_TASK_ID=1, SLURM_CPUS_ON_NODE=1)

library(argparse)
library(parallel)
library(tidyverse) #; theme_set(theme_minimal())
library(seqinr)


parser = ArgumentParser(description = "Make predictions of all distance entries given fitted effect sizes.")
parser$add_argument('distver')
parser$add_argument('fitstem')
parser$add_argument('--fasta', nargs = "+", action = "extend")
parser$add_argument('--ifold', default = Sys.getenv("SLURM_ARRAY_TASK_ID"))
parser$add_argument('--outfile', default = '')
parser$add_argument('-isSummary', action = 'store_true')


# test
inputArg = parser$parse_args(c(
    "thai_map"
    , '03-fits/nnls/thai_map/AA_NS2A_beyondE/adj_envelope/aasub/subset'
    , '--fasta', '02-processedData/AA_Envelope/1.fasta'
    , '--fasta', '02-processedData/AA_NS2A_beyondE/subset.fasta'
))
inputArg = parser$parse_args(c(
    "thai_map"
    , '03-fits/nnls/thai_map/AA/adj_none/aasub/envelope_protein_E'
    , '--fasta', '02-processedData/AA_Envelope/1.fasta'
))

    

# parse arguments
inputArg = parser$parse_args()

# set output file for predictions
if(inputArg$outfile == ''){
    outfile = with(inputArg, file.path(fitstem,ifold,'prediction.txt'))
} else {
    outfile = inputArg$outfile
}


#   Read in sequences
#   .................

Seqs = lapply(inputArg$fasta, read.fasta)

# reduce to viruses that have sequences in all listed FASTA
# (just a safety check)
# and concatenate them
Seqnames = lapply(Seqs,names) %>% 
    Reduce(f = intersect) 
Seqs = Seqnames %>%
    lapply(function(seqname){
        lapply(Seqs, function(x) x[[seqname]]) %>%
        do.call(what = c)
    }); names(Seqs) = Seqnames



#   Read in the fitted effect sizes
#   ...............................

getFilename = function(x){
    ifelse(inputArg$isSummary
        , file.path(inputArg$fitstem, paste0(x,'_summary.csv'))
        , file.path(inputArg$fitstem, inputArg$ifold, paste0(x,'.csv'))
    )
}

d = with(inputArg, getFilename('d')) %>%
    read_csv( col_types = cols(
            Feature = col_character(),
            effect = col_double(),
            clustNum = col_double(),
            pos = col_double()
        )
    )%>%
    mutate(
        clustNum = str_extract(Feature, '[0-9]+')
        , serumAA = substring(Feature, 1, 1)
        , virusAA = str_extract(Feature, '[a-z]$')
    )
if(inputArg$isSummary){
    d$effect = d$Median
}

# determine average sum of measurement noise between two identical antigens
v = with(inputArg, getFilename('v')) %>%
    read_csv %>%
    with(mean(avidity))    
p = with(inputArg, getFilename('p')) %>%
    read_csv %>%
    with(mean(potency))
m = v + p


# vector of cluster numbers (groups of collinear sites)
d.clustNum = d$clustNum
names(d.clustNum) = with(d, paste0(serumAA,pos,virusAA))

# vector of effect sizes
d = d$effect
names(d) = names(d.clustNum)



#   Read in virus pairs to predict
#   ..............................

DLong = read_csv(paste0("02-processedData/antigenicDist/", inputArg$distver, ".csv"))

apply(DLong, 1, function(x){
    mutated = which(Seqs[[ x[['serum']] ]] != Seqs[[ x[['virus']] ]])
    subs = paste0(
        Seqs[[x[['serum']]]][mutated]
        , mutated
        , Seqs[[x[['virus']]]][mutated]
    )
    sum(by(d[subs], d.clustNum[subs], unique)) + m
}) %>%
as.character %>%
writeLines(outfile)


