#!/usr/bin/env Rscript --vanilla


    #   Generate AA fasta file by
    #   subsampling sites within each protein to meet the specified length
    #   : this sensitivity analysis would only make sense for proteins that have
    #   ; sequence longer than the specified length as redundant sites would emerge
    #   : and would be removed on colinearity reduction step
    
inputArgs = commandArgs(trailing = T)

outdir = inputArgs[1]
outLength = as.integer(inputArgs[2]) # 495 for E
nSamples = as.integer(inputArgs[3])  # number of samples to generate
protein_name = inputArgs[4] # name of protein to subsample residues from
aaver = "AA"
source("Scripts/sequence_load/aaver.R") # AA sequence



lapply(1:nSamples, function(iSample){

    aa_proteins = aa[[protein_name]]

    if(length(aa_proteins[[1]]) <= outLength){
        # return nothing if the protein is not longer than the specified length
        return(NULL)
    }
    
    set.seed(outLength * iSample)
    # scramble the site order and subset just up till the length specified
    i = sample.int(length(aa_proteins[[1]]), replace = F)[1:outLength]
    
    # generate the sampled sequence
    aa.sampled = lapply(aa_proteins, function(aa_protein){
        aa_protein[i]
    })
    
    # write to FASTA
    dir.create(file.path(outdir,iSample), recursive = T)
    write.fasta(aa.sampled
        , names = names(aa.sampled)
        , file.out = paste0(
            file.path(outdir, iSample, protein_name)
            , ".fasta"
        )
    )
})

