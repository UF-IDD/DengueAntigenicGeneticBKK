#!/usr/bin/env Rscript --vanilla
# rm(list=ls()); Sys.setenv(SLURM_ARRAY_TASK_ID=1, SLURM_CPUS_ON_NODE=1)

library(argparse)
library(parallel)
library(tidyverse); theme_set(theme_minimal())


parser = ArgumentParser(description = "fit substitution effect sizes using nnls")
parser$add_argument('distver')
parser$add_argument('aaver')
parser$add_argument('--adjfasta')
parser$add_argument('--adjname', default = "adj_none")
parser$add_argument('--featureTypes', default = "aasub")


iProtein = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

inputArg = parser$parse_args()
attach(inputArg)

outdir = file.path("03-fits", "nnls", distver, aaver, adjname)
ncpu = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
ncpu = ifelse(is.na(ncpu), 1, ncpu)


    #   load data

dat = new.env()
source("Scripts/sequence_load/aaver.R", local=dat) # AA sequence
dat$DLong = read_csv(paste0("02-processedData/antigenicDist/", distver, ".csv")) # antigenic distance



    #   Define a package environment
    #   that can be easily passed on to the cluster

pack = new.env()
source("Scripts/configs/generic_functions.R", local=pack)
pack$outdir = outdir
pack$no_save = Sys.getenv("NOT_SAVE_ESTIMATES") == "TRUE"




    #   reduce data to ones with both
    #   sequence and antigenic data available
    #   .....................................


strains = dat$DLong[ , c('virus','serum')] %>%
    unlist %>%
    unique %>%
    intersect(names(dat$aa[[iProtein]]))

pack$aa = dat$aa[[iProtein]][strains]
pack$protein_name = names(dat$aa)[iProtein]
pack$DLong = dat$DLong %>%
    filter(virus %in% strains, serum %in% strains)

if(inputArg$adjname != "adj_none"){
    message("loading FASTA of adjustment protein...")
    pack$aa_adj = seqinr::read.fasta(inputArg$adjfasta)[strains]
} else {
    pack$aa_adj = list()
}


    #   compute genetic feature changes
    #   ...............................

cl = makeCluster(ncpu)
clusterExport(cl, ls(pack), envir = pack)

getImportance = function(protein_name, featureType, nFolds = 100){

    if(!(featureType %in% featureTypes)){ stop("Unknown feature type.") }
    
    attach(pack)

    # append adjustment proteins
    strains = names(aa)
    aa = lapply(strains, function(x){
        c(aa_adj[x], aa[x]) %>%
            unlist(use.names = F)
    })
    names(aa) = strains

    # cluster positions with exactly the same information
    subClust = data.frame(
            position = seq_along(aa[[1]])
            , pattern = aa %>%
                do.call(what = rbind) %>%
                apply(2, paste, collapse = ";")
        ) %>%
        group_by(pattern) %>%
        summarize(pos = list(position)) %>%
        mutate(clustNum = seq_along(pos))

    # remove clusters with no diversity
    subClust = subClust %>%
        filter(
            (
                strsplit(subClust$pattern %>% as.character, ";") %>%
                sapply(function(a){ unique(a) %>% length })
            ) > 1
        ) %>%
        select(-pattern) %>%
        # representative position to use when extracting features
        # in pair-wise comparison
        mutate(pos1 = sapply(pos, function(a) a[1]))

    # compute pair-wise changes
    featLong = DLong %>%
        select(virus, serum) %>%
        mutate( sub = mapply( function(virus, serum){
                featVirus = aa[[virus]][subClust$pos1]
                featSerum = aa[[serum]][subClust$pos1]
                changed = featVirus != featSerum

                # store the change;
                # non-changed will be filled with zero later
                if(featureType == "aachange"){
                    out = subClust$clustNum[changed]
                } else if(featureType == "aasub"){
                    out = paste0(
                        featSerum[changed]
                        , subClust$clustNum[changed]
                        , featVirus[changed]
                    )
                } else if(featureType == "polardiff"){
                    out = aapolar[featVirus] - aapolar[featSerum]
                }
                out
            }
            ,virus = virus
            ,serum = serum
            ,SIMPLIFY = FALSE
        ))

    library(Matrix)
    if(featureType != "polardiff"){
        # make genetic feature matrix
        uniqueFeatures = featLong$sub %>%
            unlist %>%
            unique
        iFeat = seq_along(uniqueFeatures)
        names(iFeat) = uniqueFeatures
        features = sparseMatrix(
            i = rep(seq_along(featLong$sub), sapply(featLong$sub, length))
            , j = iFeat[unlist(featLong$sub)]
            , dims = c(nrow(featLong), length(uniqueFeatures))
        )
        colnames(features) = uniqueFeatures

    } else {
        features = featLong$sub %>%
            do.call(what = rbind)
        colnames(features) = subClust$clustNum    
    }



        #   Run nnls
        #   ........

    clusterExport(cl, c(
        "featureType"
        , "protein_name"
        , "subClust"
        , "features"
    ), envir = environment())

    parLapply(cl, 1:nFolds, function(iFold){
        library(tidyverse)
        library(Matrix)
    
        # identify test sets
        set.seed((iFold + 100)*6^4)
        iTest = sample.int( nrow(features)
            , size = round(nrow(features)/10)
            , replace = F
        )

        distTrain = DLong[-iTest, ]
        # substitutions found in the training set
        subTrain = colnames(features)[ (features[-iTest, ] %>% colSums) >= 2 ]

        dNames = lapply(featLong %>% select(virus, serum), unique)
        A.make = function(feat, d){
            feat[  , subTrain ] %>%
                cbind(
                    lapply(dNames$virus, function(x) d$virus == x ) %>%
                    do.call(what = cbind)
                ) %>%
                cbind(
                    lapply(dNames$serum, function(x) d$serum == x ) %>%
                    do.call(what = cbind)
                )
        }




        #   Fit substitution effects with nnls

        #   create objects to conform with my re-derived casting
        #   of the quadratic problem in Neher,2016 in which
        #   Q = (t(A) %*% A) + M
        #   q = -HA + (1/2)*[ lambda, 0, 0 ]
        #   where M is the diagonal matrix of 0, delta, gamma
        #   ....................................................

        H = matrix(distTrain$D, nrow=1) # vector of measured titers

        A = A.make(features[ -iTest, ], distTrain)

        M = Matrix(0, ncol(A), ncol(A), sparse = T)
        diag(M) = c(
             rep(0  , length(subTrain))
            ,rep(0.6, length(dNames$virus))    # kappa = 0.6
            ,rep(1.2, length(dNames$serum))    # delta = 1.2
        )
        Q = (t(A) %*% A) + M

        q = (- H %*% A) + c(
             rep(3, length(subTrain))      # lambda = 3
            ,rep(0, length(dNames$virus) + length(dNames$serum))
        )

        library(lsei)
        # Lawson-Hanson implementation of an algorithm for 
        # non-negative least-squares
        # : allowing the combination of non-negative and non-positive constraints
        # : solves min || Ax − b ||2 with the constraint x ≥ 0
        fit = pnnqp(
            q = Q %>% as.matrix
            ,p = t(q) %>% as.matrix
        )

        est = split(fit$x, c(
            rep('d', length(subTrain))
            ,rep('v', length(dNames$virus))
            ,rep('p', length(dNames$serum))
        ))
        names(est$d) = subTrain
        names(est$v) = dNames$virus
        names(est$p) = dNames$serum
        
    
        predictDist = A %*% do.call(c,est)
        rmseTrain = mean((predictDist - DLong$D[-iTest])^2, na.rm=T) %>% sqrt

        A.test = A.make(features[iTest, ], DLong[iTest, ])
        predictDist = A.test %*% do.call(c,est)
        rmseTest = mean((predictDist - DLong$D[iTest])^2) %>% sqrt
    

        # write estimates
        if(!no_save){
            outstem = paste0(gsub(" ", "_", file.path( outdir, featureType, protein_name, iFold )))
            data.frame(Feature = names(est$d), effect = est$d) %>%
                mutate(clustNum = gsub("[^0-9]", "", Feature) %>% as.integer) %>%
                left_join(
                    subClust %>%
                    select(-pos1) %>%
                    unnest(cols = pos)
                    , by = "clustNum"
                ) %>%
                arrange(desc(effect)) %>%
                write_csv_mkdir(file.path(outstem,"d"))
            data.frame(virus = names(est$v), avidity = est$v) %>%
                arrange(desc(avidity)) %>%
                write_csv(file.path(outstem, "v.csv"))
            data.frame(serum = names(est$p), potency = est$p) %>%
                arrange(desc(potency)) %>%
                write_csv(file.path(outstem, "p.csv"))
        }

        # return
        return(data.frame(
            featureType = featureType
            , protein_name = protein_name
            , fold = iFold
            , train_rmse = rmseTrain %>% round(4)
            , test_rmse = rmseTest %>% round(4)
        ))
    }) %>%
    do.call(what = rbind)
}


# run model fitting
lapply(featureTypes, function(featureType){
    protein_name = names(dat$aa)[iProtein]
    outfile = paste0(gsub(" ", "_", file.path( outdir, featureType, protein_name, "rmse" )))
    if(file.exists(outfile)){ return(NULL) }
    rmse = getImportance(protein_name, featureType = featureType ) %>%
        as.data.frame %>%
        pack$write_csv_mkdir(outfile)
})

stopCluster(cl)

