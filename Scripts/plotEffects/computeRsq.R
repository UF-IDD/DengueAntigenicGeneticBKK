
library(tidyverse)

instem = "02-processedData/antigenicDist/thai_map"


# read in mapped distances
DLong = read_csv(paste0(instem, ".csv"), show_col_types = FALSE) %>%
    mutate(
        virusType = substring(virus, 1, 5)
        , serumType = substring(serum, 1, 5)
    ) %>%
    mutate(pairType = ifelse(virusType == serumType, 'intra', 'inter'))


# load in predictions from serotype centroid distances
Dist.avg = readLines(paste0(instem,'.dist')) %>%
    as.numeric
sse = DLong %>%
    mutate(centroid = Dist.avg) %>%
    group_by(pairType) %>%
    summarize(sse = sum((centroid - D)^2)) %>%
    arrange(pairType)
sse.all = sum(sse$sse)


#   2) Rsq using only E
#   ...................

# function to calculate Rsq of a particular fold
calcRsq = function(predictionFile){
    pred = predictionFile %>%
        readLines %>%
        as.numeric

    Corr = cor(DLong$D, pred) %>% round(2)
    coef = lm(D ~ pred, data = DLong %>% mutate(pred = pred)) %>%
        coefficients

    Rsq.subset = DLong %>%
        mutate(pred = pred) %>%
        group_by(pairType) %>%
        summarize(fit = sum( (D - pred)^2 )) %>%
        right_join(sse, by = 'pairType') %>%
        with(round((1 - fit/sse) * 100, 1)) 
    Rsq.all = DLong %>%
        mutate(pred = pred) %>%
        summarize(fit = sum( (D - pred)^2 )) %>%
        with(round((1 - fit/sse.all) * 100, 1)) 

    out = c(coef, Corr, Rsq.all, Rsq.subset)
    names(out) = c('intercept','slope','corr','Rsq.all','Rsq.inter','Rsq.intra')
    out
}
calcRsq.fold = function(iFold, fitdir){
    file.path(fitdir
            , iFold
            , "prediction.txt"
        ) %>%
        calcRsq
}
calcRsq.allFolds = function(x, folds = 1:100){
    lapply(folds
        , calcRsq.fold
        , fitdir = x
    ) %>%
    do.call(what = rbind) %>%
    as.data.frame %>%
    lapply(quantile, c(0, 0.025, 0.5, 0.975, 1)) %>%
    do.call(what = cbind)
}



#   2) Rsq using only E
#   ...................

# function to get HPD of predictions
calcPred = function(fitdir, folds = 1:100){
    file.path(fitdir
            , folds
            , "prediction.txt"
        ) %>%
        lapply(function(predictionFile){
            predictionFile %>%
            readLines %>%
            as.numeric
        }) %>%
        do.call(what = rbind) %>%
        apply(2, quantile, c(0.025, 0.975)) %>%
        t
}