

library(tidyverse); theme_set(theme_classic())

# import titer table
titer = read.table('00-RawData/titers/MAY15-2020-combined-set1set2set3set4set5set6set7-sciencem3-exp-conditions-reasonable-adjust.txt')


# import coordinates inferred from these titers
a3d = read.table("00-RawData/cartography/reasonable-adjust-3d-coords.txt"
    , col.names = c(
        "Group"
        , "Strain"
        , "x"
        , "y"
        , "z"
    )
    , stringsAsFactors = F
) %>%
filter(Group == 'AG')


genTiterNoise = function(x, log10Variance = 0.13){
    if(x == '*'){ return(x) }
    below = grepl('^<', x)
    xvalue = str_extract(x, '[0-9]+') %>% as.integer
    out = xvalue * 10^rnorm(1, sd = sqrt(log10Variance))
    if(below & (out<xvalue)){ return(x) }
    return(round(out))
}
genTiterNoiseVector = function(x){
    sapply(x, genTiterNoise, USE.NAMES = F)
}



# generate titer datsets with noise
outdir = '02-processedData/noisyTiter'
dir.create(outdir, recursive = T)

set.seed(817919)
lapply(1:100, function(iSet){
    # sample 4x4 viruses covering varies areas of the antigenic space
    v = split(a3d, substring(a3d$Strain, 1, 5)) %>%
    lapply(function(a){
        by(a, a %>%
            mutate(far =
                (x > median(x)) &
                (y > median(y)) &
                (z > median(z))
            ) %>%
            select(far)
            , sample_n
            , size = 2
        ) %>%
        do.call(what = rbind)
    }) %>%
    do.call(what = rbind) %>%
    with(Strain)

    # generate noisy titer entries for these viruses
    noisy = lapply(1:4, function(iNoise){
        out = titer[v, ] %>% apply(2, genTiterNoiseVector)
        rownames(out) = paste(rownames(out), iNoise, sep = '_noise')
        out
    }) %>%
    do.call(what = rbind)

    titer[ setdiff(rownames(titer), v), ] %>%
        rbind(noisy) %>%
        write.table(file.path(outdir, paste0(iSet,'.txt')))
})
