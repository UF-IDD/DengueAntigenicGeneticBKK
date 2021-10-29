# Load the Racmacs package
library(Racmacs)
library(tidyverse)

iSet = Sys.getenv("SLURM_ARRAY_TASK_ID")

outdir = '03-fits/noisyCoords'
dir.create(outdir, recursive = T)

# Read in the titer table
titerFile = paste0('02-processedData/noisyTiter/', iSet,'.txt')
titer_table = read.titerTable(titerFile)


# Create the acmap object, specifying the titer table
map <- acmap(
  titer_table = titer_table
)

# Perform some optimization runs on the map object to try and determine a best map
map <- optimizeMap(
  map                     = map,
  number_of_dimensions    = 3,
  number_of_optimizations = 10000,
  minimum_column_basis    = "none",
  options = list(num_cores = Sys.getenv("SLURM_CPUS_ON_NODE") %>% as.integer)
)


# Extract distances between noisy entries of the same virus
iNoise = grep('_noise[0-9]+$', rownames(titer_table))
Noise = split(iNoise, gsub('_noise[0-9]+$', '', rownames(titer_table)[iNoise])) %>%
    lapply(function(i){
        d = dist(agCoords(map)[i, ]) %>% as.matrix
        d[lower.tri(d)]
    }) %>%
    do.call(what = c)


# Summarize the noise
Summary = function(x){
    c(
        quantile(x, c(0,0.025, 0.25, 0.5, 0.75, 0.975, 1))
        , Mean = mean(x)
        , Sd = sd(x)
    )
}

# by serotype
noiseType = by(Noise, substring(names(Noise), 1, 5), Summary)

data.frame(serotype = names(noiseType)) %>%
cbind(do.call(rbind, noiseType)) %>%
rbind(c('All', Summary(Noise))) %>% # All serotypes combined
write_csv(file.path(outdir, paste0(iSet,'.csv')))

