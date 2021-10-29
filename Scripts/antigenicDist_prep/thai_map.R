#!/usr/bin/env Rscript --vanilla

library(tidyverse); theme_set(theme_minimal())
library(rgl)


outstem = "02-processedData/antigenicDist/thai_map"
dir.create(dirname(outstem), recursive = T)

config = new.env()
source("Scripts/configs/plots.R", local = config)

    #   import 3d antigenic cartography coordinates
    #   from Katzelnick, 2021 (Science)
    #   ...........................................
    
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
mutate(
    Year = (strsplit(Strain,"/") %>% do.call(what=rbind))[ ,3] %>% as.integer
    , Serotype = (strsplit(Strain,"/") %>% do.call(what=rbind))[ ,1]
) %>%
filter(
    # only include viruses (not sera)
    Group == "AG"
    # limit to Thai viruses from 1994 on
    , grepl('THAILAND', Strain)
    , Year >= 1994
) %>%
select(-Group)

# plot 3D scatterplot
open3d()
par3d(windowRect = c(100, 100, 612, 612))
with(a3d, plot3d(
    x, y, z
    , col = config$Colors$Serotype[Serotype]
    , type = "s"
    , radius = 0.2
    , aspect = TRUE
    , xlab = "" , ylab = "" , zlab = ""
))
legend3d("topright"
    , legend = names(config$Colors$Serotype)
    , pch = 16
    , col = config$Colors$Serotype
    , cex=1
    , inset=c(0.02)
)

# play the 3D rotation
#play3d( spin3d( axis = c(0, 0, 1), rpm = 5), duration = 20 )

# save to GIF
movie3d(
    movie = basename(outstem)
    , spin3d( axis = c(0, 0, 1), rpm = 5)
    , duration = 20
    , dir = dirname(outstem)
    , type = "gif"
    , clean = TRUE
)


    #   Convert to long format
    #   ......................
    
a3d.mat = a3d %>%
    select(x, y, z) %>%
    as.matrix
rownames(a3d.mat) = a3d$Strain
DLong = dist(a3d.mat) %>% 
    as.matrix %>%
    as.data.frame %>%
    mutate(virus = a3d$Strain) %>%
    gather(serum, D, -virus) %>%
    filter(virus != serum)

write_csv( DLong
    , paste0(outstem, ".csv")
)




# antigenic distance predictions using just serotype centroids
avgDist = a3d %>%
    group_by(Serotype) %>%
    summarize(
        x = mean(x)
        , y = mean(y)
        , z = mean(z)
    ) %>%
    arrange(Serotype) %>%
    select(-Serotype) %>%
    dist %>%
    as.matrix
rownames(avgDist) = colnames(avgDist) = paste0('DENV',1:4)

DLong %>%
    mutate(
        virusType = substring(virus, 1, 5)
        , serumType = substring(serum, 1, 5)
    ) %>%
    with(mapply( function(v,s){ avgDist[v,s] }
        , v = virusType
        , s = serumType
    )) %>%
    as.character %>%
    writeLines( paste0(outstem,'.dist') )
    