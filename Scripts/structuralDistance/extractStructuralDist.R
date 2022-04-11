
library(tidyverse); theme_set(theme_classic())
library(bio3d)


# Structure (PDB object) downloaded from
# : https://www.rcsb.org/structure/3j27
pdb = read.pdb("00-RawData/structure/3j27.pdb")


# https://rdrr.io/cran/bio3d/man/dm.html
# function to extract structural distance from PDB object


# select one atom per residue
x = dm(pdb
    , inds = atom.select(pdb, "protein"
        , eleno = pdb$atom %>%
            filter(resno %in% 1:495) %>%
            group_by(resno) %>%
            filter(eleno == min(eleno)) %>%
            with(eleno)
        , grp = T
        , grpby = pdb$atom[,"resno"]
        , verbose=T
    )
)
image(x, axes = F)
lapply(1:2, function(i){
    axis(i
        , at = seq(0, 1, length.out = 495/15 + 1)
        , labels = seq(0, 495, by = 15)
    )    
})




# process one residue at a time to allow closest atom of other residues
Dist = lapply(1:494, function(i){
    lapply((i+1):495, function(j){
        if(i==j){ return(NULL) }
        eleno.self = pdb$atom %>%
            filter(resno %in% i) %>%
            with(min(eleno))
        Atoms = atom.select(pdb, "protein"
            , eleno = c(
                eleno.self
                , pdb$atom %>%
                    filter(resno == j) %>%
                    with(eleno)
            )
        )
        Dist = dm(pdb, inds =  Atoms)[1, ] %>% min(na.rm = T)
        data.frame(i, j, Dist)
    }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)

# plot to check
Dist %>%
    ggplot(aes(x = i, y = j, fill = Dist))+
    geom_tile()+
    scale_fill_gradient2("Distance (angstrom)"
        , low = "red"
        , mid = "#dddddd"
        , high = "black"
        , midpoint = 40
    )+
    coord_fixed(ratio = 1)

# show that without accounting for proximity to other dimers, proximity may be over estimated
Dist %>%
    filter(i==1) %>%
    arrange(j) %>%
    mutate(D = x[1,-1]) %>%
    ggplot(aes(x = Dist, y = D))+
    geom_point()+
    xlab('Distance, All dimers considered')+
    ylab('Distance, Residue grouping')


outdir = "02-processedData/structure"
dir.create(outdir, recursive = T)

Dist %>%
    write_csv(file.path(outdir, "ResidueDist.csv"))

