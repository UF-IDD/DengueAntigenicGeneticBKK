
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

outdir = "05-coevolution/fastcov"
dir.create(outdir, recursive = T)


posNS2A = readLines("02-processedData/AA_NS2A_beyondE/subset.txt") %>% as.integer
eff = read_csv("04-plots/fits_adjusted/AA_NS2A_beyondE_nonzero.csv", col_types = "ciidddi")
eff.E = read_csv("03-fits/nnls/thai_map/AA/adj_none/aasub/d_summary.csv", col_types = "ccidddddii")


# load database of known antigenically relevant substitutions
dat = new.env()
source("Scripts/antibodyDB_prep/db_load.R", local = dat, verbose = F)
aaver = "AA"
source("Scripts/sequence_load/aaver.R", local = dat, verbose = F)


# define functions that given counts observed AA at the position
varaa = new.env()
varaa$wukabat = function(x){
    (sum(x) * length(x)) / max(x)
}

# select non-overlapping proteins that forms the polyprotein
protein_i = c(
    # 1   # anchored capsid
    # , 5 # prM
    'E' = 3 # envelope
    # , 6 # NS1
    , 'NS2A' = 7 # NS2A
    # , 8 # NS2B
    # , 9 # NS3
    # , 10 # NS4A
    # , 12 # 2K
    # , 11 # NS4B
    # , 14 # NS5
)
variations = lapply(protein_i, function(i){
    tibble(
        wukabat = do.call(rbind, dat$aa[[i]]) %>%
            apply(2, function(aai){
                counts = table(aai)
                varaa$wukabat(counts)
            })
    ) %>%
    mutate(pos = seq_along(wukabat))
})





indir = "05-coevolution/fastcov"
f = "nonstructural protein NS2A.pairs.txt"
coev = file.path(indir, f) %>%
    read_tsv(col_types = "cccdd") %>%
    mutate(
          posX = str_extract(SiteX, '^[0-9]+') %>% as.integer
        , posY = str_extract(SiteY, '^[0-9]+') %>% as.integer
        , X = ifelse(posX <= 495, 'E', 'NS2A')
        , Y = ifelse(posY <= 495, 'E', 'NS2A')
    ) %>%
    filter(X == "E", Y != "E") %>%
    mutate(posY = posY - 495) %>%
    mutate(res = substring(Pair,2))
    
# eff %>%
#     left_join(coev, by = c(
#         'pos' = 'posY'
#     )) %>%
#     filter(res == res1 | res == res2) %>%
#     filter(pos==2)



gVar = lapply(variations, function(x){
    g = ggplot(x, aes(x = pos, y = wukabat))+
        geom_col(width = 1)+
        scale_x_continuous(expand = c(0,0), limits = c(0.5, max(x$pos)+0.5))+
        ylab('Wukabat\ncoef.')
    return(g)
})




# measures = c(
#     'shannon' = "Shannon entropy"
#     , 'simpson' = "Simpson diversity index"
#     , 'wukabat' = "Wu-kabat variability coefficient"
# )


####################

voidPlot = ggplot()+theme_void()+theme(
    plot.margin = margin(-10,-10,-10,-10, 'lines')
    , axis.title = element_blank()
)
#library(patchwork)

gCoev = coev %>%
    ggplot(aes(x = posY, y = posX))+
    geom_density_2d_filled(contour_var = "ndensity", adjust = 1/4)+
    scale_y_continuous(expand = c(0,0)
        , breaks = seq(10,495, by= 20)
        , position = 'right'
    )+
    scale_x_continuous(expand = c(0,0)
        , limits = c(1,219)
        , position = 'top'
        , breaks = seq(10,210, by= 10)
    )+
    # scale_x_continuous(expand = c(0,0), position = 'top')+
    xlab('AA position in NS2A')+
    ylab('AA position in E')+
    theme(
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0)
    )

g = ggarrange(
    gVar[['NS2A']] + 
        scale_y_continuous(expand = c(0,0), position = 'right')+
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
        # scale_y_continuous(position = 'right')
    , voidPlot
    , voidPlot
    , eff %>%
        ggplot(aes(x= pos))+
        # geom_histogram(binwidth = 1, fill = 'red')+
        geom_bar(width = 1, fill = 'red')+
        scale_x_continuous(expand = c(0,0), limits = c(0.5,219.5))+
        scale_y_continuous(expand = c(0,0), position = 'right')+
        ylab('Num.\nsubs.')+
        theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    , voidPlot
    , voidPlot
    , gCoev
    , eff.E %>%
        ggplot(aes(x= pos))+
        geom_histogram(binwidth = 1, fill = 'red')+
        scale_x_continuous(expand = c(0,0), limits = c(1,495))+
        scale_y_continuous(expand = c(0,0), position = 'right')+
        ylab('Num.\nsubs.')+
        theme(axis.title.y = element_blank(), axis.text.y = element_blank())+
        coord_flip()
    , gVar[['E']]+coord_flip() +
        scale_y_continuous(expand = c(0,0), position = 'right')+
        theme(
            axis.title.y = element_blank()
            , axis.text.y = element_blank()
            ,     plot.margin = margin(-10,-10,-10,-10, 'lines')
        )
    , nrow = 3
    , ncol = 3
    , align = 'hv'
    , heights = c(2,2,6)
    , widths = c(6,2,2)
    , common.legend = T
    , legend = "right"
) 
ggsave(g
    , filename = file.path(outdir, paste0("coev_", "hv",".pdf"))
    , height = 7
    , width = 9
)
