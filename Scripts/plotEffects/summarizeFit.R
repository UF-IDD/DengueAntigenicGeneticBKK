#!/usr/bin/env Rscript --vanilla


library(tidyverse)
library(parallel)

inputArgs = commandArgs(trailing = T)
fitdir = inputArgs[1]

ncpu = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
ncpu = ifelse(is.na(ncpu), 1, ncpu)


    #   Effects

if(Sys.getenv("NOT_SAVE_ESTIMATES") != "TRUE"){
    # identify files storing estimated effects "d.csv"
    dFiles = list.files(fitdir
        , pattern = "d\\.csv$"
        , recursive = T
        , full.names = T
    )

    # derive metadata from file path
    meta = strsplit(dFiles, "/") %>%
        do.call(what = rbind)
    meta = meta[ , ncol(meta) - (7:1) ]
    meta = data.frame(
        fitmethod = meta[ ,1]
        , antigenicdist = meta[ ,2]
        , aaver = meta[ ,3]
        , adjver = meta[ ,4]
        , featureType = meta[ ,5]
        , protein_name = meta[ ,6]
        , iFold = meta[ ,7] %>% as.integer
    )

    cl = makeCluster(ncpu)

    # process and write summarized effects
    d = clusterMap(cl, function(m,f){
            library(tidyverse)
            # read in the effects
            m$d = lapply(f, read_csv, col_types = list(
                col_character() # Feature
                , col_number() # effect
                , col_integer() # clustNum
                , col_integer() # pos
            ))
            # summarize
            m %>%
                unnest(cols = d) %>%
                group_by(protein_name, Feature, pos) %>%
                summarize(
                    lower = quantile(effect, 0.025)
                    , Median = median(effect)
                    , upper = quantile(effect, 0.975)
                    , Min = min(effect)
                    , Max = max(effect)
                    , n_nonzero = sum(effect > 0)
                    , n_total = n()
                )
        }
        , split(meta, meta %>% select(protein_name))
        , split(dFiles, meta %>% select(protein_name))
        , SIMPLIFY = F
    ) %>%
    do.call(what = rbind)

    d %>%
        write_csv(file.path(fitdir, "d_summary.csv"))
    
    d %>%
        # retain only entries with non-zero effects
        filter(lower > 0) %>%
        arrange(protein_name, desc(Median)) %>%
        write_csv(file.path(fitdir, "d_nonzero.csv"))
        
    # get uncertainties around expected (v+p)
    vFiles = list.files(fitdir
        , pattern = "v\\.csv$"
        , recursive = T
        , full.names = T
    )
    pFiles = list.files(fitdir
        , pattern = "p\\.csv$"
        , recursive = T
        , full.names = T
    )
    # distances arising from measurement uncertainties 
    dm = list.files('03-fits/noisyCoords', full.names = T) %>%
        lapply(read_csv, col_types = cols(
          serotype = col_character(),
          `0%` = col_double(),
          `2.5%` = col_double(),
          `25%` = col_double(),
          `50%` = col_double(),
          `75%` = col_double(),
          `97.5%` = col_double(),
          `100%` = col_double(),
          Mean = col_double(),
          Sd = col_double()
        )) %>%
        do.call(what = rbind) %>%
        filter(serotype == 'All') %>%
        arrange(Mean) %>%
        mutate(i = seq_along(Mean))
    
    plotoutdir = file.path(fitdir, 'intercept')
    dir.create(plotoutdir, recursive = T)
    clusterMap(cl, function(vf,pf,od,dm){
            library(tidyverse); theme_set(theme_classic())
            library(ggpubr)
            source('Scripts/configs/plots.R')
            # meta data
            meta = strsplit(vf, "/") %>% do.call(what = rbind)
            meta = meta[ , ncol(meta) - (7:1) ]
            protein_name = meta[ ,6][1]
            meta = data.frame(iFold = meta[ ,7] %>% as.integer)
            
            # read in the effects
            meta$v = lapply(vf, read_csv, col_types = list(
                col_character() # virus 1
                , col_number() # intercept
            ))
            meta$p = lapply(pf, read_csv, col_types = list(
                col_character() # virus 2
                , col_number() # intercept
            ))
            # plot
            v = meta %>% select(iFold, v) %>% unnest(cols = v)
            p = meta %>% select(iFold, p) %>% unnest(cols = p)
            getYear = function(x){
                x = strsplit(x,'/') %>% do.call(what = rbind)
                as.integer(x[ ,3])
            }
            vp = left_join(v, p, by = c('iFold', 'virus' = 'serum')) %>%
                rename(v = avidity, p = potency) %>%
                mutate(Year = getYear(virus))
                
            vp.avg = vp %>%
                group_by(iFold) %>%
                summarize(mid = mean((v+p)/2))
            g1 = vp %>%
                group_by(virus, Year) %>%
                summarize(
                    vp = (v+p)/2
                    , lower = min(vp), mid = mean(vp), upper = max(vp)
                ) %>%
                ungroup %>%
                arrange(mid) %>%
                mutate(
                    Serotype = substring(virus, 1, 5)
                    , virus = factor(virus, levels = unique(virus))
                ) %>%
                ggplot(aes(x = virus, ymin = lower, ymax = upper
#                     , color = Serotype
                    , color = Year
                ))+
                geom_hline(data = vp.avg
                    , aes(yintercept = mid), color = 'grey'
                )+
                geom_linerange(size = 0.5)+
#                 scale_color_manual(values = Colors$Serotype)+
                scale_color_gradient2(
                    low = '#202585'
                    , mid = '#bce6cf'
                    , high = 'red'
                    , midpoint = 2005
                )+
                ylim(c(0,4))+
                xlab('Virus')+
                theme(
                    axis.ticks.x = element_blank()
                    , axis.text = element_blank()
                    , axis.title.y = element_blank()
                    , legend.position = c(0.05, 0.95)
                    , legend.justification = c(0,1)
                )
            g2 = vp %>%
                ggplot(aes(y = (v+p)/2, group = iFold))+
                geom_hline(data = vp.avg
                    , aes(yintercept = mid), color = 'grey'
                )+
                geom_density()+
                scale_x_reverse()+
                scale_y_continuous('Virus-specific intercepts'
                    , limits = c(0,4)
                    , position = 'right'
                )+
                xlab('Density')+
                theme(
                    axis.ticks.x = element_blank()
                    , axis.text.x = element_blank()
                    , axis.title.y.right = element_text(angle = 90)
                )
            g3 = dm %>%
                ggplot(aes(x = i))+
                geom_boxplot(aes(
                    group = i
                    , ymin = `2.5%`/2
                    , lower = `25%`/2
                    , middle = `50%`/2
                    , upper = `75%`/2
                    , ymax = `97.5%`/2
                ), stat = 'identity', color = '#ebe660', fill = '#999999', size = 0.1)+
                geom_step(aes(y = `50%`/2), direction = 'h', color = 'orange')+
                geom_step(aes(y = Mean/2), direction = 'h')+
                ylab('Distance from self / 2')+
                xlab('Synthetic\ndataset')+
                ylim(c(0,4))+
                theme(
                    axis.ticks.x = element_blank()
                    , axis.text.x = element_blank()
                )
            g = ggarrange(g2, g1, g3
                , ncol = 3
                , nrow = 1
                , widths = c(1,2,1)
                , align = 'h'
                , labels = 'auto'
                , vjust = 1
                , hjust = c(-0.5, -1.5, -0.5)
            )
            ggsave(filename = file.path(od, paste0(protein_name,'.pdf'))
                , width = 6
                , height = 4.5
            )
            vp.avg %>%
                summarize(
                    Mean = mean(mid)
                    , lower = quantile(mid, 0.025)
                    , upper = quantile(mid, 0.975)
                )
        }
        , split(vFiles, meta %>% select(protein_name))[1]
        , split(pFiles, meta %>% select(protein_name))[1]
        , od = plotoutdir
        , dm = list(dm)
        , SIMPLIFY = F
    ) 
    
}

    #   RMSE
    
rmseFiles = list.files(fitdir
    , pattern = "rmse\\.csv$"
    , recursive = T
    , full.names = T
)

rmse = lapply(rmseFiles, read_csv, col_types = list(
    col_character() # featureType
    , col_character() # protein_name
    , col_integer() # fold
    , col_number() # train_rmse
    , col_number() # test_rmse
))

rmse %>%
    do.call(what = rbind) %>%
    gather(dataset, rmse, ends_with("_rmse")) %>%
    group_by(featureType, protein_name, dataset) %>%
    summarize(
        lower = quantile(rmse, 0.025)
        , Median = median(rmse)
        , upper = quantile(rmse, 0.975)
    ) %>%
    write_csv(file.path(fitdir, "rmse_summary.csv"))
    