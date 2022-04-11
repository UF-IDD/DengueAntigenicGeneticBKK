#' ---
#' title: Effects of substitutions adjusted for E
#' author: "Angkana T. Huang"
#' output:
#'   html_document:
#'      toc: true
#' params:
#'     fitdir: ""
#'     nulldir: ""
#'     aadir: ""
#'     lengthAdj: 495
#'     exportAdjuster: false
#' ---


#+ include = FALSE
knitr::opts_chunk$set(include = TRUE, echo = FALSE, warning = FALSE)

library(tidyverse); theme_set(theme_minimal())
library(ggpubr)

#' # Inputs
params
attach(params)

plotdir = "04-plots/fits_adjusted"
dir.create(plotdir, recursive = T)

aaver = basename(aadir)

source("Scripts/configs/plots.R")


dfit = read_csv(file.path(fitdir, "d_nonzero.csv")) %>%
    mutate(pos = pos - lengthAdj) %>%
    rename(iSample = protein_name) %>%
    mutate(subClust = str_extract(Feature, "[0-9]+") %>% as.integer)

if(exportAdjuster){
    # effects of substitutions on the adjuster protein
    fit.adjuster = dfit %>%
        # positions in adjuster (E)
        filter(pos <= 0) %>%
        mutate(pos = pos + lengthAdj) %>%
        mutate(
            Feature = paste0(
                substring(Feature,1,1)
                , pos
                , str_extract(Feature, "[a-z]$")
            )
        )
    
    # plot
    g = fit.adjuster %>%
        ggplot(aes(x = pos, y = Median, ymin = lower, ymax = upper))+
        geom_errorbar()+
        geom_point()+
        labs(
            x = "AA position"
            , y = "Effect size"
            , alpha = "Number of estimates"
        )+
        ggtitle("Effect size of adjuster protein")
    plot(g)

    # export for future use
    fit.adjuster %>%
        arrange(desc(Median)) %>%
        write_csv(file.path(plotdir, paste0(aaver,"_nonzero_adjuster.csv")))
}



dfit = dfit %>%
    # remove positions from E
    filter(pos > 0) %>%
    # remove positions that were collinear with E
    anti_join(
        dfit %>%
            filter(pos <= 0) %>%
            select(iSample, subClust)
        , by = c("iSample","subClust")
    )

rmse.fit = read_csv(file.path(fitdir, "rmse_summary.csv")) %>%
    filter(dataset == "test_rmse") %>%
    rename(iSample = protein_name)

rmse.null = read_csv(file.path(nulldir, "rmse_summary.csv")) %>%
    filter(dataset == "test_rmse") %>%
    rename(iSample = protein_name)



# convert back to actual position on protein of interest
posActual = dfit$iSample %>%
    unique %>%
    lapply(function(iSample){
        readLines(paste0(aadir,iSample,".txt")) %>%
        as.integer
    })
names(posActual) = dfit$iSample %>% unique

# how often were each position sampled?
posFreq = posActual %>%
    unlist %>%
    table

dfit = dfit %>%
    mutate(
        pos = mapply(function(iSam, p){
                posActual[[as.character(iSam)]][p]
            }
            , iSam = iSample
            , p = pos
        )
    )




#' How frequent were these positions estimated to have non-zero effect size
#' given they are included in the analysis
posFreqDf = data.frame(pos = as.integer(names(posFreq)), nSampled = as.integer(posFreq)) %>%
    left_join(
        dfit %>%
            select(iSample, pos) %>%
            unique %>%
            group_by(pos) %>%
            summarize(nEffect = n())
        , by = "pos"
    ) %>%
    mutate(nEffect = ifelse(is.na(nEffect), 0, nEffect)) %>%
    mutate(freqEffect = nEffect/nSampled)

ggarrange(
    posFreqDf %>%
        ggplot(aes(x = nSampled, y = freqEffect))+
        geom_point(alpha = 0.2, size = 3)+
        theme_classic()+
        xlab("Number of times sampled")+
        ylab("Frequency of having\nnon-zero effect")+
        ylim(c(0,max(posFreqDf$freqEffect)))
    , posFreqDf %>%
        ggplot(aes(x = pos, y = freqEffect))+
        geom_col(width = 1)+
        theme_classic()+
        ylim(c(0,max(posFreqDf$freqEffect)))
    , ncol = 2
    , nrow = 1
    , align = "h"
)

hist(posFreqDf$freqEffect, 100, col = "black"
    , xlab = "Frequency of having non-zero effect when included"
    , ylab = "Number of positions"
    , main = NULL
)


#' Number of nonzero effect substitutions identified across a range of significance thresholds.
dfit %>%
    mutate(prop_nonzero = n_nonzero / n_total) %>%   
    arrange(desc(prop_nonzero)) %>%
    mutate(i = seq_along(prop_nonzero)) %>%
    ggplot(aes(x = i, y = prop_nonzero, color = lower > 0))+
    geom_step(direction = 'h')+
    scale_color_manual(values = c('FALSE'='#444444', 'TRUE'='red'), guide = 'none')+
    xlab('Number of substitutions')+
    ylab('Proportion of estimations\nshowing nonzero effect')+
    ylim(c(0,1))+
    theme_classic()
ggsave(
    filename = file.path(plotdir, paste0('significanceThresh_',aaver,'.pdf'))
    , width = 5
    , height = 3
)
  


#' RMSE: sampling AA from NS2A vs null
rmseRange = range(rmse.null$upper, rmse.fit$lower)
rmse.null = rmse.null %>%
    arrange(lower) %>%
    mutate(iSample = factor(iSample, levels = unique(iSample)) %>% as.integer)
rmse.fit = rmse.fit %>%
    arrange(lower) %>%
    mutate(iSample = factor(iSample, levels = unique(iSample)) %>% as.integer)

ggplot(mapping = aes(x = iSample))+
    # null
    geom_errorbar(data = rmse.null, aes(ymin = lower, ymax = upper), color = "grey")+
    geom_point(data = rmse.null, aes(y = Median, color = "Other proteins"))+
    # fit
    geom_errorbar(data = rmse.fit, aes(ymin = lower, ymax = upper), color = "yellow", alpha = 0.5)+
    geom_point(data = rmse.fit, aes(y = Median, color = "NS2A"))+
    ylab("RMSE")+
    ylim(rmseRange)+
    theme(
        panel.grid = element_blank()
        , axis.text.x = element_blank()
    )+
    labs(color = "AA sampled from")+
    scale_color_manual(
        values = c("NS2A" = "red", "Other proteins" = "black")
    )


#' Effect sizes of substitutions at positions which showed non-zero effect >99% of the times they were included
dfit.99 = dfit %>%
    semi_join(  posFreqDf %>% filter(freqEffect > .99)
        , by = "pos"
    ) %>%
    mutate(
        Feature = paste0(
            substring(Feature,1,1)
            , pos
            , str_extract(Feature, "[a-z]$")
        )
    ) %>%
    group_by(Feature, pos, subClust) %>%
    summarize(
        Median = median(Median)
        , lower = min(lower)
        , upper = max(upper)
        , count = n()
    ) 
dfit.99 %>%
    ggplot(aes(x = pos, y = Median, ymin = lower, ymax = upper, alpha = count))+
    geom_errorbar()+
    geom_point()+
    labs(
        x = "AA position"
        , y = "Effect size"
        , alpha = "Number of estimates"
    )

# Interactive table of these substitutions
DT::datatable(dfit.99)

# export for future use
posFreqDf %>%
    write_csv(file.path(plotdir, paste0(aaver,"_freqNonzero.csv")))
dfit.99 %>%
    arrange(desc(Median)) %>%
    write_csv(file.path(plotdir, paste0(aaver,"_nonzero.csv")))



posPresence = data.frame(iSample = seq_along(posActual))
posPresence$pos = posActual
posPresence = posPresence %>% unnest(cols = pos)
posPresence = by(posPresence, posPresence$pos, function(x){
    rmse.fit %>%
        mutate(included = iSample %in% x$iSample) %>%
        group_by(included) %>%
        summarize(
            Median = median(Median)
            , lower = quantile(Median, 0.025)
            , upper = quantile(Median, 0.975)
        ) %>%
        mutate(pos = x$pos[1])
}) %>%
do.call(what = rbind)
posPresence = posPresence %>%
    filter(included) %>%
    left_join( posPresence %>% filter(!included)
        , by = "pos"
        , suffix = c('.include','.exclude')
    )

# if not all positions are included
if(!all(is.na(posPresence$included.exclude))){
    posPresence %>%
        ggplot(aes(x = Median.include, y = Median.exclude))+
        geom_abline(slope = 1, intercept = 0, color = "blue", lty = 2)+
        geom_errorbar(aes(ymin = lower.exclude, ymax = upper.exclude))+
        geom_errorbarh(aes(xmin = lower.include, xmax = upper.include))+
        geom_point(size = 0.2)+
        coord_fixed(ratio = 1)+
        xlab("RMSE when included")+
        ylab("RMSE when excluded")

    posPresence %>%
        ggplot(aes(x = pos, y = Median.include - Median.exclude))+
        geom_col(width = 1)+
        ylab("RMSE include - exclude")

    posPresence %>%
        ggplot(aes(x = Median.include - Median.exclude))+
        geom_histogram(bins = 60)+
        xlab("RMSE include - exclude")+
        ylab("# positions")
}