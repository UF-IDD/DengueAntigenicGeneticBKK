#' ---
#' title: Compare identified substitutions with effects beyond E (30 vs 60AA from NS2A)
#' author: "Angkana T. Huang"
#' output:
#'   html_document:
#'      toc: true
#' params:
#'     effectdir: "04-plots/fits_adjusted"
#' ---


#+ include = FALSE
knitr::opts_chunk$set(include = TRUE, echo = FALSE, warning = FALSE)



library(tidyverse); theme_set(theme_minimal())
library(ggpubr)
attach(params)


plotdir = "04-plots/compare_adjusted"
dir.create(plotdir, recursive = T)

source("Scripts/configs/plots.R")


# Import identified substitutions
effectFiles = list.files(effectdir, pattern = "\\_[0-9]+_nonzero.csv$", full.names = T)
vers = basename(effectFiles) %>%
            gsub(pattern = "_nonzero\\.csv$", replacement = "")
subs = lapply(effectFiles, read_csv, col_types = cols(
        Feature = col_character(),
        pos = col_double(),
        subClust = col_double(),
        Median = col_double(),
        lower = col_double(),
        upper = col_double(),
        count = col_double()
    )) %>%
    c(list(
        by = c("Feature","pos")
        , suffix = paste0(".",vers)
    )) %>%
    do.call(what = full_join)

#' 2x2 table of nonzero effect sites in the two versions
#' shows that all nonzero effect sites captured when sampling 30AA at a time
#' were captured when sampling 60AA at time.
pos.nonzero = subs %>%
    group_by(pos) %>%
    summarize(
        present30 = any(!is.na(Median.AA_NS2A_30))
        , present60 = any(!is.na(Median.AA_NS2A_60))
    )
pos.nonzero %>%
    with(table(present30, present60))

#' Comparison of their effect sizes at sites that were captured
#' in both site sampling schemes.
subs %>%
    ggplot(aes(x = Median.AA_NS2A_30, y = Median.AA_NS2A_60))+
    geom_abline(slope = 1, color = "blue", linetype = 2)+
    geom_errorbar(aes(ymin = lower.AA_NS2A_60, ymax = upper.AA_NS2A_60), color = "#dddddd")+
    geom_errorbarh(aes(xmin = lower.AA_NS2A_30, xmax = upper.AA_NS2A_30), color = "#dddddd")+
    geom_point(color = "black")+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank())+
    xlab('Effect size (sampling 30 AA each time)')+
    ylab('Effect size (sampling 60 AA each time)')



#' To be conservative, we rely on sites identified by
#' downsampling NS2A to 60AA (includes less sites).
subs = read_csv(effectFiles[which(vers == "AA_NS2A_60")], show_col_types = FALSE)

# write unique positions to file for later use
subs$pos %>%
    unique %>%
    sort %>%
    as.character %>%
    writeLines(file.path(plotdir,"positions_NS2A_beyondE.txt"))
