#' ---
#' title: Relevant substitutions beyond E
#' author: "Angkana T. Huang"
#' output:
#'   html_document:
#'      toc: true
#' params:
#'     plotdir: "04-plots/beyondE"
#' ---


#+ include = FALSE
knitr::opts_chunk$set(include = TRUE, echo = FALSE, warning = FALSE)

library(tidyverse); theme_set(theme_minimal())
library(ggpubr)
attach(params)

# plot setups
dir.create(plotdir, recursive = T)
source("Scripts/configs/plots.R")


# extract protein lengths
dat = new.env()
aaver = "AA"
source("Scripts/sequence_load/aaver.R", local = dat, verbose = F)
proteinLengths = sapply(dat$aa, function(x){
    x[[1]] %>% length
}) %>%
sort


# load short names for the proteins
x = read_csv("02-processedData/protein_shortname.csv", show_col_types = FALSE)
short = x$short
names(short) = x$protein_name



#'
#' # Protein with excess signals
#'

    #   Import RMSE
    #   ...........
    
# Actual protein
fitdir = "03-fits/nnls/thai_map/AA/adj_none/aasub"
rmse = read_csv(file.path(fitdir, "/rmse_summary.csv"), show_col_types = FALSE) %>%
    mutate(
        i = 0
        , aaSource = "Actual"
    )

# Random sites on the polyprotein (aside NS2A and E)
rmseRandomPoly = lapply(names(proteinLengths), function(pname){
    pLength = proteinLengths[[pname]]
    paste0(
        "03-fits/nnls/thai_map/AA_random2length_"
        , pLength
        , "/adj_none/aasub/rmse_summary.csv"
    ) %>%
    read_csv(show_col_types = FALSE) %>%
    arrange(Median) %>%
    mutate(
        i = seq_along(protein_name)
        , protein_name = pname
        , aaSource = "Polyprotein"
    )
}) %>%
do.call(what = rbind)

# Random sites on the E protein
rmseRandomE = lapply(names(proteinLengths), function(pname){
    pLength = proteinLengths[[pname]]
    indir = "03-fits/nnls/thai_map/AA_sample2length_"
    I = paste0(indir, pLength) %>% list.files

    lapply(I, function(i){
        paste0(
            indir
            , pLength
            , "/"
            , i
            , "/adj_none/aasub/rmse_summary.csv"
        ) %>%
        read_csv(show_col_types = FALSE) %>%
        arrange(Median) %>%
        mutate(
            i = i
            , protein_name = pname
            , aaSource = "Envelope"
        )
    }) %>%
    do.call(what = rbind)
}) %>%
do.call(what = rbind)



    #   Compare signals against the null comparators
    #   ............................................

legendTitle = 'Origin of sites'
rmseExcess = list(rmseRandomPoly, rmseRandomE, rmse) %>%
    Reduce(f = rbind) %>%
    filter(dataset == "test_rmse") %>%
    mutate(
        protein_name = factor(protein_name, levels = names(proteinLengths))
        , aaSource = factor(aaSource, levels = unique(aaSource))
    ) %>%
    group_by(protein_name) %>%
    mutate(
        i = seq_along(i)
        , i = as.integer(protein_name) + (i/(n()*2)+ 0.25)
    )
gExcess = rmseExcess %>%
    ggplot(aes(
        x = Median, y = i
        , color = aaSource
        , size = aaSource
    ))+
    geom_hline(yintercept = seq_along(dat$aa), size = 0.7, linetype = 1)+
    geom_linerange(aes(xmin = lower, xmax = upper), alpha = 0.7)+
    geom_point()+
    annotate("text"
        , x = 0.9
        , y = seq_along(proteinLengths) + 0.4
        , label = paste0(
            short[names(proteinLengths)]
            , ' (', proteinLengths, ')'
        )
        , hjust = 1
        , vjust = 0#-0.1
    )+
    scale_size_manual(legendTitle, values = c(0.2, 0.2, 1.2))+
    scale_color_manual(legendTitle, values = Colors$aaSource)+
    xlab('RMSE')+
    theme(
        panel.grid = element_blank()
        , axis.text.y = element_blank()
        , axis.title.y = element_blank()
        , axis.ticks.x = element_line(size = 1)
        , legend.position = c(0, 0.1)
        , legend.justification = c(0,0)
        , legend.background = element_rect(fill = "white")
    )




#' Comparing actual performance to the best performing comparator of each set
#' reveals NS2A having more predictive performance than expected.
avgPerformance = function(x){
    x %>%
        filter(dataset == "test_rmse") %>%
        group_by(protein_name) %>%
        summarize(
            lower = mean(lower)
            , upper = mean(upper)
            , Median = mean(Median)
        )
}
lapply( list(
        "polyprotein" = avgPerformance(rmseRandomPoly)
        , "E" = avgPerformance(rmseRandomE)
    ), function(comparator){
        rmse %>%
            filter(dataset == "test_rmse") %>%
            select(protein_name, Median, lower, upper) %>%
            left_join(comparator
                , by = 'protein_name'
                , suffix = c('', '.comparator')
            ) %>%
            mutate(
                deltaMedian = Median - Median.comparator
                , deltaMin = upper - lower.comparator
                , deltaMax = lower - upper.comparator
            ) %>%
            select(protein_name, starts_with('delta')) %>%
            mutate(lowerErrorThanComparator = deltaMin < 0)
    }
)




#'
#' # Identifying sites with excess signals on the NS2A
#'


# Domains of NS2A
domain = read_csv("00-RawData/domains/protein_NS2A.csv", comment = "#", na = "NA"
        , show_col_types = FALSE
    ) %>%
    # fill in names for domains without names
    mutate(domain = ifelse(domain=="", paste0('unnamed',seq_along(domain)), domain)) %>%
    # account for extra residue of DENV-4 at AA position 189
    mutate(
        posStart = ifelse(posStart > 189, posStart+1, posStart)
        , posEnd = ifelse(posEnd >= 189, posEnd+1, posEnd)
    )
domain$pos = mapply( seq, domain$posStart, domain$posEnd, SIMPLIFY = F)
lenMax = max(domain$posEnd)

# adjusted effect sizes
d = "04-plots/fits_adjusted/AA_NS2A_beyondE_nonzero.csv" %>%
    read_csv(show_col_types = FALSE) %>%
    left_join( domain %>%
        select(-posStart, -posEnd) %>%
        unnest(cols = pos)
        , by = "pos"
    )

Colors$location = c(
    "ER lumen" = "#27913e" # green
    , "transmembrane" = "#ffe100" # yellow
    , "cytosol" = "#1c9dff" # blue
)


plotSubsetRmse = function(nAA){
    fitdir = paste0('03-fits/nnls/thai_map/AA_NS2A_',nAA,'/adj_envelope/aasub')
    nulldir = paste0('03-fits/nnls/thai_map/AA_polyNonE_',nAA,'/adj_envelope/aasub')

    rmse.fit = read_csv(file.path(fitdir, "rmse_summary.csv")
            , show_col_types = FALSE
        ) %>%
        filter(dataset == "test_rmse") %>%
        arrange(Median, lower, upper) %>%
        mutate(nAA = nAA) %>%
        mutate(iSample = seq_along(nAA))

    rmse.null = read_csv(file.path(nulldir, "rmse_summary.csv")
            , show_col_types = FALSE
        ) %>%
        filter(dataset == "test_rmse") %>%
        arrange(Median, lower, upper) %>%
        mutate(nAA = nAA) %>%
        mutate(iSample = seq_along(nAA))

    legendTitle = paste(nAA, "sites sampled from")
    sourceColors = c("NS2A" = "red", "Other proteins" = "grey")
    ggplot(mapping = aes(x = iSample))+
        # null
        geom_ribbon(data = rmse.null, aes(ymin = lower, ymax = upper, fill = "Other proteins"), color = NA, alpha = 0.5)+
        geom_line(data = rmse.null, aes(y = Median, color = "Other proteins"))+
        # fit
        geom_ribbon(data = rmse.fit, aes(ymin = lower, ymax = upper, fill = "NS2A"), alpha = 0.3, color = NA)+
        geom_line(data = rmse.fit, aes(y = Median, color = "NS2A"))+
        scale_x_continuous("Subsample", expand = c(0,0))+
        ylab("RMSE")+
        theme(
            panel.grid = element_blank()
            , axis.text.x = element_blank()
        )+
        scale_color_manual(legendTitle
            , values = sourceColors
        )+
        scale_fill_manual(legendTitle
            , values = sourceColors
        )+
        theme_classic()+
        theme(
            legend.position = c(1,0.05)
            , legend.justification = c(1,0)
        )

}


#' Percentile in number of times sites were sampled
#' (when 30AA sampled each time, vs 60AA).
#' Comparing sites results from both sampling schemes revealed 62 sites
#' which consistently showed excess signals.
posFreq = lapply(c(30,60), function(nAA){
    paste0("04-plots/fits_adjusted/AA_NS2A_",nAA,"_freqNonzero.csv") %>%
    read_csv(show_col_types = FALSE) %>%
    mutate(nAA = nAA)
})
gFreq = posFreq %>%
    do.call(what = rbind) %>%
    ggplot(aes(x = pos, y = freqEffect, shape = factor(nAA)))+
    geom_line(aes(group = pos), size = 0.5, color = "grey")+
    geom_point(size =2)+
    scale_shape_manual("Number of sites\nper subsample"
        , values = c(16,2)
    )+
    xlab("Position on NS2A")+
    ylab("Freq.\nnonzero effect")+
    theme_classic()

#+ include = T
lapply(posFreq, function(x){
    x$nSampled %>%
    quantile(seq(0,1,by=0.1))
})

#' To see whether the 62 sites mapped to antigenic signals just because
#' they were more diverse (random association),
#' we permuted residues at each site across viruses to
#' dissolve the association with antigenic signals, if any,
#' while conserving the diversity at the sites.
aavers = c(
    "AA_NS2A_beyondE" = "NS2A sites"
    , "AA_NS2A_beyondE_permuted" = "NS2A sites (permuted)"
    , "AA_polyNonE_62" = "Sites on other proteins"
)
Colors$aavers = c(
    "NS2A sites" = "red"
    , "NS2A sites (permuted)" = "yellow"
    , "Sites on other proteins" = "grey"
)

rmse = lapply(names(aavers), function(aaver){
    file.path(
        "03-fits/nnls/thai_map/"
        , aaver
        , "adj_envelope/aasub/rmse_summary.csv"
    ) %>%
    read_csv(show_col_types = FALSE) %>%
    filter(dataset == "test_rmse") %>%
    arrange(lower, Median, upper) %>%
    head(100) %>%
    select(-protein_name) %>%
    mutate(siteSource = aavers[aaver]) %>%
    mutate(i = seq_along(siteSource))
})

legendTitle = "Source of 62AA"
gPermute = rmse[sapply(rmse, nrow) > 1] %>%
    do.call(what = rbind) %>%
    ggplot(aes(x = i, y = Median, color = siteSource, fill = siteSource))+
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA)+
    geom_line(size = 1.2)+
    geom_rect(data = rmse[sapply(rmse, nrow) == 1] %>% 
        do.call(what = rbind)
        , aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper)
        , alpha = 0.3
        , color = NA
    )+
    geom_hline(data = rmse[sapply(rmse, nrow) == 1] %>% 
        do.call(what = rbind)
        , aes(yintercept = Median, color = siteSource)
    )+
    scale_x_continuous(expand = c(0,0))+
    scale_color_manual( legendTitle
        , values = Colors$aavers
    )+
    scale_fill_manual( legendTitle
        , values = Colors$aavers
    )+
    theme_classic() %+%
    theme(
        panel.grid = element_blank()
        , legend.position = c(1,0.3)
        , legend.justification = c(1,0)
        , legend.background = element_rect(fill = NA, color = NA)
        , plot.margin = margin(5.5, 11, 3, 1)
    )+
    xlab("Subsample/permutation")+
    ylab("RMSE")

gExcessSite = ggarrange(
    plotSubsetRmse(60)
    , plotSubsetRmse(30)
    , gFreq
    , ncol = 1
    , nrow = 3
    , labels = letters[2:4]
)

#' Root-mean-squared-error (RMSE)
lapply(rmse, function(x){
    x %>%
        group_by(siteSource) %>%
        summarize(
            lower = mean(lower)
            , upper = mean(upper)
            , Median = mean(Median)
        )
}) %>%
do.call(what = rbind)


# Effect sizes of the 62 sites
gNS2A = ggarrange(
    domain %>%
        mutate(domain = if_else(grepl("^unnamed", domain), "", domain)) %>%
        ggplot(aes(y = 1))+
        geom_linerange(aes(
            xmin = posStart, xmax = posEnd
            , color = location
            , size = domain == ''
        ))+
        geom_text(aes(label = domain, x = posEnd)
            , hjust = -0.1, vjust = 0, nudge_y = 0.1, size = 2.5, angle = 90)+
        scale_size_manual(values = c(3,1), guide = F)+
        scale_color_manual(values = Colors$location, guide = F)+
        xlim(c(0,lenMax))+
        ylim(c(0.8,nrow(domain)+0.5))+
        theme(panel.grid = element_blank()
            , axis.title = element_blank()
            , axis.text = element_blank()
            , plot.margin = margin(0,0,0,1)
        )

    , d %>%
        ggplot(aes(x = pos, y = Median, color = location))+
        geom_linerange(aes(ymin = lower, ymax = upper), size = 0.4)+
        geom_point(size = 0.9)+
        scale_y_continuous("Effect size"
            , expand = c(0,0)
        )+
        scale_x_continuous("AA position (NS2A)"
            , limits = c(0, lenMax)
            , breaks = seq(5, lenMax, by = 50)
        )+
        scale_color_manual(values = Colors$location)+
        theme(
            panel.grid = element_blank()
            , axis.ticks = element_line(size = 0.8)
            , axis.ticks.length = unit(0.5, "lines")
            , panel.border = element_rect(color = "grey", fill = NA)
            , plot.margin = margin(0,0,0,1)
            , legend.position = 'bottom'
            , legend.justification = c(1, 1)
        )

    , ncol = 1
    , nrow = 2
    , heights = c(1,4)
    , align = "v"
    , common.legend = T
    , legend = "bottom"
)


# NS2A topology
posTransit = domain %>%
    unnest(cols = pos) %>%
    with(head(location, -1) != tail(location, -1)) %>%
    which


domain.topo = domain %>%
    group_by(location) %>%
    mutate(
        step.scale = if_else( location == "transmembrane"
            , max(posEnd - posStart + 1) / (posEnd - posStart + 1)
            , 1
        )
        , step.y = step.y * step.scale
    ) %>%
    ungroup %>%
    unnest(cols = pos) %>%
    mutate(
        x = cumsum(step.x) - step.x * (pos %in% posTransit)
        , y = cumsum(step.y) - step.y * (pos %in% posTransit)
    )
 
posFirst = domain.topo %>%
    filter(grepl("TMS", domain)) %>%
    group_by(domain) %>%
    filter(pos == min(pos))

domain.topo.membrane = domain.topo %>%
    mutate(xmin = min(x) - 5, xmax = max(x) + 5) %>%
    filter(location == "transmembrane") %>%
    group_by(xmin, xmax) %>%
    summarize(ymin = min(y), ymax = max(y)) %>%
    ungroup

gTopo = domain.topo %>%
    ggplot()+
    geom_rect(data = domain.topo.membrane
        , aes(xmin = xmin, xmax = xmax, ymin = ymin + 1.5, ymax = ymax - 0.2 )
        , fill = Colors$location[['transmembrane']]
    )+
    annotate("text", label = c("ER lumen","ER membrane","Cytosol")
        , x = domain.topo.membrane$xmin + 3
        , y = with(domain.topo.membrane, c(
              ymax + 1
            , (ymin+ymax)/2
            , ymin - 6
        ))
        , vjust = 0
        , hjust = 0
        , size = 2.5
        , color = c(
            Colors$location[['ER lumen']]
            , 'black'
            , Colors$location[['cytosol']]
        )
    )+
    geom_point(data = domain.topo %>%
            filter(grepl("TMS", domain))
        , aes(x = x, y = y)
        , shape = 15
        , size = 1.5
        , color = "#aaaaaa"
    )+
    geom_line(aes(x = x, y = y), size = 0.3)+
    geom_point(data = domain.topo %>% 
            filter(pos %in% unique(d$pos))
        , aes(x = x, y = y)
        , shape = 16
        , size = 1
    )+
    geom_text(data = posFirst %>% filter(step.y == 0)
        , aes(x = x, y = y + 8, label = domain)
        , hjust = 0
        , size = 2.5
    )+
    coord_fixed(ratio = 0.45)+
    theme_void()+
    theme(
        plot.margin = margin(0,0,-1,2, 'lines')
    )+
    ylim(with(domain.topo, range(y) + c(-5, 10)))




g = ggarrange(
    gExcess
    , gExcessSite
    , ggarrange( gPermute, gTopo, gNS2A
        , ncol = 1
        , nrow = 3
        , heights = c(2.5,1,3)
        , labels = letters[5:7]
    )
    , ncol = 3
    , nrow = 1
    , labels = c('a', '', '')
    , widths = c(1, 1.3, 1.3)
)
ggsave( g
    , filename = file.path(plotdir, "beyondE.jpg")
    , width = 11
    , height = 7
)



#'
#' # Distribution of sites on the protein
#'

# locations in the cell
domain %>%
    mutate(nSites = posEnd - posStart + 1) %>%
    group_by(location) %>%
    summarize(nSites = sum(nSites)) %>%
    left_join(    
        d %>%
            select(pos, location) %>%
            unique %>%
            group_by(location) %>%
            summarize(nSiteNonzero = n())
        , by = "location"
    ) %>%
    mutate(prop = nSiteNonzero/nSites)

x = domain %>%
    mutate(domain.label = if_else(
            grepl("^unnamed", domain)
            , ""
            , domain
        ) %>%
        paste0("\n(",posStart," to ",posEnd,")")
    )
domain.labels = x$domain.label
names(domain.labels) = x$domain


# Are they more associated with some region than expected?
posVariable = dat$aa[["nonstructural protein NS2A"]] %>%
    do.call(what = rbind) %>%
    apply(2, function(x){ length(unique(x)) > 1 }) %>%
    which

dProp = x %>%
    select(domain, domain.label, pos) %>%
    mutate(domain.label = gsub('\n',' ',domain.label)) %>%
    mutate(domain.label = gsub('^ +','',domain.label)) %>%
    mutate(domain.label = factor(domain.label, levels = rev(domain.label))) %>%
    unnest(cols = pos) %>%
    group_by(domain, domain.label) %>%
    summarize(
        countSiteAll = n()
        , countSite = sum(pos %in% posVariable)
    ) %>%    
    left_join( d %>%
        select(pos, domain) %>%
        unique %>%
        group_by(domain) %>%
        summarize(countEffect = n())
        , by = 'domain' 
    ) %>%
    mutate(
        countEffect = ifelse(is.na(countEffect), 0, countEffect)
    ) %>%
    ungroup %>%
    mutate( pvalue = 1 - pbinom(countEffect - 1, countSite, sum(countEffect)/sum(countSite)))

g = ggarrange(
    dProp %>%
        ggplot(aes(y = domain.label))+
        geom_col(mapping = aes(x = countSite), fill = 'black')+
        geom_col(mapping = aes(x = countEffect), fill = 'red')+
        geom_col(mapping = aes(x = countSiteAll), fill = NA, color = 'black')+
        scale_x_reverse('Number of sites', expand = c(0,0))+
        scale_y_discrete(position = 'right')+
        theme_classic()+
        theme(
            axis.text.y = element_blank()
            , axis.title.y = element_blank()
        )
    , dProp %>%
        ggplot(aes(y = domain.label, x = pvalue))+
        geom_col(fill = 'red')+
        geom_text(aes(label = round(pvalue,4) %>% format(digits = 3))
            , x = 0.02, hjust = 0, color = 'white', size = 2
        )+
        scale_x_continuous('Probability', expand = c(0,0))+
        theme_classic()+
        theme(
            axis.text.y = element_text(hjust = 0.35, vjust = 0.5, size = rel(0.8)) 
            , axis.title.y = element_blank()
            , axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        )

    , ncol = 2
    , nrow = 1
    , widths = c(2,2)
    , align = 'h'
    , labels = 'auto'
    , label.x = c(0,0.9)
)
ggsave( g
    , filename = file.path(plotdir, 'NS2A_domainAssoc.pdf')
    , width = 6
    , height = 5
)


#'
#' # Table of estimated effect sizes these 62 NS2A sites
#'

d %>%
    mutate(Position = pos
        , Substitution = Feature
        , Effect.median = round(Median, 2)
        , Effect.lower = round(lower, 2)
        , Effect.upper = round(upper, 2)
    ) %>%
    select(Position, Substitution, starts_with('Effect')) %>%
    arrange(Position, Effect.median) %>%
    DT::datatable()


#' 
#' # Improvement in % variance explained
#'


Rsq = new.env()
source("Scripts/plotEffects/computeRsq.R", local = Rsq)

R.beyond = '03-fits/nnls/thai_map/AA_NS2A_beyondE/adj_envelope/aasub/subset' %>%
    Rsq$calcRsq.allFolds()

R.E = '03-fits/nnls/thai_map/AA/adj_none/aasub/envelope_protein_E' %>%
    Rsq$calcRsq.allFolds()


# Improvement
(R.beyond[c("0%","100%"), ] - R.E[c("100%","0%"), ])[ , c("Rsq.all", "Rsq.inter", "Rsq.intra")]

