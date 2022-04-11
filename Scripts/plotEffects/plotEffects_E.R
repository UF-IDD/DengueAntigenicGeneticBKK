#' ---
#' title: Effects of E protein substitutions
#' author: "Angkana T. Huang"
#' output:
#'   html_document:
#'      toc: true
#' params:
#'     outdir: "04-plots"
#' ---


#+ include = FALSE
knitr::opts_chunk$set(include = TRUE, echo = FALSE, warning = FALSE)

# params = list(outdir = "04-plots")
library(tidyverse); theme_set(theme_minimal())
library(ggpubr)
attach(params)

source("Scripts/configs/plots.R")

# load database of known antigenically relevant substitutions
dat = new.env()
source("Scripts/antibodyDB_prep/db_load.R", local = dat, verbose = F)
aaver = "AA"
source("Scripts/sequence_load/aaver.R", local = dat, verbose = F)


# load short names for the proteins
x = read_csv("02-processedData/protein_shortname.csv", show_col_types = FALSE)
short = x$short
names(short) = x$protein_name

# labels for datasets
datasetLabels = c(
    "test_rmse" = "Out of sample"
    , "train_rmse" = "Training samples"
)
# colors for datasets
datasetColors = c(
    "test_rmse" = "red"
    , "train_rmse" = "black"
)

setPlotdir = function(x){
    plotdir <<- file.path(outdir, x)
    dir.create(plotdir)
}

# z-score of 95%CI
zCI = c(`2.5%` = -1.96, `50%` = 0, `97.5%` = 1.96)


#'
#' # AA diversity along the polyprotein
#  ....................................

setPlotdir('variations')

# define functions that given counts observed AA at the position
varaa = new.env()
# compute Shannon entropy
varaa$shannon = function(x){
    -sum(log2(x/sum(x)) * x/sum(x))   
}
# compute Simpson diversity index
varaa$simpson = function(x){
    1 - sum(x * (x-1))/(sum(x) * (sum(x)-1))
}
# compute Wu-kabat variability coefficient
varaa$wukabat = function(x){
    (sum(x) * length(x)) / max(x)
}

# select non-overlapping proteins that forms the polyprotein
protein_i = c(
    1   # anchored capsid
    , 5 # prM
    , 3 # envelope
    , 6 # NS1
    , 7 # NS2A
    , 8 # NS2B
    , 9 # NS3
    , 10 # NS4A
    , 12 # 2K
    , 11 # NS4B
    , 14 # NS5
)
x = dat$aa[protein_i]
# concatenate into polyprotein
dat$aapoly = lapply(names(x[[1]]), function(strain){
    lapply(x, function(aap){
        aap[[strain]]
    }) %>%
    unlist(use.names = F)
})
names(dat$aapoly) = names(x[[1]])

# compute variation measures
variations = do.call(rbind, dat$aapoly) %>%
    apply(2, function(aai){
        counts = table(aai)
        data.frame(
            shannon = varaa$shannon(counts)
            , simpson = varaa$simpson(counts)
            , wukabat = varaa$wukabat(counts)
        )
    }) %>%
    do.call(what = rbind) %>%
    as.data.frame %>%
    mutate(pos = seq_along(shannon))

measures = c(
    'shannon' = "Shannon entropy"
    , 'simpson' = "Simpson diversity index"
    , 'wukabat' = "Wu-kabat variability coefficient"
)
protein_regions = data.frame(
        protein_name = names(dat$aa)[protein_i]
        , protein_length = sapply(dat$aa[protein_i], function(aap){
            aap[[1]] %>% length
        })
    ) %>%
    mutate(
        protein_name = factor(protein_name, levels = protein_name)
        , protein_end = cumsum(protein_length)
        , protein_start = protein_end - protein_length + 1
    )
    
#+ include = T
lapply(names(measures), function(measure){
    variations %>%
        mutate(
            protein_name = factor(
                short[as.character(protein_regions$protein_name[ 
                    cut(pos, c(0,protein_regions$protein_end)) %>%
                    as.integer
                ])]
                , levels = short #[protein_regions$protein_name]
            )
        ) %>%
        ggplot(aes_(x = as.name(measure)))+
        geom_freqpoly(aes(y = ..density..), bins = 50)+
        facet_wrap( ~ protein_name, ncol = 3)+
        xlab(measures[measure])+
        ylab("Density")+
        theme_classic()+
        theme(
#             strip.text = element_text(size = 5)
            axis.text.x = element_text(angle = 90, hjust = 1)
        )
    ggsave(
        filename = file.path(plotdir, paste0("variation_",measure,".pdf"))
        , width = 4
        , height = 5
    )
}) %>%
invisible

#' Variations are low in most positions. Below shows the
#' distribution of variations along the polyprotein.
#+ include = T
ggarrange(
    protein_regions %>%
        ggplot(aes(y = protein_name))+
        geom_linerange(aes(xmin = protein_start, xmax = protein_end), size = 3)+
        geom_vline(aes(xintercept = protein_end), color = "red", size = 0.1)+
        scale_y_discrete(
            breaks = names(short)
            , labels = paste0(names(short), " (", short, ")")
            , position = 'right'
        )+
        theme(
            panel.grid = element_blank()
            , axis.title.y = element_blank()
            , axis.text.x = element_blank()
        )+
        facet_grid( "" ~ ., switch = "y")

    , variations %>%
        gather(measure, value, -pos) %>%
        ggplot(aes(x = pos, y = value))+
        geom_vline(data = protein_regions
            , aes(xintercept = protein_end), color = "red", size = 0.1
        )+
        geom_col(width = 1.1)+
        facet_grid( measure ~ ., scale = "free_y"
            , labeller = labeller(measure = gsub(" ", "\n", measures))
            , switch = "y"
        )+
        xlab("Amino acid position\non the polyprotein")+
        scale_y_continuous(position = 'right')+
        theme_classic()+
        theme(
            axis.title.y = element_blank()
        )
        
    , ncol = 1
    , nrow = 2
    , heights = c(0.8,2)
    , align = 'v'
)
#+ include = F
ggsave(filename = file.path(plotdir, "variation.pdf"))



#'
#' # E protein
#  ...........

setPlotdir('E')

fitdir = "03-fits/nnls/thai_map/AA/adj_none/aasub"
rmse = read_csv(file.path(fitdir, "rmse_summary.csv"), show_col_types = FALSE)
dE.all = read_csv(file.path(fitdir, "d_summary.csv"), show_col_types = FALSE) %>%
    filter(protein_name == "envelope_protein_E")
dE = dE.all %>%
    filter(lower > 0)

# length of E protein
lenE = dat$aa[['envelope protein E']][[1]] %>% length


#' Number of nonzero effect substitutions identified across a range of significance thresholds.
dE.all %>%
    mutate(prop_nonzero = n_nonzero / n_total) %>%   
    arrange(desc(prop_nonzero)) %>%
    mutate(i = seq_along(prop_nonzero)) %>%
    ggplot(aes(x = i, y = prop_nonzero, color = lower > 0))+
    geom_step(direction = 'h')+
    scale_color_manual(values = c('#444444', 'red'), guide = 'none')+
    xlab('Number of substitutions in E')+
    ylab('Proportion of estimations\nshowing nonzero effect')+
    theme_classic()
ggsave(
    filename = file.path(plotdir, 'significanceThresh_E.pdf')
    , width = 5
    , height = 3
)


#' 
#' `r nrow(dE)` nonzero effect substitutions, residing on
#' `r unique(dE$pos) %>% length` sites, showed nonzero effect
#' (95%IQR of effects of substitutions across the 100 estimations excluded zero).
#' Substitutions meeting our significance criteria 
#' were estimated as having nonzero effects in at least the following proportion of estimations they were involved in.
dE.all %>%
    mutate(prop_nonzero = n_nonzero / n_total) %>%   
    filter(lower > 0) %>%
    with(min(prop_nonzero))
#' Number of substitutions showing nonzero effect every time:
dE.all %>%
    mutate(prop_nonzero = n_nonzero / n_total) %>%   
    with(sum(prop_nonzero == 1))
#' Number of positions involved:
dE.all %>%
    filter(n_nonzero == n_total) %>%   
    select(pos) %>%
    unique %>%
    nrow
#' Positions involved:
dE.all %>%
    filter(n_nonzero == n_total) %>%   
    with(unique(pos) %>% sort)


#' Distribution of % Variance explained across the 100 estimations
#' are shown below.

Rsq = new.env()
source("Scripts/plotEffects/computeRsq.R", local = Rsq)

'03-fits/nnls/thai_map/AA/adj_none/aasub/envelope_protein_E' %>%
    Rsq$calcRsq.allFolds()

Pred = '03-fits/nnls/thai_map/AA/adj_none/aasub/envelope_protein_E' %>% Rsq$calcPred()


Rsq$DLong %>%
    cbind(Pred) %>%
    ggplot(aes(x = D))+
    geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`)
        , size = 0.5
        , alpha = 0.4
    )+
    geom_abline(slope = 1, color = 'red', linetype = 2)+
    coord_fixed(ratio = 1)+
    xlab("Observed distance")+
    ylab("Predicted distance")
ggsave(
    filename = file.path(plotdir, 'pred_E.pdf')
    , width = 5
    , height = 5
)



#'
#' ## Inferred effects vs known epitopes
#'

#' Number of monoclonal antibodies compiled in the DENVab database.
dat$ab %>%
    select(mAb) %>%
    unique %>%
    nrow
#' Of these, a subset were shown to have neutralization activity.
dat$ab %>%
    filter(`Assay Type` == "Neutralization") %>%
    select(mAb) %>%
    unique %>%
    nrow


# (known epitopes vs not) x diversity
dat$ab.E = dat$ab %>%
    filter(Protein == "E") %>%
    filter(`Epitope Type` == "Residue") %>%
    filter(`Assay Type` == "Neutralization") %>%
    unnest(cols = positions) %>%
    rename(pos = x1)

# Domains of E
domain.E = read_csv("00-RawData/domains/protein_envelope.csv"
    , comment = "#"
    , show_col_types = FALSE
)

variations = variations %>%
    filter(with(protein_regions %>% filter(protein_name == "envelope protein E")
        , between(pos, protein_start, protein_end)
    )) %>%
    # convert to position on E protein
    mutate(pos = pos - min(pos) + 1) %>%
    left_join( dat$ab.E %>%
        group_by(pos) %>%
        summarize(
            nReport = n()
            , nReportHuman = sum(Host == "Human")
        )
        , by = 'pos'
    ) %>%
    mutate(
        nReport = ifelse(is.na(nReport), 0, nReport)
        , nReportHuman = ifelse(is.na(nReportHuman), 0, nReportHuman)
    ) 

gPotentMab = dat$ab.E %>%
    filter(mAb %in% c(
        '1F4', '14C10', '2D22', '5J7', 'EDE1-2B2', 'EDE1-2C8'
    )) %>%
    ggplot(aes(x = pos, y = mAb, color = pos %in% dE$pos))+
    geom_point(shape = 1)+
    scale_color_manual(
        values = c("TRUE" = 'red', 'FALSE' = "#aaaaaa")
    )+
    guides(color = F)+
    xlim(c(0,lenE))+
    theme(
        panel.grid = element_blank()
        , axis.title = element_blank()
        , axis.text.x = element_blank()
        , axis.text.y = element_text(size = rel(0.65))
    )


#+ include = T
ggarrange(
    domain.E %>%
        mutate(i = as.integer(factor(domain, levels = unique(domain)))) %>%
        ggplot(aes(y = i))+
        geom_linerange(aes(xmin = posStart, xmax = posEnd), size = 3, color = "#cccccc")+
        geom_text(aes(label = domain, x = posStart), hjust = 0, nudge_y = 0.4, size = 3)+
        xlim(c(0,lenE))+
        ylim(c(0.8,nrow(domain.E)+0.5))+
        theme(panel.grid = element_blank()
            , axis.title = element_blank()
            , axis.text = element_blank()
            , plot.margin = margin(0,0,0,0)
        )

    # Where are the non-zero effect substitutions?
    , dE %>%
        ggplot(aes(x = pos, y = Median))+
        # known relevant sites, without variation in Thai viruses
        geom_segment(data = variations %>%
                    filter(shannon == 0, nReportHuman > 0)
            , aes(x = pos, xend = pos, y = max(dE$upper) * 0.90, yend = max(dE$upper)*1.02)
            , color = "#aaaaaa"
            , size = 0.2
            , linetype = 1
        )+
        # known relevant sites, with variations in Thai viruses
        geom_vline(data = variations %>%
                    filter(shannon > 0, nReportHuman > 0)
            , aes(xintercept = pos, alpha = nReportHuman)
            , color = "#aaaaaa"
            , size = 0.2
        )+
        scale_alpha(range = c(0.4,1), guide = F)+
        geom_linerange(aes(y = Median, ymin = lower, ymax = upper), size = 0.4)+
        geom_point(size = 0.7)+
        # known antigenic relevant substitutions found in dataset
        geom_point(data = dE %>%
                semi_join( dat$ab.E
                    , by = "pos"
                )
            , color = "red"
        )+
        scale_y_continuous("Effect size"
            , expand = c(0,0)
        )+
        scale_x_continuous("AA position (Envelope)"
            , limits = c(0, lenE)
            , breaks = seq(5, lenE, by = 50)
        )+
        theme(
            panel.grid = element_blank()
            , axis.ticks = element_line(size = 0.8)
            , axis.ticks.length = unit(0.5, "lines")
            , panel.border = element_rect(color = "grey", fill = NA)
            , plot.margin = margin(0,0,0,0)
        )
    , gPotentMab
    , ncol = 1
    , nrow = 3
    , align = 'v'
    , heights = c(1,6,1.2)
)
ggsave( filename = file.path(plotdir, "dE.jpg")
    , width = 6 * 1.1
    , height = 4 * 1.1
)

#' Contigency table of known epitopes vs observed variations in the dataset.
#+ include = T
tabE = dE %>%
    select(pos) %>%
    unique %>%
    mutate(nonzero = 1) %>%
    full_join( variations
        , by = 'pos'
    )
tabE %>%
    with(table(
        effect = ifelse(!is.na(nonzero), "non-zero", "zero")
        , mAb_database = ifelse(nReport > 0, "known epitopes", "not epitopes/unknown")
        , diversity = ifelse(shannon > 0, "has diversity", "no diversity")
    ))
tabE %>%
    with(table(
        effect = ifelse(!is.na(nonzero), "non-zero", "zero")
        , mAb_database = ifelse(nReportHuman > 0, "known human epitopes", "not human epitopes/unknown")
        , diversity = ifelse(shannon > 0, "has diversity", "no diversity")
    ))



# export to plot on E protein structure
tabE %>%
    mutate(
        variable = as.integer(shannon > 0)
        , epitope = as.integer(nReportHuman > 0)
        , nonzero = as.integer(!is.na(nonzero))
    ) %>%
    select(pos, nonzero, variable, epitope) %>%
    write_csv(file.path(plotdir, 'effect_diversity.csv'))



#' Odds ratio of being picked up by the model
fit = glm(nonzero ~ epitope
    , data = tabE %>%
        filter(shannon > 0) %>%
        mutate(
            nonzero = !is.na(nonzero)
            , epitope = nReport > 0
        )
    , family = "binomial"
) %>%
summary
#+ echo = T
exp(coef(fit)[2,1] + zCI * coef(fit)[2,2]) %>% round(2) # odds ratio, 95%CI


#' Odds ratio of being picked up by the model (human mAb)
fit = glm(nonzero ~ epitope
    , data = tabE %>%
        filter(shannon > 0) %>%
        mutate(
            nonzero = !is.na(nonzero)
            , epitope = nReportHuman > 0
        )
    , family = "binomial"
) %>%
summary
#+ echo = T
exp(coef(fit)[2,1] + zCI * coef(fit)[2,2]) %>% round(2) # odds ratio, 95%CI




#+ echo = F
# fetch the E protein sequences
dat$aaE = dat$aa[["envelope protein E"]]


#' ## Off diagonals of the contingency table
posEpitopes = tabE %>%
    filter(nReport > 0) %>%
    with(pos)
posEpitopesHuman = tabE %>%
    filter(nReportHuman > 0) %>%
    with(pos)
posNonzero = tabE %>%
    filter(!is.na(nonzero)) %>%
    with(pos)

#' ### Neighbors by linear distance
posVaryWithinWindow = function(x, wsize = 0, posExclude = integer(0)){
    tabE %>%
        filter(!(pos %in% posExclude)) %>%
        filter(shannon > 0) %>%
        with( intersect(pos
            , sapply(x, function(a) a + -wsize:wsize) %>%
                Reduce(f = union) %>%
                Filter(f = function(a) between(a,1,lenE))
        ))
}
#' Number of epitopes that were missed (although there was variability in the data)
#' that were within N sites from nonzero effect sites.
getNProximal = function(posEp){
    lapply(0:5, function(w){
        sitesMissed = setdiff(posVaryWithinWindow(posEp, 0),posNonzero)
        data.frame(
            windowSize = w
            , nMissed = length(sitesMissed)
            , nProximal = sum(sitesMissed %in% posVaryWithinWindow(posNonzero, w))
        ) %>%
        mutate( propProximal = nProximal / nMissed)
    }) %>%
    do.call(what = rbind) %>%
    rename(N = windowSize)    
}

getNProximal(posEpitopes)
getNProximal(posEpitopesHuman)


#' Probability of observing greater than or equal to this number of overlap
#' between nonzero effect sites and epitope regions (with the given the expansion window size).
getProbGe =  function(w, posEp, posExclude = integer(0)){
    data.frame(
        windowSize = w
        , nIntersect = intersect(posNonzero, posVaryWithinWindow(posEp, w)) %>%
            setdiff(posExclude) %>%
            length
        , nEpitopeSite = posVaryWithinWindow(posEp, w) %>%
            setdiff(posExclude) %>%
            length
    ) %>%
    mutate(
        propEpitopeSite = nEpitopeSite / 
            nrow(tabE %>% filter(shannon>0) %>% filter(!(pos %in% posExclude)))
        , propGe = 1 - pbinom(
            q = nIntersect - 1
            , posNonzero %>% setdiff(posExclude) %>% length
            , propEpitopeSite
        )
    )
}

# human + murine epitopes
probGe = lapply(0:10, getProbGe, posEp = posEpitopes) %>%
    do.call(what = rbind) 

# human epitopes
probGeHuman = lapply(0:10, getProbGe, posEp = posEpitopesHuman) %>%
    do.call(what = rbind) 

# nonhuman epitopes
probGeNonhuman = lapply(0:10, getProbGe, posEp = setdiff(posEpitopes, posEpitopesHuman)) %>%
    do.call(what = rbind) 


plotPropProximal = function(x){
    g = x %>%
        select(windowSize, propGe, propEpitopeSite) %>%
        gather(Type, prop, -windowSize) %>%
        ggplot(aes(x = windowSize, y = prop, fill = Type))+
        geom_col(position = 'dodge')+
        geom_text(data = x %>%
                mutate(Type = 'propGe')
            , aes(y = propGe, label = round(propGe,3) %>% sprintf(fmt = '%.3f'))
            , angle = 90
            , hjust = ifelse(max(x$propGe) > max(x$propEpitopeSite), 1, 0)
            , vjust = 1
            , nudge_x = 0.025
            , nudge_y = ifelse(max(x$propGe) > max(x$propEpitopeSite), -0.01, 0.01)
            , color = ifelse(max(x$propGe) > max(x$propEpitopeSite), 'white', 'red')
        )+
        scale_fill_manual(
            values = c('grey','red')
            , labels = c(
                'Proportion of variable sites within N sites of known epitopes'
                ,'Probability of the overlap being >= the observed by coincidence'
            )
        )+
        scale_x_continuous("N"
            , breaks = 0:10
        )+
        ylab("")+
        theme_classic()+
        theme(
            legend.title = element_blank()
            , legend.position = 'bottom'
            , legend.direction = 'vertical'
            # , axis.title.y = element_blank()
        )
    return(g)
}
ggarrange(
      plotPropProximal(probGe)
    , plotPropProximal(probGeHuman)
    , plotPropProximal(probGeNonhuman)
    , ncol = 1
    , nrow = 3
    , common.legend = T
    , legend = "bottom"
    , labels = letters[2:4]
    , font.label = list(size = 18)
)
ggsave( filename = file.path(plotdir, "propEpitopeOverlap.pdf")
    , width = 4 * 1.5
    , height = 4 * 1.5
)

probGe %>%
    rename(N = windowSize)
probGeHuman %>%
    rename(N = windowSize)







#' Number of non-epitope sites with nonzero effect
#' that were close to epitope sites that were not detected
tabE %>%
    filter(nReport == 0, shannon > 0, !is.na(nonzero)) %>%
    mutate(
        deltaEpitopeZero = sapply(pos, function(pe){
            (pe - setdiff(posEpitopes, posNonzero)) %>% abs %>% min
        })
        , deltaEpitope = sapply(pos, function(pe){
            (pe - posEpitopes) %>% abs %>% min
        })
    ) %>%
    with( lapply(1:max(deltaEpitopeZero), function(x){
        data.frame(
            delta = x
            , numNonzero.fromEpitopeZero = sum(deltaEpitopeZero <= x)
            , numNonzero.fromEpitope = sum(deltaEpitope <= x)
        )
    })) %>%
    do.call(what = rbind) %>%
    head(3)


#' ### Neighbors by structural distance

# import structural distance
dat$sdist = read_csv("02-processedData/structure/ResidueDist.csv", col_types = 'iid')

dat$sdist %>%
    ggplot(aes(x = Dist))+
    geom_histogram(binwidth = 1, fill = 'black')+
    xlab('Angstrom')+
    ylab('Position pairs')


# function to get positions within a structural distance window from epitope positions that has variation
posVaryWithinWindow = function(x, wsize = 0, posExclude = integer(0)){
    tabE %>%
        filter(!(pos %in% posExclude)) %>%
        filter(shannon > 0) %>%
        with( intersect(pos
            , dat$sdist %>%
                filter(Dist <= wsize) %>%
                filter(i %in% x | j %in% x) %>%
                with(union(i,j)) %>%
                union(x)
        ))
}

window.vector = seq(1,8, by = 0.5)

# human + murine epitopes
probGe = lapply(window.vector, getProbGe, posEp = posEpitopes) %>%
    do.call(what = rbind) 

# human epitopes
probGeHuman = lapply(window.vector, getProbGe, posEp = posEpitopesHuman) %>%
    do.call(what = rbind) 

# nonhuman epitopes
probGeNonhuman = lapply(0:10, getProbGe, posEp = setdiff(posEpitopes, posEpitopesHuman)) %>%
    do.call(what = rbind) 




#########

plotPropProximal = function(x){
    g = x %>%
        select(windowSize, propGe, propEpitopeSite) %>%
        gather(Type, prop, -windowSize) %>%
        ggplot(aes(x = windowSize, y = prop, fill = Type))+
        geom_col(position = 'dodge')+
        geom_text(data = x %>%
                mutate(Type = 'propGe')
            , aes(y = propGe, label = round(propGe,3) %>% sprintf(fmt = '%.3f'))
            , angle = 90
            , hjust = ifelse(max(x$propGe) > max(x$propEpitopeSite), 1, 0)
            , vjust = 1
            , nudge_x = 0.025
            , nudge_y = ifelse(max(x$propGe) > max(x$propEpitopeSite), -0.01, 0.01)
            , color = ifelse(max(x$propGe) > max(x$propEpitopeSite), 'white', 'red')
        )+
        scale_fill_manual(
            values = c('grey','red')
            , labels = c(
                'Proportion of variable sites within N sites of known epitopes'
                ,'Probability of the overlap being >= the observed by coincidence'
            )
        )+
        scale_x_continuous("N"
            , breaks = 0:10
        )+
        ylab("")+
        theme_classic()+
        theme(
            legend.title = element_blank()
            , legend.position = 'bottom'
            , legend.direction = 'vertical'
            # , axis.title.y = element_blank()
        )
    return(g)
}


#######





plotPropProximal = function(x){
    g = x %>%
        select(windowSize, propGe, propEpitopeSite) %>%
        gather(Type, prop, -windowSize) %>%
        ggplot(aes(x = windowSize, y = prop, fill = Type))+
        geom_col(position = 'dodge')+
        geom_text(data = x %>%
                mutate(Type = 'propGe')
            , aes(y = propGe, label = round(propGe,3) %>% sprintf(fmt = '%.3f'))
            , angle = 90
            , hjust = ifelse(max(x$propGe) > max(x$propEpitopeSite), 1, 0)
            , vjust = 1
            , nudge_x = 0.025
            , nudge_y = ifelse(max(x$propGe) > max(x$propEpitopeSite), -0.01, 0.01)
            , color = ifelse(max(x$propGe) > max(x$propEpitopeSite), 'white', 'red')
        )+
        scale_fill_manual(
            values = c('grey','red')
            , labels = c(
                'Proportion of variable sites within X angstroms of known epitopes'
                ,'Probability of the overlap being >= the observed by coincidence'
            )
        )+
        scale_x_continuous("X"
            , breaks = 0:max(window.vector)
        )+
        ylab("")+
        theme_classic()+
        theme(
            legend.title = element_blank()
            , legend.position = 'bottom'
            , legend.direction = 'vertical'
        )
    return(g)
}            

ggarrange(
      plotPropProximal(probGe)
    , plotPropProximal(probGeHuman)
    , plotPropProximal(probGeNonhuman)
    , ncol = 1
    , nrow = 3
    , common.legend = T
    , legend = "bottom"
    , labels = letters[5:7]
    , font.label = list(size = 18)

)
ggsave( filename = file.path(plotdir, "propEpitopeOverlap_angstrom.pdf")
    , width = 4 * 1.5
    , height = 6 * 1.5
)

probGe %>%
    rename(N = windowSize)
probGeHuman %>%
    rename(N = windowSize)




#' Where are the nonzero effect sites that were in the anchor/stem domain?
tabE %>%
    filter(nReport == 0, shannon > 0, !is.na(nonzero)) %>%
    filter(with( domain.E %>% filter(domain == "Stem/anchor")
        , between(pos, posStart, posEnd)
    )) %>%
    with(pos)



#'
#' ## Estimated effect sizes
#'

dE.out = dE %>%
    mutate(
        Position = pos
        , Substitution = mapply(gsub,'[0-9]+', pos, Feature) %>% toupper
        , Effect.median = Median
        , Effect.lower = lower
        , Effect.upper = upper
    ) %>%
    arrange(Position, Effect.median) %>%
    select(Position, Substitution, starts_with('Effect')) %>%
    mutate_at(vars(matches('^Effect')), round, digits = 3)
    
dE.out %>%
    DT::datatable()
    
dE.out %>%
    write_csv(file.path(plotdir, 'Subs_nonzero_E.csv'))