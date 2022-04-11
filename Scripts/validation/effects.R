library(Matrix)
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

outdir = "05-validateEstimates"
dir.create(outdir, recursive = T)

# import nonzero effect substitutions identified
source("Scripts/validation/load_effectsize.R")

# import available titer measurements
titers = read.table("00-RawData/titers/MAY15-2020-combined-set1set2set3set4set5set6set7-sciencem3-exp-conditions-reasonable-adjust.txt")
# keep just the 20 main ones
titers = titers[ , 1:20]
# recode titers <x
processTiters = function(x){
    mult = ifelse(grepl('^<', x), 0.5, 1) # multiplier for titer
    mult = ifelse(grepl('^>', x),   2, mult)
    mult * as.integer(str_extract(x, '[0-9]+'))
}
titers = titers %>%
    mutate_all(processTiters)

# import antigenic distances
dMap = read_csv("02-processedData/antigenicDist/thai_map.csv", col_types = "ccn") %>%
    rename(v1 = serum, v2 = virus, mapDist = D)

# import mAb database
dat = new.env()
source("Scripts/antibodyDB_prep/db_load.R", local = dat, verbose = F)
dat$ab.E = 
    dat$ab %>%
    filter(Protein == "E") %>%
    filter(`Epitope Type` == "Residue") %>%
    filter(`Assay Type` == "Neutralization") %>%
    unnest(cols = positions) %>%
    rename(pos = x1)
dat$ab.E = split(dat$ab.E, dat$ab.E$Host)



# load short names for the proteins
x = read_csv.silent("02-processedData/protein_shortname.csv")
short = x$short
names(short) = x$protein_name
Long = names(short); names(Long) = short

# import sequences
aaver = "AA"
source("Scripts/sequence_load/aaver.R")
names(aa) = short[ names(aa) ]
# reduce to just the proteins with nonzero effects
aa = aa[names(posNonzero)]
# reduce to just Thai viruses that we have titers for and have been mapped
aa = lapply(aa, function(x){
    x[ intersect(names(x), union(dMap$v1, dMap$v2)) ]
})

# Domains of E
domain.E = read_csv("00-RawData/domains/protein_envelope.csv"
    , comment = "#"
    , show_col_types = FALSE
)





    #   Relationship between fold-change and mapped distance
    #   .............

getFoldChange = function(x, col.v1 = "v1", col.v2 = "v2", outcol.suffix = ""){
    fold = titers[x[[col.v1]], ] / titers[x[[col.v2]], ]
    i = fold<1
    i = !is.na(i) & i
    fold[i] = 1 / (fold[i])
    
    x[[paste0('foldAvg', outcol.suffix)]] = rowMeans(fold-1)
    x
}

dMap %>% 
getFoldChange %>%
    mutate(pairType = ifelse(substring(v1,1,5)==substring(v2,1,5), "Homotypic", "Heterotypic")) %>%
    ggplot(aes(x = mapDist, y = foldAvg, color = pairType))+
    geom_abline(slope = 1, intercept = 0)+
    geom_point(shape = 1, stroke = 0.1, size = 0.3, alpha = 0.6)+
    scale_x_continuous('Antigenic distance', breaks = 0:10)+
    scale_y_continuous('Average fold change in titers', breaks = seq(0,50, by =5))+
    scale_color_manual("Virus pair"
        , values = c('red','#aaaaaa')
    )+
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )+
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1, stroke = 1)))




    #   Identify substitutions between virus pairs
    #   .............

# non-zero effect substitutions in NS2A from the model
subs = mapply( function(eff, Seqs){
        lapply(names(Seqs), function(v1){
        lapply(names(Seqs), function(v2){
            # skip comparing the same virus
            if(v1 == v2){ return(NULL) }
            # skip comparing heterotypic pairs
            if(substring(v1,1,5) != substring(v1,1,5)){ return(NULL) }            
            
            subs = paste0(Seqs[[v1]], seq_along(Seqs[[v1]]), Seqs[[v2]])[Seqs[[v1]] != Seqs[[v2]]]
            subs = intersect(subs, names(eff))
            out = data.frame(v1, v2, D = sum(eff[subs]))
            out$sub = list(subs)
            out$nsub = length(subs)
            out
        }) %>% do.call(what = rbind)
        }) %>% do.call(what = rbind)
    }
    , eff = effectMid['NS2A']
    , Seqs = aa['NS2A']
    , SIMPLIFY = F
)

# any substitution in E
subs[['E.any']] = 
    lapply(names(aa[['E']]), function(v1){
    lapply(names(aa[['E']]), function(v2){
        # skip comparing the same virus
        if(v1 == v2){ return(NULL) }
        # skip comparing heterotypic pairs
        if(substring(v1,1,5) != substring(v1,1,5)){ return(NULL) }            

        i = which(aa[['E']][[v1]] != aa[['E']][[v2]])
        out = data.frame(v1, v2, nsub = length(i))
        out$sub = list(paste0(aa[['E']][[v1]][i], i, aa[['E']][[v2]][i]))
        out
        
    }) %>% do.call(what = rbind)
    }) %>% do.call(what = rbind)

# hash the substitution vectors to ease joining
hashSubs = function(sub){
    sapply(sub, function(sx){
        sort(sx) %>%
        paste0(collapse = ';') %>% 
        digest::digest(algo = "md5")
    })
}
subs = lapply(subs, function(x){
    x %>% mutate(sub.hash = hashSubs(sub))
})



    #   Group substitutions into mutually exclusive groups by positions
    #   ..................................................

posSet.E = list(
    # 1) Known epitopes: human-derived mAb (hmAb)
    'hmAb' = dat$ab.E[['Human']]$pos %>% sort %>% unique
    # 2) Known epitopes: murine-derived mAb (mmAb)
    , 'mmAb' = setdiff(dat$ab.E[['Murine']]$pos, dat$ab.E[['Human']]$pos) %>% sort %>% unique
    # 3) Not in known epitopes: EDI/II/III
    , 'E.nonstem' = str_extract(names(effectMid[['E']]), '[0-9]+') %>% 
        as.integer %>% 
        unique %>%
        setdiff(dat$ab.E[['Murine']]$pos) %>%
        setdiff(dat$ab.E[['Human']]$pos) %>%
        Filter(f = function(x) x < 394) %>%
        sort
    # 4) Not in known epitopes: Stem/Anchor
    , 'E.stem' = str_extract(names(effectMid[['E']]), '[0-9]+') %>% 
        as.integer %>% 
        unique %>%
        setdiff(dat$ab.E[['Murine']]$pos) %>%
        setdiff(dat$ab.E[['Human']]$pos) %>%
        Filter(f = function(x) x >= 394) %>%
        sort
)

# check that they are mutually exclusive
if(length(Reduce(posSet.E, f = union)) != length(Reduce(posSet.E, f = c))){
    stop('They are not mutually exclusive!')
}



    #   Functions to validate effects of substitutions
    #   .............................

getDist = new.env()
getDist$foldChangeAvg = function(triplet, Sub){
    triplet %>%
        getFoldChange() %>%
        getFoldChange(col.v2 = "v2.control", outcol.suffix = ".control") %>%
        rename_at(vars(matches('^foldAvg')), function(x) gsub('^foldAvg', 'd', x))
}
getDist$mapDist = function(triplet, Sub){
    triplet %>%
        left_join(dMap, by = c('v1','v2')) %>%
        left_join(dMap, by = c('v1','v2.control' = 'v2'), suffix = c('','.control')) %>%
        rename_at(vars(matches('^mapDist')), function(x) gsub('^mapDist', 'd', x))
}

# given a substitution, return how many virus pairs have these in the specified dataset
SubPresent = function(sx, subs.list){
    sapply(subs.list, function(x) sx %in% x)
}
tabSubN = function(sx, subs.full = subs[['E.any']]){
    Present = SubPresent(sx, subs.list = subs.full$sub)
    tibble(
        sub = sx
        , Nvirus.full = subs.full$v2[Present] %>% unique %>% length
        , Npair.full = Present %>% sum
    )
}

# given a substitution, extract virus triplets and perform analysis
compareTriplet = function(sx
    , sx.eff
    , subs.full
    , metric = 'mapDist'
    , allType = TRUE
    , plotSepV2 = FALSE
){
    dvar.label = c(
        "mapDist" = "antigenic distance"
        , "foldChangeAvg" = "mean fold change in titers"
    )[metric]

    triplet = subs.full %>%
        filter(SubPresent(sx, sub)) %>%
        mutate(
            sub = lapply(sub, function(sl){ setdiff(sl, sx) })
            , sub.hash.control = hashSubs(sub)
        ) %>%
        select(-sub) %>%
        inner_join(subs.full %>%
                select(v1,v2,sub.hash,sub.hash.other)
            , by = c('v1', 'sub.hash.control' = 'sub.hash', 'sub.hash.other')
            , suffix = c('', '.control')
        ) %>%
        getDist[[metric]](Sub = sx) %>%
        filter(!is.na(d), !is.na(d.control))
    
    out = tibble(
        Nvirus = unique(triplet$v2) %>% length
        , Npair = triplet %>% select(v1,v2) %>% unique %>% nrow
    )
    # perform analysis if there are ample datapoints
    if(out$Nvirus > 2 & out$Npair > 30){
        if(plotSepV2){
            g = lapply(c('v2','v2.control'), function(vgroup){
                triplet$vgroup = triplet[[vgroup]]
                triplet %>%
                    mutate_at(vars(matches('^v[12]')), function(x) substring(x,1,5)) %>%
                    group_by_at(vars(matches('^v[12]'), vgroup)) %>%
                    mutate(delta = d - d.control) %>%
                    summarize(
                        Ntriplet = n()
                        , pvalue = mean(delta <= 0)
                        , pvalue.txt = ifelse(pvalue<0.001, "p<0.001", sprintf("p=%.3f", pvalue))
                        , dLower = quantile(delta, 0.025)
                        , dMid = quantile(delta, 0.5)
                        , dUpper = quantile(delta, 0.975)
                    ) %>%
                    ggplot(aes(y = vgroup, color = pvalue <= 0.05))+
                    geom_vline(xintercept = 0, linetype = 3, size = 0.5)+
                    geom_linerange(aes(xmin = dLower, xmax = dUpper))+
                    geom_point(aes(x = dMid))+
                    geom_text(aes(x = dMid, label = pvalue.txt), nudge_y = rel(0.25), size = rel(3))+
                    facet_wrap( ~ v1 + v2, scale = 'free_y'
                        , ncol = 1
                        , labeller = label_wrap_gen(multi_line=FALSE)
                    )+
                    scale_color_manual(values = c('FALSE'='grey', 'TRUE'='red'), guide = 'none')+
                    ylab(vgroup)+
                    xlab(paste0('Difference in\n', dvar.label))
            })
            g = ggarrange(plotlist = g
                    , ncol = 2
                    , nrow = 1
                ) %>%
                annotate_figure(top = paste0(
                    toupper(sx)
                    , ' (d=', round(sx.eff,2), ')'
                ))
            plot(g)            
        }
        
        # an overal p-value for each serotype pair
        out = triplet %>%
            mutate(v1name = v1) %>%
            mutate(v2name = v2) %>%
            mutate(v2controlname = v2.control) %>%
            mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
            # whether to combine all serotypes
            mutate(
                v1 = ifelse(allType, 'All', v1)
                , v2 = ifelse(allType, 'All', v2)
            ) %>%
            group_by(v1, v2) %>%
            mutate(delta = d - d.control) %>%
            summarize(
                N.v1 = unique(v1name) %>% length
                , N.v2 = unique(v2name) %>% length
                , N.v2.control = unique(v2controlname) %>% length
                , pvalue = mean(delta <= 0)
                , dLower = quantile(delta, 0.025)
                , dMid = quantile(delta, 0.5)
                , dUpper = quantile(delta, 0.975)
            )
    } else{
        out = tibble(
                # v1 = NA, v2 = NA, Ntriplet = NA, pvalue = NA, dLower = NA, dMid = NA, dUpper = NA
                v1 = NA, v2 = NA
                , N.v1 = out$Npair/out$Nvirus
                , N.v2 = out$Nvirus
                , N.v2.control = NA
                , pvalue = NA, dLower = NA, dMid = NA, dUpper = NA
            )
    }
    out
}


# function to plot delta.d and pvalue onto phylogeny
library(treeio)
library(ggtree)
plotTree = function(v2data, subname, protein){
    tr = file.path(
        "00-RawData/phylogeny"
        , list.files(
            "00-RawData/phylogeny"
            , pattern = paste0('^',gsub('DENV', 'd', v2data$v2[1]), '.+\\.treefile$')
        )
    ) %>%
    read.tree
    tr = rename_taxa(tr
        , as_tibble(tr) %>%
            mutate(label.new = str_extract(label, "[0-9]+/[0-9]{2}$"))
        , label
        , label.new
    )

    # MRCA of all v2
    Mrca = groupOTU(tr, v2data$tip.label) %>%
        as_tibble %>%
        filter(group == 1) %>%
        with(setdiff(parent, node))
    try({
        tr = tree_subset(tr, Mrca)    
    })
    
    tr %>%
        ggtree(aes(color = dMid), size = 0.3, layout="slanted") %<+% v2data +
        geom_tippoint(aes(color = dMid, alpha = is.na(pvalue), shape = pvalue), stroke = 1, size = 1.5)+
        scale_shape_manual(values = c('[0,0.05]' = 16, '(0.05,1]' = 2))+
        scale_alpha_manual(values = c('TRUE' = 0, 'FALSE' = 1))+
        scale_color_gradient2('Difference in antigenic distance'
            , low = "purple"
            , mid = "#cccccc"
            , high = "red"
            , na.value = "#cccccc"
        )+
        scale_x_continuous(expand = c(0.1,0))+
        # geom_treescale(fontsize=3, linesize=1.5, offset=2, x = 0.01, -0.5)+
        guides(size = 'none', alpha = 'none', shape = 'none')+
        theme(
            plot.margin = unit(c(0.5,0,0.5,0), "lines")
        )+
        annotate("text", x = 0, y = 0
            # , vjust = 1.2, hjust = 0
            , vjust = 0, hjust = 0, angle = 90
            , size = rel(2.5)
            , label = paste0(protein, ': ', toupper(subname), '\n(',v2data$v1[1],',', v2data$v2[1], ')')
        )
}

plotTreeFromSub = function(sx, protein, metric = 'mapDist'){
    subs.full = left_join(subs[[c('E' = 'E.any', 'NS2A' = 'NS2A')[protein]]]
        , subs[[c('E' = 'NS2A', 'NS2A' = 'E.any')[protein]]] %>% select(v1,v2,sub.hash) %>% rename(sub.hash.other = sub.hash)
        , by = c('v1','v2')
    )

    triplet = subs.full %>%
        filter(SubPresent(sx, sub))
    out = tibble(
        Nvirus = unique(triplet$v2) %>% length
        , Npair = triplet %>% select(v1,v2) %>% unique %>% nrow
    )
    # perform analysis if there are ample datapoints
    if(out$Nvirus > 2 & out$Npair > 30){
        triplet = triplet %>%
                mutate(
                    sub = lapply(sub, function(sl){ setdiff(sl, sx) })
                    , sub.hash.control = hashSubs(sub)
                ) %>%
                select(-sub) %>%
                inner_join(subs.full %>%
                        select(v1,v2,sub.hash,sub.hash.other)
                    , by = c('v1', 'sub.hash.control' = 'sub.hash', 'sub.hash.other')
                    , suffix = c('', '.control')
                ) %>%
                getDist[[metric]](Sub = sx) %>%
                filter(!is.na(d), !is.na(d.control)) %>%
                mutate(delta = d - d.control)

        # compute pvalue for each v2
        v2data = triplet %>%
            mutate(v2name = v2) %>%
            mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
            group_by(v1,v2,v2name) %>%
            mutate(delta = d - d.control) %>%
            summarize(
                , pvalue = mean(delta <= 0)
                , dMid = quantile(delta, 0.5)
                , .groups = 'drop'
            ) %>%
            mutate(pvalue = cut(pvalue, breaks = c(0, 0.05, 1), include.lowest = T)) %>%
            mutate(tip.label = str_extract(v2name, "[0-9]+/[0-9]{2}$"), .before = 1)

        # plot results onto phylogeny
        g = split(v2data, with(v2data, paste(v1,v2))) %>%
            lapply(plotTree, subname = sx, protein = protein)
        return(g)
    } else {
        return(NULL)
    }
}


    #   Compute p-values
    #   ................

# Positions in E
# ...................
subs.E = Reduce(subs[['E.any']]$sub, f = union)

lapply(c('allType', 'pairType'), function(ver){
    lapply(names(posSet.E), function(prefix, exportdir = outdir){
        # all substitutions happening in these positions
        Subs = subs.E[sapply(str_extract(subs.E, '[0-9]+') %>% as.integer, function(x) x %in% posSet.E[[prefix]])]
        # effect size of these substitutions
        eff = sapply(Subs, function(sub) ifelse(sub %in% names(effectMid[['E']]), effectMid[['E']][[sub]], 0))
        
        dir.create(exportdir, recursive = T)
        mapply(function(sub.g, eff.g, effgroup){
                lapply(ls(getDist)[2], function(metric){
                    outstem = file.path(exportdir, paste(prefix, effgroup, metric, sep = "_"))
                    pdf(
                        paste0(outstem, '.pdf')
                        , height = 8
                        , width = 12
                    )
                    out = lapply(sub.g, function(sub){
                            cbind(
                                tabSubN(sub)
                                , compareTriplet(sub
                                    , eff.g
                                    , subs.full = left_join(subs[['E.any']]
                                        , subs[['NS2A']] %>% select(v1,v2,sub.hash) %>% rename(sub.hash.other = sub.hash)
                                        , by = c('v1','v2')
                                    )
                                    , metric = metric
                                    , allType = ver == 'allType'
                                )
                            )
                        }) %>%
                        do.call(what = rbind)
                    dev.off()
                    out %>%
                        # arrange(desc(Nvirus), desc(Npair), desc(Nvirus.full), desc(Npair.full)) %>%
                        write_csv(paste0(outstem, '.csv'))
                })
            }
            , sub.g = split(Subs, eff==0)
            , eff.g = split(eff , eff==0)
            , effgroup = c('nonzero', 'zero')
            , SIMPLIFY = F
        )
    }, exportdir = file.path(outdir, ver))
}) %>%
invisible



# Positions in NS2A
# ...................

lapply(c('allType', 'pairType'), function(ver){
    # ver = 'allType'
    exportdir = file.path(outdir, ver)
    dir.create(exportdir, recursive = T)
    lapply(ls(getDist)[2], function(metric){
        outstem = file.path(exportdir, paste('NS2A', "nonzero", metric, sep = "_"))
        pdf(
            paste0(outstem, '.pdf')
            , height = 8
            , width = 12
        )
        out = lapply(names(effectMid[['NS2A']]), function(sub){
                cbind(
                    tabSubN(sub, subs.full = subs[['NS2A']])
                    , compareTriplet(sub
                        , effectMid[['NS2A']][[sub]]
                        , subs.full = left_join(subs[['NS2A']]
                            , subs[['E.any']] %>% select(v1,v2,sub.hash) %>% rename(sub.hash.other = sub.hash)
                            , by = c('v1','v2')
                        )
                        , metric = metric
                        , allType = ver == 'allType'
                    )
                )
            }) %>%
            do.call(what = rbind)
        dev.off()
        out %>%
            write_csv(paste0(outstem, '.csv'))
    })
}) %>%
invisible



    #   Factors of sample size for triplet analysis
    #   ....................................


GroupLabels = c(
    'hmAb' =  'hmAb'
    , 'mmAb' = 'mmAb only'
    , 'E.nonstem' = 'EDI/II/III'
    , 'E.stem' = 'E Stem/Anchor'
    , 'NS2A' = 'NS2A'
)
GroupLabels.inv = names(GroupLabels)
names(GroupLabels.inv) = GroupLabels


P = list.files(file.path(outdir, 'allType'), pattern = '_nonzero_mapDist\\.csv$') %>%
    lapply(function(f){
        meta = strsplit(gsub('\\.csv$','',f), '_') %>% do.call(what = rbind) %>% as.data.frame
        colnames(meta) = c('Group', 'Effect', 'Metric')
        p = read_csv(file.path(outdir, 'allType', f), col_types = "ciicciiinnnn")
        meta %>%
            mutate(
                p = list(p)
                , Group = factor(Group, levels = names(GroupLabels), labels = GroupLabels)
            )
    }) %>%
    do.call(what = rbind) %>%
    unnest(cols = p)

# Substitutions not listed here are ones that did not find triplets at all
split(P, ifelse(P$Group=='NS2A', 'NS2A', 'E')) %>%
    lapply(dim)

# Nvirus.full - N.v2 = number of viruses that failed to find a matching control
# : listed because some control was found but not enough
P %>%
    filter(Nvirus.full > 2) %>%
    filter(N.v2 <= 2) %>%
    dim

# Substitutions limited by number of virus harboring that substitution
P %>%
    filter(Nvirus.full <= 2) %>%
    dim

# Reasons limiting assessment of large effect substitutions
lapply(effect, function(x){
     x %>% 
        filter(lower > 0, Median > 0.5) %>%
        left_join(P %>% 
                mutate(protein_name = ifelse(Group=='NS2A', 'NS2A', 'E'))
            , by = c('protein_name', 'sub')
        ) %>%
        summarize(
            total.sub.large = n()
            , num.sub.no.control = sum(
                (Nvirus.full > 2 & N.v2 <= 2) |
                (Nvirus.full > 2 & (N.v1*N.v2) <= 30)
            )
        )
})




    #   Trends in p-values and number of data points
    #   ..................


# All serotypes combined
# delta.d vs estimated effect size
P = list.files(file.path(outdir, 'allType'), pattern = '_nonzero_mapDist\\.csv$') %>%
    lapply(function(f){
        meta = strsplit(gsub('\\.csv$','',f), '_') %>% do.call(what = rbind) %>% as.data.frame
        colnames(meta) = c('Group', 'Effect', 'Metric')
        p = read_csv(file.path(outdir, 'allType', f), col_types = "ciicciiinnnn")
        meta %>%
            mutate(
                Nsub = p$sub %>% unique %>% length
                , Nsub.eval = p %>% filter(!is.na(pvalue)) %>% with(sub %>% unique %>% length)
                , p = list(p %>% filter(!is.na(pvalue)))
                , Group = factor(Group, levels = names(GroupLabels), labels = GroupLabels)
            )
    }) %>%
    do.call(what = rbind) %>%
    unnest(cols = p) %>%
    mutate(
        sub = factor(sub)
        , protein_name = ifelse(Group == 'NS2A', 'NS2A', 'E')
    ) %>%
    left_join( effect %>%
        lapply(function(x){
            x %>%
            select(protein_name, sub, lower, Median, upper) %>%
            rename(
                effLower = lower
                , effMid = Median
                , effUpper = upper
            )
        }) %>%
        do.call(what = rbind)
        , by = c('protein_name', 'sub')
    ) 


# correspondence between estimated and observable effects
P %>% with(cor(effMid, dMid))  # correlation overall
P %>%
    ggplot(aes(x = dMid, y = effMid))+
    geom_abline(slope = 1, linetype = 3)+
    geom_point(size = 0.8)+
    geom_linerange(aes(xmin = dLower, xmax = dUpper))+
    geom_linerange(aes(ymin = effLower, ymax = effUpper))+
    facet_wrap( ~ Group, nrow = 1)+
    xlab('Observed difference in antigenic distance')+
    ylab('Estimated effect size of substitution')
ggsave(filename = file.path(outdir, 'effectRelationship.pdf')
    , width = 8
    , height = 4
)


# delta.d
# Overall vs Separated by serotype pairs
P.pair = list.files(file.path(outdir, 'pairType'), pattern = '\\.csv$') %>%
    lapply(function(f){
        meta = strsplit(gsub('\\.csv$','',f), '_') %>% do.call(what = rbind) %>% as.data.frame
        colnames(meta) = c('Group', 'Effect', 'Metric')
        p = read_csv(file.path(outdir, 'pairType', f), col_types = "ciicciiinnnn")
        meta %>%
            mutate(
                Nsub = p$sub %>% unique %>% length
                , Nsub.eval = p %>% filter(!is.na(pvalue)) %>% with(sub %>% unique %>% length)
                , p = list(p %>% filter(!is.na(pvalue)))
            )
    }) %>%
    do.call(what = rbind) 


# more functions for this plot
computePvalue = function(triplet){
    triplet %>%
        summarize(pvalue = mean(delta <= 0)) %>%
        mutate(pvalue = ifelse(pvalue<0.001, "p<0.001", sprintf("p=%.3f", pvalue)))
}
getBreaks = function(x, digits = 1){
    seq(round(min(x-10^-digits), digits), round(max(x+10^-digits), digits), by =0.1)
}
plotPvalueDistribution = function(sx, protein, metric = 'mapDist'){
    subs.full = left_join(subs[[c('E' = 'E.any', 'NS2A' = 'NS2A')[protein]]]
        , subs[[c('E' = 'NS2A', 'NS2A' = 'E.any')[protein]]] %>% select(v1,v2,sub.hash) %>% rename(sub.hash.other = sub.hash)
        , by = c('v1','v2')
    )

    triplet = subs.full %>%
        filter(SubPresent(sx, sub)) %>%
            mutate(
                sub = lapply(sub, function(sl){ setdiff(sl, sx) })
                , sub.hash.control = hashSubs(sub)
            ) %>%
            select(-sub) %>%
            inner_join(subs.full %>%
                    select(v1,v2,sub.hash,sub.hash.other)
                , by = c('v1', 'sub.hash.control' = 'sub.hash', 'sub.hash.other')
                , suffix = c('', '.control')
            ) %>%
            getDist[[metric]](Sub = sx) %>%
            filter(!is.na(d), !is.na(d.control)) %>%
            mutate(delta = d - d.control)

        g = triplet %>% 
            mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
            mutate(pairType = paste(v1,v2, sep = ', ')) %>%
            ggplot(aes(x = delta, fill = delta > 0))+
            geom_vline(xintercept = 0, linetype = 3)+
            geom_histogram(
                aes(group = pairType, fill = pairType)
                , breaks = getBreaks(triplet$delta)
                , alpha = 0.4
                , position = "identity"
            )+
            scale_fill_brewer(palette = "Dark2")+
            scale_x_continuous(
                limits = c(min(triplet$delta-0.1), max(triplet$delta+0.1))
            )+
            ylab('Num. triplets')+
            theme(
                axis.title.x = element_blank()
                , legend.title = element_blank()
                , legend.position = "top"
                , legend.justification = c(1,1)
                , legend.key.size = unit(0.5, 'lines')
                , legend.text = element_text(size = rel(0.8), hjust = 1)
            )+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            ggtitle(paste0(
                protein, ': ',toupper(sx)
                , ' (', triplet %>% 
                    computePvalue %>%
                    with(pvalue)
                , ')'
            ))
        return(g)
}


split(P, P$Group) %>%
lapply(function(Pg){
    SubGroup = Pg$Group[1]
    P.pair = P.pair %>%
        mutate(Group = GroupLabels[Group]) %>%
        filter(
            Group == SubGroup
            , Effect == "nonzero"
        ) %>%
        unnest(cols = p)
    if(nrow(P.pair)==0){ return(NULL) }
    P.pair = P.pair %>%
        select(sub, v1, v2, pvalue, dLower, dMid, dUpper)
    
    Pg = Pg %>%
        filter(Effect == "nonzero") %>%
        select_at(vars(colnames(P.pair))) %>%
        rbind(P.pair) %>%
        mutate(
            pos = str_extract(sub, '[0-9]+') %>% as.integer
            , v1v2 = paste(v1, v2, sep = ',')
        ) %>%
        arrange(pos, sub, v1v2) %>%
        mutate(
            sub = factor(sub, levels = unique(sub))
            , v1v2 = factor(v1v2, levels = unique(v1v2) %>% sort)
        ) %>%
        group_by(sub) %>%
        mutate(v1v2.nudge = ifelse(v1v2 == 'All,All'
            , 0
            , (rank(v1v2)-1) / (length(unique(v1v2))+2)
        ))
    gAll = Pg %>%
        mutate(pvalue = cut(pvalue, breaks = c(0,0.05, 0.1, 1), include.lowest = T)) %>%
        ggplot(aes(
            x = dMid, y = as.integer(sub) - v1v2.nudge
            , size = v1v2 == 'All,All'
            , shape = v1v2 == 'All,All'
            , color = pvalue
        ))+
        geom_vline(xintercept = 0, linetype = 3)+
        geom_linerange(aes(xmin = dLower, xmax = dUpper))+
        geom_point(size = 2)+
        scale_size_manual(values = c(0.3, 1), guide = 'none')+
        scale_shape_manual(values = c(1, 16), guide = 'none')+
        scale_y_continuous(
            breaks = function(x) setdiff(ceiling(x[1]):floor(x[2]), 0)
            , labels = levels(Pg$sub) %>% toupper
        )+
        scale_color_manual('P-value'
            , values = c('[0,0.05]'='red', '(0.05,0.1]'='orange', '(0.1,1]'='grey')
        )+
        ylab(paste0('Substitution in ', ifelse(SubGroup=="NS2A", "NS2A", "E")))+
        xlab('Difference in antigenic distance')
        
    # Delta distributions for subs with p<=0.1
    sub.list = Pg %>%
        filter(pvalue <= 0.1) %>%
        arrange(pvalue, v1v2) %>%
        with(unique(sub))
    if(length(sub.list) > 0){
        gDist = lapply(sub.list, plotPvalueDistribution, protein = ifelse(SubGroup=="NS2A", "NS2A", "E"))
        g = ggarrange(plotlist = gDist
                , nrow = ceiling(length(gDist)/3)
                , ncol = min(3, length(gDist))
            ) %>%
            annotate_figure(bottom = "Difference in antigenic distance")        
    } else {
        g = ggplot() + theme_void()
        gDist = list()
    }
    g = ggarrange(gAll, g
        , nrow = 2
        , ncol = 1
        , heights = c(1.8, ceiling(length(gDist)/3))
        , labels = ifelse(length(gDist)>0, 'auto', '')
        , font.label = list(size = 20)
    )
    ggsave(g
        , filename = file.path(outdir, paste0('delta_', GroupLabels.inv[SubGroup],'.pdf'))
        , width = 10
        , height = max(3, (nlevels(Pg$sub) + ceiling(length(gDist)/3)) * 0.6)
    )
})


    #   Distribution of p-value on phylogeny
    #   ....................................

plotdir = file.path(outdir, 'tree')
dir.create(plotdir, recursive = T)
mapply(function(x, Grp, Eff){
        if(nrow(x)==0){ return(NULL) }
        prot = ifelse(Grp=='NS2A', 'NS2A', 'E')
        g = x %>%
            filter(!is.na(pvalue)) %>%
            filter(N.v2 > 5) %>%
            with(sub) %>%
            unique %>%
            lapply(plotTreeFromSub, protein = prot) %>%
            do.call(what = c) %>%
            Filter(f = function(x) !is.null(x))
        if(length(g)==0){ return(NULL)}
        g = ggarrange(plotlist = g
            , ncol = 3
            , nrow = ceiling(length(g)/3)
            , common.legend = T
            , legend = 'bottom'
        )
        ggsave(g
            , filename = file.path(plotdir, paste0(Grp,'_',Eff,'.pdf'))
            , width = 8
            , height = 8 * ceiling(length(g)/3)/3
        )
    
    }
    , x = P$p[i]
    , Grp = P$Group[i]
    , Eff = P$Effect[i]
    , SIMPLIFY = F
)






    #   Example of how the p-value varies by v2
    #   .......................................

plotPvalueDistributionDetailed = function(sx, protein, metric = 'mapDist'){
    subs.full = left_join(subs[[c('E' = 'E.any', 'NS2A' = 'NS2A')[protein]]]
        , subs[[c('E' = 'NS2A', 'NS2A' = 'E.any')[protein]]] %>% select(v1,v2,sub.hash) %>% rename(sub.hash.other = sub.hash)
        , by = c('v1','v2')
    )

    triplet = subs.full %>%
        filter(SubPresent(sx, sub)) %>%
            mutate(
                sub = lapply(sub, function(sl){ setdiff(sl, sx) })
                , sub.hash.control = hashSubs(sub)
            ) %>%
            select(-sub) %>%
            inner_join(subs.full %>%
                    select(v1,v2,sub.hash,sub.hash.other)
                , by = c('v1', 'sub.hash.control' = 'sub.hash', 'sub.hash.other')
                , suffix = c('', '.control')
            ) %>%
            getDist[[metric]](Sub = sx) %>%
            filter(!is.na(d), !is.na(d.control)) %>%
            mutate(delta = d - d.control)
            
    g = ggarrange(
        triplet %>% 
            mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
            ggplot(aes(x = delta, fill = delta > 0))+
            geom_histogram(breaks = seq(min(triplet$delta-0.1), max(triplet$delta+0.1), by =0.1))+
            geom_text(data = triplet %>%
                mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
                group_by(v1, v2) %>%
                summarize(pvalue = mean(delta <= 0)) %>%
                mutate(pvalue = ifelse(pvalue<0.001, "p<0.001", sprintf("p=%.3f", pvalue))) %>%
                mutate(delta = min(triplet$delta))
                , x = -0.1, y = 11
                , aes(label = pvalue)
                , hjust = 1
            )+
            scale_fill_manual(values = c('black', 'red'), guide = 'none')+
            scale_x_continuous(
                position = 'top'
                , limits = c(min(triplet$delta-0.1), max(triplet$delta+0.1))
            )+
            ylab('Num. triplets')+
            theme(
                axis.title.x = element_blank()
            )
        
        , triplet %>%
            mutate(v2name = v2) %>%
            mutate_at(vars(matches('^v[12]$')), function(x) substring(x,1,5)) %>%
            group_by(v1,v2,v2name) %>%
            mutate(delta = d - d.control) %>%
            summarize(
                Ntriplet = n()
                , pvalue = mean(delta <= 0)
                , pvalue.txt = ifelse(pvalue<0.001, "p<0.001", sprintf("p=%.3f", pvalue))
                , dLower = quantile(delta, 0.025)
                , dMid = quantile(delta, 0.5)
                , dUpper = quantile(delta, 0.975)
            ) %>%
            ggplot(aes(y = v2name, color = pvalue <= 0.05))+
            geom_vline(xintercept = 0, linetype = 3, size = 0.5)+
            geom_linerange(aes(xmin = dLower, xmax = dUpper))+
            geom_point(aes(x = dMid))+
            geom_text(aes(x = dMid, label = pvalue.txt), nudge_y = rel(0.35), size = rel(2.3))+
            scale_color_manual(values = c('FALSE'='grey', 'TRUE'='red'), guide = 'none')+
            scale_x_continuous(
                limits = c(min(triplet$delta-0.1), max(triplet$delta+0.1))
            )+
            xlab('Difference in antigenic distance')+
            ylab('Virus with substitution of interest')+
            theme(
                axis.text.y = element_blank()
            )

        , nrow = 2
        , ncol = 1
        , align = 'v'
        , heights = c(1,2)
        , labels = 'auto'
    ) %>%
    annotate_figure(top = paste0(protein, ': ',toupper(sx)))
    return(g)
}

# histogram
plotPvalueDistributionDetailed(
    sx = 'm160k'  # in hmAb
    , protein = 'E'
)
ggsave(
    filename = file.path(outdir, 'pvalue_withinTypePair.pdf')
    , width = 6
    , height = 8
)

# tree
plotTreeFromSub(
    sx = 'm160k'  # in hmAb
    , protein = 'E'
)
ggsave(
    filename = file.path(outdir, 'pvalue_withinTypePair_tree.pdf')
    , width = 6
    , height = 4
)
