
library(argparse)   
library(tidyverse); theme_set(theme_classic())
library(ggpubr)

parser = ArgumentParser(description = "Plot SpydrPick output.")
parser$add_argument('infile')
parser$add_argument('outbase')

# test
inputArg = parser$parse_args(c(
    "05-coevolution/spydrpick/all.filtered_ge001maf_le015gf_gt01-lt05states.L5262n1946.spydrpick_couplings.1-based.521601edges"
    , "05-coevolution/spydrpick/all"
))

inputArg = parser$parse_args()
dir.create(dirname(inputArg$outbase), recursive = T)


Site = list(
    "E" = 841 + (1:(495*3)) - 1
    , "NS2A" = 3382 + (1:(219*3)) - 1
)


dat = read_delim(inputArg$infile
        , delim = " "
        , col_names = c("pos1", "pos2", "genome_distance", "ARACNE", "MI")
        , col_types = "iiiid"
    ) %>%
    filter(MI >= quantile(MI, 0.99)) 


posMax = max(dat$pos2)
g = ggarrange(
    dat %>% mutate(posB = pos2, posA = pos1) %>%
        rbind(dat %>% mutate(posA = pos2, posB = pos1)) %>%
        ggplot(aes(x = posA, y = posB))+
        geom_density_2d_filled(contour_var = "ndensity", adjust = 1/4)+
        geom_rect(
            xmin = min(Site$NS2A), xmax = max(Site$NS2A)
            , ymin = min(Site$E), ymax = max(Site$E)
            , fill = NA
            , color = "#a6ff00"
            , size = 1.2
        )+
        geom_rect(
            xmin = 0, xmax = posMax
            , ymin = min(Site$E), ymax = max(Site$E)
            , fill = NA
            , color = "#a6ff00"
            , size = 0.1
        )+
        scale_x_continuous(
            'NN position on DENV genome'
            , expand = c(0,0)
            , limits = c(0, posMax) + 0.5
        )+
        scale_y_continuous(
            'NN position on DENV genome'
            , expand = c(0,0)
            , limits = c(0, posMax) + 0.5
        )+
        labs(fill = "Scaled density")+
        theme(legend.position = 'top')+
        guides(fill = guide_legend(
            direction = "horizontal"
            , title.position = "bottom"
            , label.position="top"
            , label.hjust = 0.5
            , label.vjust = 0.5
            , label.theme = element_text(angle = 90)
            , nrow = 1
        )) 
    , dat %>%
        # filter for interactions involving E
        filter(pos1 %in% Site$E | pos2 %in% Site$E) %>%
        mutate(
            posE = ifelse(pos1 %in% Site$E, pos1, pos2)
            , posOther = ifelse(pos1 %in% Site$E, pos2, pos1)
        ) %>%
        filter(posOther %in% Site$NS2A) %>%
        # mutate(
        #     Gene = 2*(posOther %in% Site$NS2A) + (posOther %in% Site$E)
        #     , Gene = c('Others', 'E', 'NS2A')[Gene+1]
        # ) %>%
        mutate(posE = posE - min(Site$E) + 1) %>%
        mutate(posNS2A = posOther - min(Site$NS2A) + 1) %>%
        ggplot(aes(x = posNS2A, y = posE))+
        geom_density_2d_filled(contour_var = "ndensity", adjust = 1/4)+
        scale_x_continuous(
            'NN position on NS2A gene'
            , expand = c(0,0)
            , position = "top"
            , sec.axis = sec_axis(
                trans = function(x) x/3
                , name = "AA position on NS2A gene"
            )
        )+
        scale_y_continuous(
            'NN position on E gene'
            , expand = c(0,0)
            , sec.axis = sec_axis(
                trans = function(x) x/3
                , name = "AA position on E gene"
            )
        
        )+
        guides(fill = "none")
        # labs(fill = "Scaled density")
    , nrow = 1
    , ncol = 2
    , widths = c(1,1)
    # , align = 'v'
    # , common.legend = T
    # , legend = "right"
    , labels = 'auto'
)
ggsave(g
    , filename = paste0(inputArg$outbase,"_highestMI.pdf")
    , width = 8
    , height = 7
)
