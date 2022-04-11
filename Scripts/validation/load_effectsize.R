
library(tidyverse)

read_csv.silent = function(...){read_csv(..., show_col_types = FALSE)}

# import substitutions with nonzero effect sizes
posNonzero = list(
    E = "04-plots/E/effect_diversity.csv" %>%
        read_csv.silent %>%
        filter(nonzero==1) %>%
        with(pos)
    , "NS2A" = "04-plots/compare_adjusted/positions_NS2A_beyondE.txt" %>%
        readLines %>%
        as.integer
)
effect = "03-fits/nnls/thai_map/AA_NS2A_beyondE/adj_envelope/aasub/d_summary.csv" %>%
    read_csv.silent %>%
    mutate(
        protein_name = ifelse(pos <= 495, "E", "NS2A")
        , pos = sapply(pos, function(x){
            ifelse(x <= 495, x, posNonzero$NS2A[x - 495])
        })
        , sub = mapply(gsub
            , pattern = "[0-9]+"
            , replacement = pos
            , x = Feature
        )
    ) %>%
    select(-Feature)
effect = split(effect, effect$protein_name)


# effect size vectors, named by substitution
effectMid = lapply(effect, function(x){
    x = x %>% filter(lower > 0)
    out = x$Median
    names(out) = x$sub
    out
})
