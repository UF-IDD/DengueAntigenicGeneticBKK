
library(tidyverse)

# load monoclonal antibody data
dbdir = "00-RawData/antibodyDB"
act = read_csv(paste0(dbdir,"/activities.csv"), show_col_types = FALSE)
epitope = read_csv(paste0(dbdir,"/epitopes.csv"), show_col_types = FALSE)
mab = read_csv(paste0(dbdir,"/mAbs.csv"), show_col_types = FALSE)

ab = Reduce( full_join, list(act, epitope, mab)) %>%
    filter(grepl("^E",Protein)) %>%
    filter(`Epitope Type` != "Predicted") %>%
    filter(Infection == "Primary") %>%
    rename(mAb = `mAb name`) %>%
    select(mAb, Host, Specificity, `Assay Type`, `Epitope Type`, Protein, Domain, Sequence) %>%
    unique

ab$positions = lapply(
    strsplit(ab$Sequence,",")
    ,function(x){
        iRange = grepl("-",x)
        x = strsplit(x[iRange],"-") %>%
            lapply(as.integer) %>%
            do.call(what=rbind) %>%
            rbind(
                matrix(c(as.integer(x[!iRange]), as.integer(x[!iRange])), ncol=2)
            ) %>%
            as.data.frame
        colnames(x) = paste0('x',1:2)
        x
})        

rm(list = setdiff(ls(), c("ab","act")))
