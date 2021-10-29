
    #   Configuration to harmonize plots across all analyses



Colors = new.env()

Colors$serotype$max = c(
    "DENV1" = "#fc3d03" # red
    , "DENV2" = "#14ad00" # green
    , "DENV3" =  "#8c03fc" # purple
    , "DENV4" = "#fcba03" # yellow
)
Colors$serotype$min = c(
    "DENV1" = "#f0aaaa" # red
    , "DENV2" = "#c6f584" # green
    , "DENV3" = "#8892d1" # purple
    , "DENV4" = "#f7e092" # yellow
)

Colors$Serotype = Colors$serotype$max


Colors$aaSource = c(
    'Actual' = 'red'
    , 'Polyprotein' = '#8f8f8f' # grey
    , 'Envelope' = '#036bfc' # blue
)
