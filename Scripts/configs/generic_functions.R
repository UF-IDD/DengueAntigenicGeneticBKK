

write_csv_mkdir = function(x, outstem){
    outfile = paste0(outstem, ".csv")
    dir.create(dirname(outfile), recursive = T)
    write_csv(x, outfile)
}
