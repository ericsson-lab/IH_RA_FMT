wd = "/Users/z_mcadams/Desktop/"
setwd(wd)
list.files(wd)

library(xfun)

file$Confidence = NULL

file = 'taxonomy.tsv'
gsub_file(file, 'd__Bacteria; p', 'p')
gsub_file(file, '; c__.*;', '')
gsub_files(file, ' .*','')
gsub_files(file, ';', '')

readLines(file)

write.table(file, "aldex2.tsv", sep="\t)
