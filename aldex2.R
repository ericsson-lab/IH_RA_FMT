library(ALDEx2)

directory = "/Users/z_mcadams/Desktop/"
file = 'aldex_upload.csv'
setwd(directory)

conds = c(rep('RA_FMT',5), rep("IH_FMT",5), 
          rep("RA_FMT_PRO", 5), rep("IH_FMT_PRO", 5)
          )

## SET ROWNAMES
file = read.csv(file)
file2 <- file[,-1]
rownames(file2) <- file[,1]

x.all <- aldex(file2, conditions =  conds, mc.samples=128, test="kw", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=TRUE
)

# aldex.clr: generating Monte-Carlo instances and clr values
# operating in serial mode
# removed rows with sums equal to zero
# computing center with all features
# data format is OK
# dirichlet samples complete
# transformation complete
# aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)
# operating in serial mode

x <- aldex.clr(file2, conds, mc.samples=128, denom="all", verbose=F)
# operating in serial mode
# computing center with all features

x.kw = aldex.kw(x)
# operating in serial mode

write.table(x.kw, "kw_aldex2.tsv", sep="\t")
