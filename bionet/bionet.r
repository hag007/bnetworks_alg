library("Rgraphviz", lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")
library("BioNet", lib.loc  = "/specific/netapp5/gaga/hagailevi/evaluation/Renv")

# network.file.name <- "/media/hag007/Data/bnet/networks/dip.sif"
# deg.file.name <- "/media/hag007/Data/bnet/datasets/GWAS_2hr_glucose/data/score.tsv"
# fdr=0.05
# is.pval.score=T

##load DIP ppita
ig <- loadNetwork.sif(network.file.name)
data <- read.delim(deg.file.name, sep="\t", stringsAsFactors=FALSE, row.names = 1)
pval <- data[["pval"]]
names(pval) <- rownames(data)
combined.dataset <- data
subnet <- subNetwork(rownames(combined.dataset), ig)
subnet <- rmSelfLoops(subnet)

pval <-na.omit(pval)
# pval <- pval[!pval==1]
if (is.pval.score){
	pval[pval < 1e-150] <- 1e-150
	fb <- fitBumModel(pval, plot = FALSE, starts = 10)
	scores <- scoreNodes(subnet, fb, fdr = fdr)
} else{
	scores <- pval
}

if (max(scores)<0){
  module.genes=list()  
} else {
  # scores = -log10(pval)
  module <- runFastHeinz(subnet, scores)
  # logFC <- dataLym$diff
  # names(logFC) <- dataLym$label
  # plotModule(module, scores = scores) # , diff.expr = logFC
  module.genes <- nodes(module)
}
bg.genes <- nodes(subnet)
