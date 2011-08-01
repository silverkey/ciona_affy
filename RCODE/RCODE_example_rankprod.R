library(limma)
library(affy)
library(RankProd)
targets = readTargets("target1.txt")
data = ReadAffy(filenames=targets$FileName, cdfname='cint06a520380fcdf')
eset = rma(data)

# Load annotations
ann.table = 'CIONA_AFFY_ANNOTATION_SUMMARY_SMALL.xls'
ann = read.table(file=ann.table, quote="", sep="\t", comment.char="", header=T)
ann$design_length = nchar(as.vector(ann$design))
# Let's give a look to the length of the design sequences
design.length = nchar(as.vector(ann$design))
short = as.vector(ann[design.length<250, 'probe_id'])
long = as.vector(ann[design.length>=250, 'probe_id'])

# Aspecific filters
eset = eset[long, ]
library(genefilter)
f1 <- pOverA(0.5, median(exprs(eset)))
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
eset <- eset[selected, ]

# The rank product analysis for single-origin 
# data can be performed by either RP or RPadvance.
# Differentially expressed genes between class 2 (class lable=1)and class 1 (class lable=0).
exprs = exprs(eset)
cl = c(0,0,1,1,1)
RP.out = RP(exprs,cl,num.perm=100,logged=T,gene.names=rownames(exprs),rand=123)
res = topGene(RP.out,cutoff=0.05,method="pfp",logged=TRUE,logbase=2,gene.names=rownames(exprs))
res = data.frame(rbind(res$Table1,res$Table2))
res = res[,c(3,4,5)]
colnames(res) = c('FC','FDR','P.value')
x = res$FC
res$FC = ifelse(x>1, x, -(1/x))
res$ID = rownames(res)

res.ann = merge(res[,c(4,1,2,3)], ann[,c(1,8,2,3,4,5,6,7)], by.x='ID', by.y='probe_id', all.x=T)
res.ann.sort = res.ann[order(res.ann$FDR, decreasing=F),]

# To avoid Excel complaints if the table starts with 'ID'
cn = colnames(res.ann.sort)
cn[1] = 'probe_id'
colnames(res.ann.sort) = cn

# Write outputs
write.table(res.ann.sort, file='RES.xls', sep="\t", quote=F, row.names=F)

# Make a table with the expression values
exprs = exprs(rma(data))
exprs = as.data.frame(exprs)
exprs$probe_id = rownames(exprs)
cord = c(length(colnames(exprs)), 2:length(colnames(exprs))-1)
write.table(exprs[,cord], file='expression.xls', sep="\t", quote=F, row.names=F)

