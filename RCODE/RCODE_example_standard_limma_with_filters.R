library(limma)
library(affy)
targets = readTargets("targets.txt")
data = ReadAffy(filenames=targets$FileName, cdfname='cint06a520380fcdf')
eset = rma(data)
exprs = exprs(eset)

# Make a table with the expression values
exp = as.data.frame(exprs)
exp$probe_id = rownames(exp)
cord = c(length(colnames(exp)), 2:length(colnames(exp))-1)
write.table(exp[,cord], file='expression.xls', sep="\t", quote=F, row.names=F)

# Load annotations
ann.table = 'CIONA_AFFY_ANNOTATION_SUMMARY_SMALL.xls'
ann = read.table(file=ann.table, quote="", sep="\t", comment.char="", header=T)
ann$design_length = nchar(as.vector(ann$design))

# Let's give a look to the length of the design sequences
design.length = nchar(as.vector(ann$design))
short = as.vector(ann[design.length<250, 'probe_id'])
long = as.vector(ann[design.length>=250, 'probe_id'])
pdf(file='probe_design_impact.pdf', paper='a4r', width=8.3, height=11.7, pointsize=8)
boxplot(exprs[long,], boxwex=0.25, at=1:ncol(exprs)-0.2, col=rep('darkblue',ncol(exprs)), outline=F,
        main='Signals of long (darkblue) and short (lightblue) probe design')
boxplot(exprs[short,], boxwex=0.25, at=1:ncol(exprs)+0.2, col=rep('lightblue',ncol(exprs)), outline=F, add=T, names=F)
dev.off()

# Aspecific filters
eset = eset[long, ]
library(genefilter)
f1 <- pOverA(0.5, median(exprs(eset))*2)
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(eset, ff)
sum(selected)
eset <- eset[selected, ]

# CHANGE HERE!!!!!!!!! #
design = cbind(WT=c(1,1,1,0,0),DW=c(0,0,0,1,1))
rownames(design) = targets$FileName 

fit = lmFit(eset, design)

# CHANGE HERE!!!!!!!!! #
cont.matrix = makeContrasts(DW-WT, levels=design)
fit = contrasts.fit(fit, cont.matrix)
fit = eBayes(fit)

res = topTable(fit,number=Inf)

# Count the numbers of genes
nrow(res)
nrow(res[res$P.Value<0.05,])
nrow(res[(res$logFC<(-log2(1.5)) | res$logFC>log2(1.5)) & res$P.Value<0.05,])
nrow(res[(res$logFC<(-log2(1.5)) | res$logFC>log2(1.5)) & res$adj.P.Val<0.05,])

############## OUTPUT GENERATION ##############
# Add unlogged FC for lazy users....
x = res$logFC
res$linearFC = ifelse(x>0, 2^x, -(2^-x))

# Merge annotations and results
res.ann = merge(res[,c(1,8,2,5,6)], ann[,c(1,8,2,3,4,5,6,7)], by.x='ID', by.y='probe_id', all.x=T)
res.ann.sort = res.ann[order(res.ann$P.Value, decreasing=F),]

# Make tables with results

# To avoid Excel complaints if the table starts with 'ID'
cn = colnames(res.ann.sort)
cn[1] = 'probe_id'
colnames(res.ann.sort) = cn

# Write outputs
write.table(res.ann.sort, file='TOTAL_RES.xls', sep="\t", quote=F, row.names=F)
write.table(res.ann.sort[1:100,], file='TOP_100_RES.xls', sep="\t", quote=F, row.names=F)
write.table(res.ann.sort[(res.ann.sort$logFC<(-log2(1.5)) | res.ann.sort$logFC>log2(1.5)) & res.ann.sort$P.Value<0.05,],
            file='FILTERED_FC1.5_P0.05_RES.xls', sep="\t", quote=F, row.names=F)
write.table(res.ann.sort[(res.ann.sort$logFC<(-log2(1.5)) | res.ann.sort$logFC>log2(1.5)) & res.ann.sort$adj.P.Val<0.05,],
            file='FILTERED_FC1.5_Padj0.05_RES.xls', sep="\t", quote=F, row.names=F)

# Make a table with the expression values
exprs = exprs(eset)
exprs = as.data.frame(exprs)
exprs$probe_id = rownames(exprs)
cord = c(length(colnames(exprs)), 2:length(colnames(exprs))-1)
write.table(exprs[,cord], file='expression.xls', sep="\t", quote=F, row.names=F)


