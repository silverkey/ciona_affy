#!/usr/bin/R

# Finally GOstat on Ciona intestinalis!!!

# Write the table GO_significant.xls in output

# The file of the selected genes. Must contains a single column with
# the affymetrix probe id of the differentially expressed genes
selected.filename = 'affy_selected.txt'

# The file containing the association of GO id to the probe id.
# It must contains 2 columns, the first with the affymetrix id and 
# the second with the GO id.
# This file is coming out using the script associate_GO_to_affy_id.pl
# on the annotation database
go.annotation.filename = 'GO_array_annotation.txt'

# The file with the definition of the GO id. It contains 2 columns
# the first containing the GO id and the second with the definitions
# This file come out using the script parse_standard_go_table.pl
# on the official GO file GO.terms_alt_ids
go.definition.filename = 'go_definition.txt'

sel = read.table(file=selected.filename,sep="\t",header=F,quote="")
ass = read.table(file=go.annotation.filename,sep="\t",header=F,quote="")
def = read.table(file=go.definition.filename,sep="\t",header=F,quote="")

sel.ass = ass[ass[,1] %in% sel[,1],]
classes = as.vector(unique(sel.ass[,2]))

tot.uni = length(as.vector(unique(ass[,1])))
tot.sel = length(as.vector(unique(sel[,1])))

res = NULL

for(i in 1:length(classes)) {

  class = classes[i]

  n.uni = nrow(ass[ass[,2]==class,])
  n.sel = nrow(sel.ass[sel.ass[,2]==class,])

  p.uni = n.uni/tot.uni*100
  p.sel = n.sel/tot.sel*100

  if(p.sel > p.uni) {
    probes = paste(sel.ass[sel.ass[,2]==class,1],collapse=',')
    class.def = as.character(def[def[,1]==class,2])
    p.value = prop.test(c(n.sel,n.uni),c(tot.sel,tot.uni),alternative='g')$p.value
    res = rbind(res,data.frame(GO.id=class,n.uni=n.uni,n.sel=n.sel,p.uni=p.uni,p.sel=p.sel,p.value=p.value,adj.p=1,definition=class.def,probes=probes))
  }
}

adj.pval = p.adjust(res$p.value)
res$adj.p = adj.pval

sig = res[res$adj.p<=0.05,]
sig = sig[order(sig$adj.p),]

write.table(sig,file='GO_significant.xls',sep="\t",quote=F,row.names=F)
