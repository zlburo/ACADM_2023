###Figure 2###
##ACADM diff
library(limma)
library(ggplot2)
library(ggpubr)

gene="ACADM"              
expFile="symbol.txt"    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=t(data[gene,,drop=F])

group=sapply(strsplit(rownames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       
treatNum=length(group[group==0])     
Type=c(rep(1,conNum), rep(2,treatNum))

exp=cbind(data, Type)
exp=as.data.frame(exp)
colnames(exp)=c("gene", "Type")
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
exp$gene=log2(exp$gene+1)

outTab=exp
colnames(outTab)=c(gene, "Type")
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp.txt", sep="\t", quote=F, row.names=F)

group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(exp, x="Type", y="gene", color="Type",
                  xlab="",
                  ylab=paste0(gene, " expression"),
                  legend.title="Type",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

##paire diff
library(limma)
library(ggpubr)
expFile="geneExp.txt"      

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
geneName=colnames(rt)[1]

normalData=rt[rt$Type=="Normal",1,drop=F]
normalData=as.matrix(normalData)
rownames(normalData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(normalData))
normalData=avereps(normalData)
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
tumorData=avereps(tumorData)
sameSample=intersect(row.names(normalData), row.names(tumorData))
data=cbind(normalData[sameSample,,drop=F], tumorData[sameSample,,drop=F])
colnames(data)=c("Normal", "Tumor")
data=as.data.frame(data)

pdf(file="pairDiff.pdf", width=5, height=4.5)
ggpaired(data, cond1="Normal", cond2="Tumor", fill="condition",
         xlab="", ylab=paste0(geneName, " expression"),
         legend.title="Type",
         palette=c("blue","red"))+
  #stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
  stat_compare_means(paired = TRUE, symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",label.x = 1.35)
dev.off()

##data combat
library(limma)
library(sva)
rt=read.table("GEO15641.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo15641=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo15641=avereps(geo15641)

#log2
qx=as.numeric(quantile(geo15641, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo15641[geo15641<0]=0
  geo15641=log2(geo15641+1)}
geo15641=normalizeBetweenArrays(geo15641)

##36895
rt=read.table("GEO36895.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo36895=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo36895=avereps(geo36895)

#log2
qx=as.numeric(quantile(geo36895, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo36895[geo36895<0]=0
  geo36895=log2(geo36895+1)}
geo36895=normalizeBetweenArrays(geo36895)

###GEO46699
rt=read.table("GEO46699.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo46699=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo46699=avereps(geo46699)

#log2
qx=as.numeric(quantile(geo46699, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo46699[geo46699<0]=0
  geo46699=log2(geo46699+1)}
geo46699=normalizeBetweenArrays(geo46699)

###GEO53000
rt=read.table("GEO53000.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo53000=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo53000=avereps(geo53000)

#log2
qx=as.numeric(quantile(geo53000, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo53000[geo53000<0]=0
  geo53000=log2(geo53000+1)}
geo53000=normalizeBetweenArrays(geo53000)

###GEO53757
rt=read.table("GEO53757.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo53757=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo53757=avereps(geo53757)

#log2
qx=as.numeric(quantile(geo53757, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo53757[geo53757<0]=0
  geo53757=log2(geo53757+1)}
geo53757=normalizeBetweenArrays(geo53757)

###ICGC
rt=read.table("ICGC.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
ICGC=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
ICGC=avereps(ICGC)

#log2
qx=as.numeric(quantile(ICGC, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  ICGC[ICGC<0]=0
  ICGC=log2(ICGC+1)}
ICGC=normalizeBetweenArrays(ICGC)


sameGene=Reduce(intersect,list(rownames(geo15641),
                               rownames(geo36895),
                               rownames(geo46699),
                               rownames(geo53000),
                               rownames(geo53757),rownames(ICGC)))

geo15641Out=geo15641[sameGene,]
geo36895Out=geo36895[sameGene,]
geo46699Out=geo46699[sameGene,]
geo53000Out=geo53000[sameGene,]
geo53757Out=geo53757[sameGene,]
ICGCOut=ICGC[sameGene,]


all=cbind(geo15641Out,geo36895Out,geo46699Out,geo53000Out,geo53757Out,ICGCOut)
batchType=c(rep(1,ncol(geo15641Out)),rep(2,ncol(geo36895Out)),rep(2,ncol(geo46699Out)),rep(2,ncol(geo53000Out)),rep(2,ncol(geo53757Out)),rep(2,ncol(ICGCOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)

geo15641Out=outTab[,colnames(geo15641Out)]
geo15641Out[geo15641Out<0]=0

geo36895Out=outTab[,colnames(geo36895Out)]
geo36895Out[geo36895Out<0]=0

geo46699Out=outTab[,colnames(geo46699Out)]
geo46699Out[geo46699Out<0]=0

geo53000Out=outTab[,colnames(geo53000Out)]
geo53000Out[geo53000Out<0]=0

geo53757Out=outTab[,colnames(geo53757Out)]
geo53757Out[geo53757Out<0]=0

ICGCOut=outTab[,colnames(ICGCOut)]
ICGCOut[ICGCOut<0]=0

geo15641Tab=rbind(ID=colnames(geo15641Out), geo15641Out)
write.table(geo15641Tab, file="geo15641.normalize.txt", sep="\t", quote=F, col.names=F)

geo36895Tab=rbind(ID=colnames(geo36895Out), geo36895Out)
write.table(geo36895Tab,file="geo36895.normalize.txt",sep="\t",quote=F,col.names=F)

geo46699Tab=rbind(ID=colnames(geo46699Out), geo46699Out)
write.table(geo46699Tab,file="geo46699.normalize.txt",sep="\t",quote=F,col.names=F)

geo53000Tab=rbind(ID=colnames(geo53000Out), geo53000Out)
write.table(geo53000Tab,file="geo53000.normalize.txt",sep="\t",quote=F,col.names=F)


geo53757Tab=rbind(ID=colnames(geo53757Out), geo53757Out)
write.table(geo53757Tab,file="geo53757.normalize.txt",sep="\t",quote=F,col.names=F)

ICGCTab=rbind(ID=colnames(ICGCOut), ICGCOut)
write.table(ICGCTab,file="ICGC.normalize.txt",sep="\t",quote=F,col.names=F)

###Figure 3###
##cliHeatmap-Figure 3A
library(limma)
library(ComplexHeatmap)
expFile="geneExp.txt"       
cliFile="clincial without unknown.txt"     

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
Type=factor(Type, levels=c("Low","High"))
data=cbind(as.data.frame(data), Type)
data=data[order(data[,gene]),] 

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,"Type",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

bioCol=c("#0066FF","#FF9900","#FF0000","#ed1299", "#0dbc21", "#246b93", "#cc8e12", "#d561dd", "#c93f00", 
         "#ce2523", "#f7aa5d", "#9ed84e", "#39ba30", "#6ad157", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
         "#1a918f", "#7149af", "#ff66fc", "#2927c4", "#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
         "#4aef7b", "#e86502",  "#99db27", "#e07233", "#8249aa","#cebb10", "#03827f", "#931635", "#ff523f",
         "#edd05e", "#6f25e8", "#0dbc21", "#167275", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1",
         "#dd27ce", "#07a301", "#ddd53e",  "#391c82", "#2baeb5","#925bea", "#09f9f5",  "#63ff4f")
colorList=list()
colorList[[gene]]=c("Low"="blue", "High"="red")
j=0
for(cli in colnames(rt[,2:ncol(rt)])){
  cliLength=length(levels(factor(rt[,cli])))
  cliCol=bioCol[(j+1):(j+cliLength)]
  j=j+cliLength
  names(cliCol)=levels(factor(rt[,cli]))
  cliCol["unknow"]="grey75"
  colorList[[cli]]=cliCol
}

ha=HeatmapAnnotation(df=rt, col=colorList)
zero_row_mat=matrix(nrow=0, ncol=nrow(rt))
Hm=Heatmap(zero_row_mat, top_annotation=ha)

pdf(file="heatmap.pdf", width=7, height=5)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
##cliCor-Figure 3B-H
library(limma)
library(ggpubr)

expFile="geneExp.txt"       #?????????ļ?
cliFile="clinical1.txt"      #?ٴ??????ļ?

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)

for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c(gene, clinical)]
  colnames(data)=c(gene, "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                    xlab=clinical,
                    ylab=paste(gene, " expression"),
                    legend.title=clinical)+ 
    stat_compare_means(comparisons = my_comparisons)
  pdf(file=paste0("clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}

###Figure 4
library(limma)
library(ggpubr)
library(survival)
library(survminer)
risk=read.table("ACADM-Risk.txt",header=T,sep="\t",check.names=F,row.names=1)     
cli=read.table("clinical3.txt",sep="\t",check.names=F,header=T,row.names=1)      
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
data=cbind(futime=risk[,1],fustat=risk[,2],cli,risk=risk[,"Risk"])
for(i in colnames(data[,3:(ncol(data)-1)])){
  rt=data[,c("futime","fustat",i,"risk")]
  rt=rt[(rt[,i]!="unknow"),]
  colnames(rt)=c("futime","fustat","clinical","risk")
  tab=table(rt[,"clinical"])
  tab=tab[tab!=0]
  for(j in names(tab)){
    rt1=rt[(rt[,"clinical"]==j),]
    tab1=table(rt1[,"risk"])
    tab1=tab1[tab1!=0]
    labels=names(tab1)
    if(length(labels)==2){
      titleName=j
      if((i=="age") | (i=="Age") | (i=="AGE")){
        titleName=paste0("age",j)
      }
      diff=survdiff(Surv(futime, fustat) ~risk,data = rt1)
      pValue=1-pchisq(diff$chisq,df=1)
      if(pValue<0.001){
        pValue="p<0.001"
      }else{
        pValue=paste0("p=",sprintf("%.03f",pValue))
      }
      fit <- survfit(Surv(futime, fustat) ~ risk, data = rt1)
      surPlot=ggsurvplot(fit, 
                         data=rt1,
                         conf.int=F,
                         pval=pValue,
                         pval.size=6,
                         title=paste0("Patients with ",titleName),
                         legend.title="Risk",
                         legend.labs=labels,
                         font.legend=12,
                         xlab="Time(years)",
                         break.time.by = 1,
                         palette=c("red", "blue"),
                         risk.table=TRUE,
                         risk.table.title="",
                         risk.table.col = "strata",
                         risk.table.height=.25)
      j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
      pdf(file=paste0("survival.",i,"_",j,".pdf"),onefile = FALSE,
          width = 6,       
          height =5)      
      print(surPlot)
      dev.off()
    }
  }
}

####Figure 5####
library(survival)
library(regplot)
library(rms)

expFile="ACADM-Risk.txt"     
cliFile="clinical3.txt"    
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

samSample=intersect(row.names(exp), row.names(cli))
exp1=exp[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(exp1, cli)

res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[2,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(exp1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

pdf(file="calibration.pdf", width=5, height=5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

####Figure 7####
###correlated genes analysis
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

gene="ACADM"               
corFilter=0.5            
pFilter=0.05            
expFile="symbol.txt"    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=log2(data+1)
x=as.numeric(data[gene,])
outTab=data.frame()
for(j in rownames(data)){
  if(gene==j){next}
  y=as.numeric(data[j,])
  corT=cor.test(x, y, method = 'pearson')
  cor=corT$estimate
  pvalue=corT$p.value
  outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
  
}

outTab=outTab[abs(as.numeric(outTab$cor))>corFilter & as.numeric(outTab$pvalue)<pFilter,]
write.table(file="Correlated genes.txt", outTab, sep="\t", quote=F, row.names=F)

###GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05       
qvalueFilter=0.05      

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")


rt=read.table("intersectgenes.txt", header=T, sep="\t", check.names=F)     #??ȡ?????ļ?

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

showNum=5
if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO-barplot.pdf", width=9, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

###KEGG analysis
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05     
qvalueFilter=0.05    
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
rt=read.table("intersectgenes.txt", header=T, sep="\t", check.names=F)    
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=15
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGG-barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()

####Figure 8####
###CIBERSORT
library("limma")       
expFile="symbol.txt"    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

v=voom(data, plot=F, save.plot=F)
out=v$E
out=rbind(ID=colnames(out), out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #?????ļ?

source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

###Lollipop analysis
inputFile="cor.result.txt"      

data = read.table(inputFile, header=T, sep="\t", check.names=F)

p.col = c('gold','pink','orange','LimeGreen','darkgreen')
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                             ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)
}

p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                           ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)
}

points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color

points.cex = fcex(x=data$cor)
data$points.cex = points.cex
data=data[order(data$cor),]

xlim = ceiling(max(abs(data$cor))*10)/10         
pdf(file="Lollipop.pdf", width=9, height=7)    
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)

segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)

points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()

###Immune cell diff analysis
library(limma)
library(reshape2)
library(ggpubr)
library(vioplot)
library(ggExtra)
expFile="geneExp.txt"             
immFile="CIBERSORT-Results.txt"   
pFilter=0.05           
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

data=as.data.frame(data)
data$gene=ifelse(data[,gene]>median(data[,gene]), "High", "Low")

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

sameSample=intersect(row.names(immune), row.names(data))
rt=cbind(immune[sameSample,,drop=F], data[sameSample,,drop=F])

data=rt[,-(ncol(rt)-1)]
data=melt(data,id.vars=c("gene"))
colnames(data)=c("gene", "Immune", "Expression")

group=levels(factor(data$gene))
data$gene=factor(data$gene, levels=c("Low","High"))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="gene",
                  xlab="",
                  ylab="Fraction",
                  legend.title=gene,
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=gene),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="immune.diff.pdf", width=7, height=6)
print(boxplot)
dev.off()

outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-2)]){
  x=as.numeric(rt[,gene])
  y=as.numeric(rt[,i])
  if(sd(y)==0){y[1]=0.00001}
  cor=cor.test(x, y, method="spearman")
  outVector=cbind(Cell=i, cor=cor$estimate, pvalue=cor$p.value)
  outTab=rbind(outTab,outVector)
}
write.table(outTab,file="cor.result.txt",sep="\t",row.names=F,quote=F)

###TIDE analysis
library(limma)
library(ggpubr)
tideFile="TIDE.txt"         
riskFile="ACADM-Risk.txt"    

tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "Risk", drop=F]
data=cbind(tide, risk)

data$Risk=ifelse(data$Risk=="High", "High-ACADM", "Low-ACADM")
group=levels(factor(data$Risk))
data$Risk=factor(data$Risk, levels=c("Low-ACADM", "High-ACADM"))
group=levels(factor(data$Risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
  gg1=ggviolin(data, x="Risk", y=i, fill = "Risk", 
               xlab="", ylab=i,
               palette=c("#0066FF","#FF0000"),
               legend.title="ACADM",
               add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
  
  pdf(file=paste0("violin.", i, ".pdf"), width=6, height=5)
  print(gg1)
  dev.off()
}
