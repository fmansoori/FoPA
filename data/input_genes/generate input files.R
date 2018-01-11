#generate input files

library("KEGGdzPathwaysGEO")
library("KEGGandMetacoreDzPathwaysGEO")

set="GSE18842"
#data(list=set,package="KEGGandMetacoreDzPathwaysGEO")
data(list=set,package="KEGGdzPathwaysGEO")

x=get(set)
#Extract from the dataset the required info
exp=experimentData(x);
dataset= exp@name
dat.m=exprs(x)
ano=pData(x)
design= notes(exp)$design
annotation= paste(x@annotation,".db",sep="")
targetGeneSets= notes(exp)$targetGeneSets

esetm=dat.m
group=ano$Group
paired=design=="Paired"
Block=ano$Block
targetgs=targetGeneSets
annotation=annotation



require(limma)
if (!annotation %in% c("hgu133a.db", "hgu133plus2.db")) {
	stopifnot(require(annotation, character.only = TRUE))
}

if (!paired) {
	G = factor(group)
	designm <- model.matrix(~0 + G)
	colnames(designm) <- levels(G)
	fit <- lmFit(esetm, designm)
	cont.matrix <- makeContrasts(contrasts = "d-c", levels = designm)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	aT1 <- topTable(fit2,coef=1, number=dim(esetm)[1])
}else{
	G = group
	block = factor(Block)
	G = factor(G)
	designm <- model.matrix(~0 + G + block)
	colnames(designm) <- substr(colnames(designm), 2, 100)
	fit <- lmFit(esetm, designm)
	cont.matrix <- makeContrasts(contrasts = "d-c", levels = designm)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	#aT1 <- topTable(fit2, coef = 1, number = dim(esetm)[1])
	aT1 <- topTable(fit2,coef=1,number=dim(esetm)[1])
}

aT1$ID = rownames(aT1)
anpack = paste(unlist(strsplit(annotation, split = ".db")), 
	"ENTREZID", sep = "")
aT1 <- aT1[aT1$ID %in% keys(get(anpack)), ]
aT1$ENTREZID = unlist(as.list(get(anpack)[aT1$ID]))
aT1 <- aT1[!is.na(aT1$ENTREZID), ]
aT1 <- aT1[!duplicated(aT1$ENTREZID), ]

#drop genes not in any geneset
#drop from esetm all duplicate genes and genes not in the genesets
esetm=esetm[rownames(esetm)%in%aT1$ID,]
rownames(esetm)<-aT1$ENTREZID[match(rownames(esetm),aT1$ID)]


if(!paired){
G=factor(group)

design <- model.matrix(~0+G)
colnames(design) <- levels(G)
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts="d-c",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1<-topTable(fit2, coef=1, number=dim(esetm)[1])

}else{
G=group
block=factor(Block)
G=factor(G)
design <- model.matrix(~0+G+block)
colnames(design)<-substr(colnames(design),2,100)
fit <- lmFit(esetm, design)
cont.matrix <- makeContrasts(contrasts="d-c",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aT1<-topTable(fit2,coef=1, number=dim(esetm)[1])
}

aT1$ID=rownames(aT1)

#make deg and write in files
#tg1<-aT1[aT1$P.Val<0.05,]
tg1 <- aT1[aT1$adj.P.Val<0.05,]
#write tg1$ID tg1$logFC in desired file

#write all genes
if(dim(tg1)[1] < 5000){
	tg_all <- aT1[1:10000,]
}else{	
	tg_all <- aT1
}
#write tg_all$ID tg_all$t in desired file



