args = commandArgs(trailingOnly=TRUE)
dataset1=args[1]
dataset2=args[2]
dkplrun1=args[3]
dkplrun2=args[4]
library(limma)
library(HTSanalyzeR)
permutationPvalue<-function(genelist,gs,permtime){
	seedn=15
        upgenepool=genelist[genelist>0]
	downgenepool=genelist[genelist<0]
	upgs=upgenepool[gs$subset]
	upgs=upgs[complete.cases(upgs)]
	bond=quantile(upgs,c(0.05,0.95))
	upgs=upgs[upgs>bond[1] & upgs<bond[2]]
	downgs=downgenepool[gs$subset]
	downgs=downgs[complete.cases(downgs)]
	bond=quantile(downgs,c(0.05,0.95))
	downgs=downgs[downgs>bond[1] & downgs<bond[2]]
	obv_up=sum(upgs);obv_down=sum(downgs)
	permv_up=c();permv_down=c()
	for (i in 1:permtime){
		permv_up=c(permv_up,sum(sample(upgenepool,length(upgs),replace=FALSE)))
		permv_down=c(permv_down,sum(sample(downgenepool,length(downgs),replace=FALSE)))
	}
	p_up=sum(permv_up>obv_up)/permtime
	p_down=sum(permv_down<obv_down)/permtime
	return (list(up=p_up,down=p_down))
}

#####################################################################
#contig based enrichment
x='log2FC'
title='enrichment of contigs'
event=read.table(paste0(dataset1),header=T,stringsAsFactors = F,sep='\t')
tags=read.table(paste0('overlap_contig'),header=T,stringsAsFactors = F)
names(tags)=c('tag','nbkmer','contig')
pdf('enrichment.pdf',wi=8,he=8)
par(mfrow=c(2,1))
event=event[order(event[x],decreasing = FALSE),]
gl=event[[x]] #generate a list rather than a dataframe
names(gl)=event$contig
gsc=list(subset=tags$contig)
row.names(event)=1:nrow(event)
idx=as.numeric(rownames(event[event$contig %in% tags$contig,]))
GSCpvalues=permutationPvalue(gl,gsc,2000)
GSCscores <- collectionGsea(collectionOfGeneSets=gsc, geneList=gl,exponent=1, nPermutations=1000, minGeneSetSize=15)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
#GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FC+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FC-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.5,main=paste0(title,'(p=',GSCpvalues,')'),xlab='log2FoldChange')
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')

###################################################################
#gene based enrichment
x='log2FoldChange'
title='enrichment of genes'
event=read.table(paste0(dkplrun1,'/gene_expression/',list.files(paste0(dkplrun1,'/gene_expression/'),pattern='*-DEGs\\.tsv$')),header=T,stringsAsFactors = F,sep='\t')
event=event[complete.cases(event),]
event=event[abs(event$log2FoldChange)>2 & event$padj<0.05,]
event=event[!is.na(event$padj),]
tags=read.table(paste0('overlap_gene'),header=T,stringsAsFactors = F)
names(tags)='id'
event$id=rownames(event)
event=event[order(event[x],decreasing = FALSE),]
gl=event[[x]] #generate a list rather than a dataframe
names(gl)=event$id
gsc=list(subset=tags$id)
row.names(event)=1:nrow(event)
idx=as.numeric(rownames(event[event$id %in% tags$id,]))
GSCpvalues <- permutationPvalue(gl,gsc,2000)
GSCscores <- collectionGsea(collectionOfGeneSets=gsc, geneList=gl,exponent=1, nPermutations=5000, minGeneSetSize=15)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
#GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FoldChange+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FoldChange-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.5,main=paste0(title,'(p=',GSCpvalues,')'),xlab=x)
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')
dev.off()

################################################################
#draw jaccard index table
dat=read.table('jaccardidx_input.txt',header=F,stringsAsFactors = F,sep=' ')
library(scales)
pdf('jaccardidx.pdf',wi=7,he=8)
plot(1:dim(dat)[1], 1:dim(dat)[1], xlim = c(0,6), ylim = c (0,dim(dat)[1]+1), type = "n",ann = F, bty = "n", xaxt = "n", yaxt ="n")
for (i in 1:dim(dat)[1]){
	v1=as.numeric(gsub(',','',dat[i,2]));v2=as.numeric(gsub(',','',dat[i,4]));jc=as.numeric(sub('%','',dat[i,3]))
  	m=jc/100*1.8;l=max(0,log(v1,10)/10*2);r=max(0,log(v2,10)/10*2)
    	rect(xleft = c(3-m-l,3-m,3+m), xright = c(3-m,3+m,3+m+r), 
	        ybottom = c(rep(dim(dat)[1]-i,3)),ytop = c(rep(dim(dat)[1]-i+1,3)),
		col=alpha(c(rainbow(12)[8],rainbow(12)[1],rainbow(12)[2]),0.7))
    	text(3-m-l-.1,dim(dat)[1]-i+.5,pos=4,dat[i,2])
        text(3+m+r+.1,dim(dat)[1]-i+.5,pos=2,dat[i,4])
        if (m>0.01){text(3,dim(dat)[1]-i+.5,dat[i,3])}
        text(3-m-l-.5,dim(dat)[1]-i+.5,(dat[i,1]),cex=1,font=1)
        text(3,14.5,('dataset1|dataset2'),cex=1.5,font=2)
}
dev.off()
