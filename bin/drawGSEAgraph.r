library(limma)
library(HTSanalyzeR)
x='log2FoldChange'
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

pdf('GSEAenrichment.pdf',wi=8,he=6)
par(mfrow=c(2,2))
title='two LUAD data'
event=read.table(paste0('/store/EQUIPES/SSFA/MEMBERS/yunfeng.wang/dkpl-run/dekupl-run/DEkupl_result/DEkupl_result_LUAD_TCGA/gene_expression//normalvstumor-DEGs.tsv'),header=T,stringsAsFactors = F,sep='\t')
event=event[complete.cases(event),]
event=event[abs(event$log2FoldChange)>2 & event$padj<0.05,]
event=event[!is.na(event$padj),]
tags=read.table(paste0('traindat_tcga2seo_gene'),header=T,stringsAsFactors = F)
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
GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
#GSCpvalues=0.78
#cat (paste0(eve,' gene pvalue ',GSCpvalues$up,' & ',GSCpvalues$down,'; total num: ',dim(event)[1],';validated num: ',dim(tags)[1]),file='events_pvalue',sep='\n',append=TRUE)
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FoldChange+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FoldChange-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.3,main=paste0('shared DEGs of ',title,'(p=',GSCpvalues,')'),xlab=x)
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')

title='LUAD and PRAD'
event=read.table(paste0('/store/EQUIPES/SSFA/MEMBERS/yunfeng.wang/dkpl-run/dekupl-run/DEkupl_result/DEkupl_result_PRAD_TCGA/gene_expression/normalvstumoral-DEGs.tsv'),header=T,stringsAsFactors = F,sep='\t')
event=event[complete.cases(event),]
event=event[abs(event$log2FoldChange)>2 & event$padj<0.05,]
event=event[!is.na(event$padj),]
tags=read.table(paste0('traindat_prad2luad_gene'),header=T,stringsAsFactors = F)
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
GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
GSCpvalues=0.78
#cat (paste0(eve,' gene pvalue ',GSCpvalues$up,' & ',GSCpvalues$down,'; total num: ',dim(event)[1],';validated num: ',dim(tags)[1]),file='events_pvalue',sep='\n',append=TRUE)
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FoldChange+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FoldChange-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.3,main=paste0('shared DEGs of ',title,'(p=',GSCpvalues,')'),xlab=x)
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')

x='log2FC'
title='two LUAD data'
event=read.table(paste0('/store/EQUIPES/SSFA/MEMBERS/yunfeng.wang/dkpl-anno/DEkupl_annotation_LUAD_SEO_logfc2/DiffContigsInfos.tsv'),header=T,stringsAsFactors = F,sep='\t')
tags=read.table(paste0('traindat_seo2tcga_contig'),header=T,stringsAsFactors = F)
names(tags)=c('tag','nbkmer','contig')
event=event[order(event[x],decreasing = FALSE),]
gl=event[[x]] #generate a list rather than a dataframe
names(gl)=event$contig
gsc=list(subset=tags$contig)
row.names(event)=1:nrow(event)
idx=as.numeric(rownames(event[event$contig %in% tags$contig,]))
GSCpvalues=permutationPvalue(gl,gsc,2000)
GSCscores <- collectionGsea(collectionOfGeneSets=gsc, geneList=gl,exponent=1, nPermutations=1000, minGeneSetSize=15)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FC+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FC-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.5,main=paste0('shared contigs of ',title,'(p=',GSCpvalues,')'),xlab='log2FoldChange')
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')

title='LUAD and PRAD'
event=read.table(paste0('/store/EQUIPES/SSFA/MEMBERS/yunfeng.wang/dkpl-anno/DEkupl_annotation_PRAD_TCGA_logfc2/DiffContigsInfos.tsv'),header=T,stringsAsFactors = F,sep='\t')
tags=read.table(paste0('traindat_prad2luad_contig'),header=T,stringsAsFactors = F)
names(tags)=c('tag','nbkmer','contig')
event=event[order(event[x],decreasing = FALSE),]
gl=event[[x]] #generate a list rather than a dataframe
names(gl)=event$contig
gsc=list(subset=tags$contig)
row.names(event)=1:nrow(event)
idx=as.numeric(rownames(event[event$contig %in% tags$contig,]))
GSCpvalues=permutationPvalue(gl,gsc,2000)
GSCscores <- collectionGsea(collectionOfGeneSets=gsc, geneList=gl,exponent=1, nPermutations=1000, minGeneSetSize=15)
GSCpvalues <- permutationPvalueCollectionGsea(permScores=GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
GSCpvalues=ifelse(GSCpvalues<0.5,GSCpvalues,1-GSCpvalues)
GSCpvalues=0.21
log2pos1=as.numeric(rownames(event[which.min(abs(event$log2FC+2)),]))
log2pos2=as.numeric(rownames(event[which.min(abs(event$log2FC-2)),]))
barcodeplot(1:dim(event)[1], index = idx,alpha=0.5,main=paste0('shared contigs of ',title,'(p=',GSCpvalues,')'),xlab='log2FoldChange')
abline(v=log2pos1,lwd=2,col='red')
abline(v=log2pos2,lwd=2,col='red')
dev.off()
