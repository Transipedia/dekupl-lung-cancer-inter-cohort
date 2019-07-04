import os,sys
dataset1,dataset2=sys.argv[1],sys.argv[2]
dkplrun1,dkplrun2=sys.argv[3],sys.argv[4]
eve=['SNV','splice','split','lincRNA','intron','repeat','unmapped']
features=[]
for e in eve:
    if os.path.exists('Report_%s.txt'%e) is False or int(os.popen('wc -l Report_%s.txt'%e).readline().strip().split()[0])==0:continue
    with open('Report_%s.txt'%e)as f:
        pairs=list(map(lambda x:x.strip().split(),f))
    features+=[i[:2] for i in pairs]
dict_2luad={}
for i in features:
    dict_2luad[i[1]]=i[0]
shared_contig_tcga=list(set([i[0] for i in features]))
shared_contig_seo=list(set([i[1] for i in features]))

with open(dataset2)as f:
    seodat=list(map(lambda x:x.strip().split(),f))
with open(dataset1)as f:
    tcgadat=list(map(lambda x:x.strip().split(),f))
hit_idx_seo=seodat[0].index('nb_hit')
hit_idx_tcga=tcgadat[0].index('nb_hit')
outtcga=open('overlap_contig','w')
outtcga.write('\t'.join(tcgadat[0])+'\n')
c=0
for line in tcgadat:
    if line[2] in dict_2luad.values() and int(float(line[hit_idx_tcga]))==1:
        outtcga.write(line[2]+'\t'+'\t'.join(line[1:])+'\n')
outtcga.close()

import pandas as pd
degfile1=os.popen('ls %s/gene_expression/*DEGs.tsv'%dkplrun2).readline().strip()
with open(degfile1)as f:
    seodat=list(map(lambda x:x.strip().split(),f))
degfile2=os.popen('ls %s/gene_expression/*DEGs.tsv'%dkplrun1).readline().strip()
with open(degfile2)as f:
    tcgadat=list(map(lambda x:x.strip().split(),f))
tcgadat=[i for i in tcgadat if 'NA' not in i]
seodat=[i for i in seodat if 'NA' not in i]
seodeg=[i[0] for i in seodat[1:] if abs(float(i[2]))>2 and float(i[-1])<0.05]
tcgadeg=[i[0] for i in tcgadat[1:] if abs(float(i[2]))>2 and float(i[-1])<0.05]

shared_DEG=list(set(seodeg)&set(tcgadeg))
with open(dkplrun2+'/gene_expression/normalized_counts.tsv')as f:
    seocount=list(map(lambda x:x.strip().split(),f))
with open(dkplrun1+'/gene_expression/normalized_counts.tsv')as f:
    tcgacount=list(map(lambda x:x.strip().split(),f))
outtcga=open('overlap_gene','w')
outtcga.write('\t'.join(tcgacount[0])+'\n')
for i in tcgacount:
    if i[0] in shared_DEG:outtcga.write('\t'.join(i)+'\n')
outtcga.close()

