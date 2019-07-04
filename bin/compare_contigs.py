import os,sys,gc
import pandas as pd
import numpy as np
#prepare Events data
dataset1,dataset2,dkplrun1,dkplrun2,genomefadir=sys.argv[1:]
if os.path.exists(os.path.dirname(dataset1)+'/Events') is False or os.path.exists(os.path.dirname(dataset2)+'/Events') is False:
    print ('Generating Events')
    os.system('Rscript parsecontigs.r %s'%dataset1)
    os.system('Rscript parsecontigs.r %s'%dataset2)
else:
    print ('Events already exist')
if os.path.exists(os.path.dirname(genomefadir)+'/genome_db.nin'):pass
else:
    if genomefadir.endswith('gz'):os.system('gunzip -c %s > %s'%(genomefadir,genomefadir.strip('.gz')))
    os.system('makeblastdb -in %s -out %s -parse_seqids -dbtype nucl'%(genomefadir.rstrip('.gz'),os.path.dirname(genomefadir)+'/genome_db'))
    print ('makeblastdb -in %s -out %s -parse_seqids -dbtype nucl'%(genomefadir.rstrip('.gz'),os.path.dirname(genomefadir)+'/genome_db'))
genomedb=os.path.dirname(genomefadir)+'/genome_db'

def run_diff(eve):
    if eve=='gene':
        deg1file=os.popen('ls %s/gene_expression/*-DEGs.tsv'%dkplrun1).readline().strip()
        dat1=pd.read_csv(deg1file,header=0,index_col=0,sep='\t')
        deg2file=os.popen('ls %s/gene_expression/*-DEGs.tsv'%dkplrun2).readline().strip()
        dat2=pd.read_csv(deg2file,header=0,index_col=0,sep='\t')
        deg1=dat1[(dat1.padj<0.05) & (abs(dat1.log2FoldChange)>2)]
        deg2=dat2[(dat2.padj<0.05) & (abs(dat2.log2FoldChange)>2)]
        overlap=len(set(deg1.index)&set(deg2.index))
        print ('Diff_genes',deg1.shape[0],str(round(100*overlap/(deg1.shape[0]+deg2.shape[0]-overlap+1),2))+'%',deg2.shape[0])
        os.system('echo Diff_genes %s %s %s > jaccardidx_input.txt'%(deg1.shape[0],str(round(100*overlap/(deg1.shape[0]+deg2.shape[0]-overlap+1),2))+'%',deg2.shape[0]))
        del dat1,dat2,deg1,deg2,overlap
        gc.collect()
    if eve=='kmer':
        deg1file=os.popen('ls %s/*_vs_*/diff-counts.tsv.gz'%dkplrun1).readline().strip()
        deg1=pd.read_csv(deg1file,header=0,index_col=0,compression='gzip',sep='\t')
        deg2file=os.popen('ls %s/*_vs_*/diff-counts.tsv.gz'%dkplrun2).readline().strip()
        deg2=pd.read_csv(deg2file,header=0,index_col=0,compression='gzip',sep='\t')
        overlap=len(set(deg1.index)&set(deg2.index))
        print ('Diff_kmers',deg1.shape[0],str(round(100*overlap/(deg1.shape[0]+deg2.shape[0]-overlap+1),2))+'%',deg2.shape[0])
        os.system('echo Diff_kmers %s %s %s >> jaccardidx_input.txt'%(deg1.shape[0],str(round(100*overlap/(deg1.shape[0]+deg2.shape[0]-overlap+1),2))+'%',deg2.shape[0]))
        del deg1,deg2,overlap
        gc.collect()

def run_noncode(eve):
    cutoff=0.8
    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset2),eve)) as f:
        SEO=list(map(lambda x:x.strip('\n').split('\t'),f))
    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset1),eve)) as f:
        LUAD=list(map(lambda x:x.strip('\n').split('\t'),f))

    luadd=dict(zip([i[-7] for i in LUAD],[i[-6] for i in LUAD]))
    srad=dict(zip([i[-7] for i in SEO],[i[-6] for i in SEO]))


    shared_contigs=[]
    out=open('Report_%s.txt'%eve,'w')
    for luad in LUAD[1:]:
        chrq,startq,endq=luad[2],int(luad[3]),int(luad[4])
        for sra in SEO[1:]:
            chrs,starts,ends=sra[2],int(sra[3]),int(sra[4])
            if eve in ['polyA','polyA_DU']:
                if chrq==chrs and abs(endq-ends)<10:
                    shared_contigs.append([luad[-6],sra[-6]])
                    out.write(luad[-6]+'\t'+sra[-6]+'\n')
            elif eve in ['lincRNA','intron','intron_DU']:
                if starts>endq or ends<startq:continue
                actlen=max(ends,endq)-min(starts,startq)
                extlen=ends-starts+endq-startq
                overlap=extlen-actlen
                if chrq==chrs and (overlap>30 or abs((endq+startq)/2-(ends+starts)/2)<30 or overlap>cutoff*(ends-starts) or overlap>cutoff*(endq-startq)):
                    shared_contigs.append([luad[-6],sra[-6]])
                    out.write(luad[-6]+'\t'+sra[-6]+'\n')#output contig
    out.close()
    shared_LUAD1,shared_SEO=len(set([i[0] for i in shared_contigs])),len(set([i[1] for i in shared_contigs]))

    print (eve,len(LUAD)-1,str(round(100*max(shared_LUAD1,shared_SEO)/(len(LUAD)+len(SEO)-2-max(shared_LUAD1,shared_SEO)+1),2))+'%',len(SEO)-1)
    os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,len(LUAD)-1,str(round(100*max(shared_LUAD1,shared_SEO)/(len(LUAD)+len(SEO)-2-max(shared_LUAD1,shared_SEO)+1),2))+'%',len(SEO)-1))
    #combineTooutput(eve)
    return (shared_LUAD1,shared_SEO)

def run_snp(eve):
    if os.path.exists('dataset1_%s_mutinfo'%eve) is False:
        os.system("python3 extractMutFromBam.py %s %s dataset1"%(eve,dataset1))
    if os.path.exists('dataset2_%s_mutinfo'%eve) is False:
        os.system("python3 extractMutFromBam.py %s %s dataset2"%(eve,dataset2))
    with open('dataset1_%s_mutinfo'%(eve))as f:
        LT=list(map(lambda x:x.strip().split(),f))
    with open('dataset2_%s_mutinfo'%(eve)) as f:
        LS=list(map(lambda x:x.strip().split(),f))


    Ntcga=int(os.popen("wc -l %s/Events/%s.tsv"%(os.path.dirname(dataset1),eve)).readline().strip().split()[0])-1
    Nseo=int(os.popen("wc -l %s/Events/%s.tsv"%(os.path.dirname(dataset2),eve)).readline().strip().split()[0])-1

    shared=0
    out=open('Report_%s.txt'%(eve),'w')
    for i in LT:
        for j in LS:
            if i[1]==j[1] and abs(int(i[2])-int(j[2]))<2:# and i[3]==j[3]:
                shared+=1
                out.write(i[0]+'\t'+j[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[5]+'\n')
    print (eve,Ntcga,str(round(100*shared/(len(LS)+len(LT)-shared+1),2))+'%',Nseo)
    os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,Ntcga,str(round(100*shared/(len(LS)+len(LT)-shared+1),2))+'%',Nseo))
    shared=0
    for i in LS:
        for j in LT:
            if i[1]==j[1] and abs(int(i[2])-int(j[2]))<2:shared+=1# and i[3]==j[3]:shared+=1
    shared=0
    out.close()
    #combineTooutput(eve)

def run_split(eve):
    hash={}
    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset2),eve)) as f:
        SEO=list(map(lambda x:x.strip('\n').split('\t'),f))
    pos=SEO[0].index('tag')
    fa=open('%s_dataset2.fa'%eve,'w')
    for i in SEO[1:]:
        fa.write('>%s\n'%i[pos])
        fa.write(i[pos+1]+'\n')
    fa.close()
    if os.path.exists('%s_dataset2'%eve) is False:
        os.system("blastn -query %s_dataset2.fa -db %s -num_threads 4 -task blastn -max_target_seqs 10 -perc_identity 95 -outfmt 6 -no_greedy -max_hsps 4 -evalue 1e-4 -out %s_dataset2 2>>/dev/null"%(eve,genomedb,eve))


    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset1),eve)) as f:
        LUAD=list(map(lambda x:x.strip('\n').split('\t'),f))
    pos=LUAD[0].index('tag')
    fa=open('%s_dataset1.fa'%eve,'w')
    for i in LUAD[1:]:
        fa.write('>%s\n'%i[pos])
        fa.write(i[pos+1]+'\n')
    fa.close()
    if os.path.exists('%s_dataset1'%eve) is False:
        os.system("blastn -query %s_dataset1.fa -db %s -num_threads 4 -task blastn -max_target_seqs 10 -perc_identity 95 -outfmt 6 -no_greedy -max_hsps 4 -evalue 1e-4 -out %s_dataset1 2>>/dev/null"%(eve,genomedb,eve))
    
    if os.path.exists('%s_dataset1'%eve) is False or os.path.exists('%s_dataset2'%eve)is False:
        os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,len(LUAD)-1,'0%',len(SEO)-1))
        return True
    with open('%s_dataset2'%eve) as f:
        LS=list(map(lambda x:x.strip().split(),f))
    with open('%s_dataset1'%eve) as f:
        LT=list(map(lambda x:x.strip().split(),f))
   # LS=[i for i in LS if float(i[2])==100] 
   # LT=[i for i in LT if float(i[2])==100]
   # PT=[i for i in PT if float(i[2])==100]

    luadd=dict(zip([i[-7] for i in LUAD],[i[-6] for i in LUAD]))
    srad=dict(zip([i[-7] for i in SEO],[i[-6] for i in SEO]))
    def creatdict(dat):
        hash={}
        for i in dat:
            if i[0] not in hash.keys():
                hash[i[0]]=[[int(i[6]),int(i[7]),(int(i[8])+int(i[9]))/2,i[0]]]
            elif len(hash[i[0]])==2:continue
            else:
                overlap=0
                for rec in hash[i[0]]:
                    newoverlap=int(i[7])-int(i[6])+rec[1]-rec[0]-(max(int(i[7]),rec[1])-min(int(i[6]),rec[0]))
                    if newoverlap>overlap:overlap=newoverlap
                if overlap<15:
                    hash[i[0]].append([int(i[6]),int(i[7]),(int(i[8])+int(i[9]))/2,i[0]])
        return (hash)
    LS_ha=creatdict(LS)
    LT_ha=creatdict(LT)
    #print (len(LT_ha),len(LS_ha),len(PT_ha))
    #out=open('shared_%s.txt'%eve,'w')
    shared_TCGA2SEO=0
    shared_SEO2TCGA=0
    out2=open('Report_%s.txt'%eve,'w')
    for lt in LT_ha.values():
        for ls in LS_ha.values():
            if len(lt)==len(ls):
                vlt=np.array([i[2] for i in lt])
                vlt.sort()
                vls=np.array([i[2] for i in ls])
                vls.sort()
                if max(abs(vlt-vls))<30:
                    shared_TCGA2SEO+=1
                    CHR=os.popen('grep %s %s/Events/%s.tsv|cut -f 3'%(lt[0][-1],os.path.dirname(dataset1),eve)).readline().strip()
     #               out.write(lt[0][-1]+'\t'+ls[0][-1]+'\t'+'\t'.join(vlt.astype(str))+'\t'+'\t'.join(vls.astype(str))+'\t'+'dataset1_vs_dataset2'+'\t'+CHR+'\n')
                    out2.write(luadd[lt[0][-1]]+'\t'+srad[ls[0][-1]]+'\n')
                    break
    out2.close()
    for ls in LS_ha.values():
        for lt in LT_ha.values():
            if len(lt)==len(ls):
                vlt=np.array([i[2] for i in lt])
                vlt.sort()
                vls=np.array([i[2] for i in ls])
                vls.sort()
                if max(abs(vlt-vls))<30:
                    shared_SEO2TCGA+=1
                    CHR=os.popen('grep %s %s/Events/%s.tsv|cut -f 3'%(ls[0][-1],os.path.dirname(dataset2),eve)).readline().strip()
      #              out.write(ls[0][-1]+'\t'+lt[0][-1]+'\t'+'\t'.join(vls.astype(str))+'\t'+'\t'.join(vlt.astype(str))+'\t'+'dataset2_vs_dataset1'+'\t'+CHR+'\n')
                    break

    print (eve,len(LUAD)-1,str(round(100*max(shared_TCGA2SEO,shared_SEO2TCGA)/(len(LUAD)+len(SEO)-2-max(shared_TCGA2SEO,shared_SEO2TCGA)+1),2))+'%',len(SEO)-1)
    os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,len(LUAD)-1,str(round(100*max(shared_TCGA2SEO,shared_SEO2TCGA)/(len(LUAD)+len(SEO)-2-max(shared_TCGA2SEO,shared_SEO2TCGA)+1),2))+'%',len(SEO)-1))
    #combineTooutput(eve)
    return (shared_TCGA2SEO,shared_SEO2TCGA)

def run_unmap_repeat(eve):
    dataset='dataset2'
    #create fa files
    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset1),eve))as f:
        event=list(map(lambda x:x.strip('\n').split('\t'),f))
    qdict={}
    q_total=len(event)-1
    pos=event[0].index('tag')
    os.system('rm -f %s_q*.fa*'%eve)
    for i in event[1:]:
        qdict[i[pos]]=i[pos+1]
        if len(i[pos+1])<100:
            os.system("echo '>%s' >> %s_q_short.fa"%(i[pos],eve))
            os.system("echo %s >> %s_q_short.fa"%(i[pos+1],eve))
        else:
            os.system("echo '>%s' >> %s_q_long.fa"%(i[pos],eve))
            os.system("echo %s >> %s_q_long.fa"%(i[pos+1],eve))
    

    with open('%s/Events/%s.tsv'%(os.path.dirname(dataset2),eve))as f:
        event=list(map(lambda x:x.strip('\n').split('\t'),f))
    sdict={}
    s_total=len(event)-1
    pos=event[0].index('tag')
    os.system('rm -f %s_s*.fa*'%eve)
    for i in event[1:]:
        sdict[i[pos]]=i[pos+1]
        if len(i[pos+1])<100:
            os.system("echo '>%s' >> %s_s_short.fa"%(i[pos],eve))
            os.system("echo %s >> %s_s_short.fa"%(i[pos+1],eve))
        else:
            os.system("echo '>%s' >> %s_s_long.fa"%(i[pos],eve))
            os.system("echo %s >> %s_s_long.fa"%(i[pos+1],eve))

    if q_total==0 or s_total==0:
        print (eve,q_total,0,s_total)
        os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,q_total,'0%',s_total))
        return True
    
    os.system('cat %s_s_*.fa >> %s_s.fa 2>>/dev/null'%(eve,eve))
    os.system('cat %s_q_*.fa >> %s_q.fa 2>>/dev/null'%(eve,eve))
    #blastn q to s
    out=open('report_q2s_%s'%eve,'w')
    try:
        os.system("blastn -query %s_q_short.fa -subject %s_s.fa -dust no -task blastn -max_target_seqs 4 -perc_identity 99 -outfmt 6 -out report_tmp_q2s_short_%s 2>/dev/null"%(eve,eve,eve))
        if os.path.exists('report_tmp_q2s_short_%s'%eve) is True:
            cutoff_short=os.popen("less report_tmp_q2s_short_%s|cut -f 11"%eve).readlines()
            cutoff_short=[float(i.strip()) for i in cutoff_short]
            cutoff_short.sort()
            ev1qs=max(0.001,min(1,cutoff_short[int(len(cutoff_short)*0.75)]))
            with open('report_tmp_q2s_short_%s'%eve)as f:
                con=list(map(lambda x:x.strip().split(),f))
            for i in con:
                if int(i[3])>30 and float(i[10])<ev1qs:
                    out.write('\t'.join(i)+'\n')
    #    print ('ev q2s cutoff for short is ',ev1qs)
    except:pass

    try:
        os.system("blastn -query %s_q_long.fa -subject %s_s.fa -dust no -task blastn -max_target_seqs 4 -perc_identity 99 -outfmt 6 -best_hit_score_edge 0.1 -best_hit_overhang 0.1 -out report_tmp_q2s_long_%s 2>/dev/null"%(eve,eve,eve))
        if os.path.exists('report_tmp_q2s_long_%s'%eve) is True:
            cutoff_long=os.popen("less report_tmp_q2s_long_%s|cut -f 11"%eve).readlines()
            cutoff_long=[float(i.strip()) for i in cutoff_long]
            cutoff_long.sort()
            ev2qs=max(0.001,min(1,cutoff_long[int(len(cutoff_long)*0.75)]))
            with open('report_tmp_q2s_long_%s'%eve)as f:
                con=list(map(lambda x:x.strip().split(),f))
            for i in con:
                if int(i[3])>30 and float(i[10])<ev2qs:
                    out.write('\t'.join(i)+'\n')
    #    print ('ev q2s cutoff for long is ',ev2qs)
    except:pass
    out.close()
    #blastn s to q
    out=open('report_s2q_%s'%eve,'w')
    try:
        os.system("blastn -query %s_s_short.fa -subject %s_q.fa -dust no -task blastn -max_target_seqs 4 -perc_identity 99 -outfmt 6 -best_hit_score_edge 0.1 -best_hit_overhang 0.1 -out report_tmp_s2q_short_%s 2>logs"%(eve,eve,eve))
        if os.path.exists('report_tmp_s2q_short_%s'%eve) is True:
            cutoff_short=os.popen("less report_tmp_s2q_short_%s|cut -f 11"%eve).readlines()
            cutoff_short=[float(i.strip()) for i in cutoff_short]
            cutoff_short.sort()
            ev1sq=max(0.001,min(1,cutoff_short[int(len(cutoff_short)*0.75)]))
            with open('report_tmp_s2q_short_%s'%eve)as f:
                con=list(map(lambda x:x.strip().split(),f))
            for i in con:
                if int(i[3])>30 and float(i[10])<ev1sq:
                    out.write('\t'.join(i)+'\n')
    #    print ('ev s2q cutoff for short is ',ev1sq)
    except:pass

    try:
        os.system("blastn -query %s_s_long.fa -subject %s_q.fa -dust no -task blastn -max_target_seqs 4 -perc_identity 99 -outfmt 6 -best_hit_score_edge 0.1 -best_hit_overhang 0.1 -out report_tmp_s2q_long_%s 2>logs"%(eve,eve,eve))
        if os.path.exists('report_tmp_s2q_long_%s'%eve) is True:
            cutoff_long=os.popen("less report_tmp_s2q_long_%s|cut -f 11"%eve).readlines()
            cutoff_long=[float(i.strip()) for i in cutoff_long]
            cutoff_long.sort()
            ev2sq=max(0.001,min(1,cutoff_long[int(len(cutoff_long)*0.75)]))
            with open('report_tmp_s2q_long_%s'%eve)as f:
                con=list(map(lambda x:x.strip().split(),f))
            for i in con:
                if int(i[3])>30 and float(i[10])<ev2sq:
                    out.write('\t'.join(i)+'\n')
    #    print ('ev s2q cutoff for long is ',ev2sq)
    except:pass
    out.close()
    
    #summarize shared contigs
    q_overlap=set([i.strip() for i in os.popen('less report_q2s_%s|cut -f 1|sort|uniq'%eve).readlines()])
    s_overlap=set([i.strip() for i in os.popen('less report_s2q_%s|cut -f 1|sort|uniq'%eve).readlines()])
    with open('report_q2s_%s'%eve)as f:
        q2s=list(map(lambda x:x.strip().split(),f))
    with open('report_s2q_%s'%eve)as f:
        s2q=list(map(lambda x:x.strip().split(),f))
    out=open('Report_%s.txt'%(eve),'w')
    for i in q2s:
        for j in s2q:
            if i[0]==j[1] and i[1]==j[0]:
                out.write(qdict[i[0]]+'\t'+sdict[i[1]]+'\n')
    out.close()
    os.system('rm -f report_* %s*.fa'%eve)
    luad,sra=int(os.popen('less Report_%s.txt|cut -f 1|wc -l'%(eve)).readline().strip()),int(os.popen('less Report_%s.txt|cut -f 2|wc -l'%(eve)).readline().strip())
    print (eve,q_total,str(round(100*max(luad,sra)/(q_total+s_total-max(luad,sra)+1),2))+'%',s_total)
    os.system("echo %s %s %s %s >> jaccardidx_input.txt"%(eve,q_total,str(round(100*max(luad,sra)/(q_total+s_total-max(luad,sra)+1),2))+'%',s_total))
    #combineTooutput(eve)

def combineTooutput(eve):
    dif1=pd.read_csv(dataset1,header=0,index_col=0,sep='\t')
    dif1.index=dif1.contig
    if os.path.exists('Report_%s.txt'%eve) is False or int(os.popen('wc -l Report_%s.txt'%eve).readline().strip().split()[0])==0:return True
    eve_report=pd.read_csv('Report_%s.txt'%eve,sep='\t',header=None)
    evedf=dif1.ix[eve_report.loc[:,0]]
    evedf.insert(0,'category',[eve]*evedf.shape[0])
    if os.path.exists('shared_contigs_dataset1_DiffContigInfo.tsv'):evedf.to_csv('shared_contigs_dataset1_DiffContigInfo.tsv',mode='a',index=False,sep='\t',header=False)
    else:evedf.to_csv('shared_contigs_dataset1_DiffContigInfo.tsv',mode='a',index=False,sep='\t',header=True)
    dif2=pd.read_csv(dataset2,header=0,index_col=0,sep='\t')
    dif2.index=dif2.contig
    evedf=dif2.ix[eve_report.loc[:,1]]
    evedf.insert(0,'category',[eve]*evedf.shape[0])
    if os.path.exists('shared_contigs_dataset2_DiffContigInfo.tsv'):evedf.to_csv('shared_contigs_dataset2_DiffContigInfo.tsv',mode='a',index=False,sep='\t',header=False)
    else:evedf.to_csv('shared_contigs_dataset2_DiffContigInfo.tsv',mode='a',index=False,sep='\t',header=True)

for eve in ['gene','kmer']:
    run_diff(eve)
for eve in ['unmapped','repeat']:
    run_unmap_repeat(eve)
    combineTooutput(eve)
for eve in ['SNV','SNV_DU']:
    run_snp(eve)
    combineTooutput(eve)
for eve in ['lincRNA','polyA','polyA_DU','intron','intron_DU']:
    run_noncode(eve)
    combineTooutput(eve)
for eve in ['split','splice','splice_DU']:
    run_split(eve)
    combineTooutput(eve)


os.system('python3 create_traindat.py %s %s %s %s'%(dataset1,dataset2,dkplrun1,dkplrun2))
os.system('Rscript GSEAgraph.r %s %s %s %s'%(dataset1,dataset2,dkplrun1,dkplrun2))
os.system('rm -f *.fa *.txt bam2text* dataset* *dataset1 *dataset2 overlap*')
