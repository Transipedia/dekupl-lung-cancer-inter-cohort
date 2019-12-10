import time
print (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),'loading packages')
from Bio import SeqIO
import sys,os,subprocess
import networkx as nx
import pandas as pd
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout
def contig_fa(dataset1,dataset2):
    subprocess.call(r'''less %s|cut -f 3|awk '{print ">all_A_"NR"\n"$1}' > %s/contig_dataset1.fa'''%(dataset1,outdir),shell=True)
    subprocess.call(r'''less %s|cut -f 3|awk '{print ">all_B_"NR"\n"$1}' > %s/contig_dataset2.fa'''%(dataset2,outdir),shell=True)
    fadict={}
    fasta_sequences = SeqIO.parse(open('%s/contig_dataset1.fa'%outdir),'fasta')
    for seq in fasta_sequences:
        fadict[seq.id]=str(seq.seq)

    fasta_sequences = SeqIO.parse(open('%s/contig_dataset2.fa'%outdir),'fasta')
    for seq in fasta_sequences:
        fadict[seq.id]=str(seq.seq)
    return fadict

def func(id):
    return len(fadict[id])

def translate(id):
    return fadict[id]

def find_cliques():
    subprocess.call("./interSeqGraph -k %s -n -A %s -B %s > %s/pairs_contigs.txt 2>/dev/null"%(ksize,'%s/contig_dataset1.fa'%outdir,'%s/contig_dataset2.fa'%outdir,outdir),shell=True)
    dat=pd.read_csv('%s/pairs_contigs.txt'%outdir,header=0,index_col=None,sep='\t')
    dat['A_L']=dat['contig_in_A'].apply(func,1)
    dat['B_L']=dat['contig_in_B'].apply(func,1)
    dat['count']+=30
    dat['weight']=dat['count']/(dat['A_L']+dat['B_L']-dat['count'])
    del dat['A_L'],dat['B_L'],dat['count']
    dat=dat.sort_values(by=['weight'],ascending=False)
    dat.to_csv('%s/graph_input.txt'%outdir,header=False,index=False,sep='\t')
    G = nx.read_edgelist("%s/graph_input.txt"%outdir,delimiter='\t',nodetype=str, data=(('weight',float),))
    comps=nx.find_cliques(G)
    out=open('%s/shared_event.tsv'%outdir,'w')
    out.write('dataset1\tdataset2\n')
    CL1share,CL2share=[],[]
    dist=[]
    CL1=int(os.popen("grep '>' %s/contig_dataset1.fa|wc -l"%(outdir)).readline().strip())
    CL2=int(os.popen("grep '>' %s/contig_dataset2.fa|wc -l"%(outdir)).readline().strip())
    def longest(L):
        length=np.array([len(i) for i in L])
        return L[np.argmax(length)]

    for clq in comps:
        dist.append(len(clq))
        bagA,bagB=[],[]
        for i in clq:
            if '_A_' in i:
                tia=translate(i)
                bagA.append(tia)
                CL1share.append(tia)
            if '_B_' in i:
                tib=translate(i)
                bagB.append(tib)
                CL2share.append(tib)
        if len(bagA)==0 or len(bagB)==0:continue
        a=longest(bagA)
        b=longest(bagB)
        out.write(a+'\t'+b+'\n')
    out.close()
    dist.sort()
    np.savetxt('%s/clique_distribution.txt'%outdir,np.array(dist),fmt='%-4d')
    CL1_share,CL2_share=len(set(CL1share)),len(set(CL2share))
    print ('shared contigs ratio in dataset1:',round(len(set(CL1share))/CL1,2),'\nshared contigs ratio in dataset2:',round(len(set(CL2share))/CL2,2),'\njaccard index is:',min(1,round(max(CL1_share,CL2_share)/(CL1+CL2-max(CL1_share,CL2_share)),2)))
    os.system('rm %s/*.fa %s/*.txt'%(outdir,outdir))
def output():
    #output DiffCountTable
    DCT1=pd.read_csv(dataset1,header=0,index_col=2,sep='\t')
    DCT2=pd.read_csv(dataset2,header=0,index_col=2,sep='\t')
    share=pd.read_csv('%s/shared_event.tsv'%outdir,header=0,index_col=None,sep='\t')
    DCT1.loc[list(set(share.iloc[:,0]))].to_csv('%s/shared_contigs_dataset1.tsv'%outdir,header=True,index=True,sep='\t')
    DCT2.loc[list(set(share.iloc[:,1]))].to_csv('%s/shared_contigs_dataset2.tsv'%outdir,header=True,index=True,sep='\t')

if __name__ == '__main__':
    if len(sys.argv)<5:
        print ('python3 clique_based_interevents.py dataset1 dataset2 outdir kmersize')
        os._exit(0)
    dataset1,dataset2,outdir,ksize=sys.argv[1:]
    os.system('mkdir %s'%outdir)
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),'create contig fasta')
    global fadict
    fadict=contig_fa(dataset1,dataset2)
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),'search for cliques')
    find_cliques()
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),'write to the output')
    output()
