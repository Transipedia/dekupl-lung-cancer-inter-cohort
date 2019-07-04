from __future__ import division
import sys,re,os
import numpy as np
import pandas as pd
import scipy.stats as stats
import itertools

def Ana_all_possible_snp(allsnp,contig):
    possible_snp=[]
    clean_contig=list(contig)
    try:
        for p_idx in range(len(allsnp)):
            insertion,deletion,insert_nb=0,0,0
            if '^' in ''.join(allsnp[:p_idx+1]):
                deletion=len(re.findall('\^([ATCG]+)',''.join(allsnp[:p_idx+1]))[0])
            elif '^' not in ''.join(allsnp[:p_idx]):
                deletion=0
            if 'I' not in allsnp[0]:
                pos=sum([int(re.search('[0-9]+',allsnp[idx]).group(0)) for idx in range(p_idx+1)])+sum([len(re.search('[ATCG]+',allsnp[idx]).group(0)) for idx in range(p_idx)])+1-deletion
            elif 'I' in allsnp[0]:
                insert_nb=len(re.findall('\d+I([ATCG]+)',allsnp[0])[0])
                insertion=int(re.findall('(\d+)I',allsnp[0])[0])+len(re.findall('\d+I([ATCG]+)',allsnp[0])[0])
                pos=sum([int(re.search('[0-9]+',allsnp[idx]).group(0)) for idx in range(p_idx+1)])+sum([len(re.search('[ATCG]+',allsnp[idx]).group(0)) for idx in range(p_idx)])+1-deletion-insertion
            if '^' in allsnp[p_idx]:
                ref,alt=clean_contig[pos-1]+re.search('[ATCG]+',allsnp[p_idx]).group(0),clean_contig[pos-1]
                clean_contig[pos-1]=ref
            elif 'I' in allsnp[p_idx]:
                pos=int(re.findall('(\d+)I',allsnp[p_idx])[0])
                ref,alt=clean_contig[pos-1],clean_contig[pos-1]+re.findall('\d+I([ATCG]+)',allsnp[p_idx])[0]
                clean_contig=np.delete(clean_contig,range(pos,pos+len(alt)-1)).astype('U8')
            else:
                if pos<insertion:insert_nb=0
                ref,alt=re.search('[ATCG]+',allsnp[p_idx]).group(0),contig[pos-1+insert_nb:pos-1+len(re.search('[ATCG]+',allsnp[p_idx]).group(0))+insert_nb] #in case there is a insertion(>1)
                clean_contig[pos-1:pos-1+len(ref)]=list(ref)
            possible_snp.append((pos,ref,alt))
        possible_snp=np.array(possible_snp)
        return (possible_snp)
    except:
        print (allsnp)
        return (possible_snp)


def parsemutation(file):
    with open(file)as f:
        content=list(map(lambda x:x.strip().split(),f))
    for line in content:
        if 'N' not in line[3] and line[1]!='MT':
            i,j=0,len(line[4])        
            soft_clip=re.findall('^(\d+S)*.*?(\d+S)*$',line[3])
            if 'S' in soft_clip[0][0]:i=int(soft_clip[0][0].strip('S'))
            if 'S' in soft_clip[0][1]:j=j-int(soft_clip[0][1].strip('S'))
            contig=line[4][i:j]
            if 'I' in line[3]:
                ins_pos=re.findall('(\d+)M(\d+)I',line[3])[0]
                ALTS=[ins_pos[0]+'I'+contig[int(ins_pos[0]):int(ins_pos[0])+int(ins_pos[1])]]
                ALTS+=re.findall('(\d+\^*[ATCG]+)',line[5])
                for (pos,ref,alt) in Ana_all_possible_snp(ALTS,contig):
                    os.system("echo '%s %s %s %s > %s %s %s' >> %s_mutinfo "%(line[-6],line[1],int(pos)+int(line[2]),ref,alt,line[3],line[5],file))
            elif 'I' not in line[3]:
                ALTS=re.findall('(\d+\^*[ATCG]+)',line[5])
                for (pos,ref,alt) in Ana_all_possible_snp(ALTS,contig):
                    os.system("echo '%s %s %s %s > %s %s %s' >> %s_mutinfo "%(line[-6],line[1],int(pos)+int(line[2]),ref,alt,line[3],line[5],file))
            else:
                print ('other',line)
if __name__ == '__main__':
    filename=sys.argv[3]
    dataset=sys.argv[2]
    eve=sys.argv[1]
    os.system('rm -f %s_%s_mutinfo'%(filename,eve))
    os.system("samtools view %s/contigs.bam|awk '{if(($15~/NM\:i\:[1-3]$/)&&($13~/NH\:i\:1$/)){print $1,$3,$4,$6,$10,$12}}' > bam2text_%s_%s"%(os.path.dirname(dataset),filename,eve))
    df=pd.read_csv('bam2text_%s_%s'%(filename,eve),header=None,index_col=0,sep=' ')
    snv=pd.read_csv('%s/Events/%s.tsv'%(os.path.dirname(dataset),eve),header=0,index_col=False,sep='\t')
    snv.index=snv['tag']
    comb=pd.merge(df,snv,left_index=True,right_index=True)
    comb.to_csv('%s_%s'%(filename,eve),header=None,index=True,sep='\t')
    parsemutation("%s_%s"%(filename,eve))

