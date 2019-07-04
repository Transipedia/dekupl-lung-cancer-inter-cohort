args = commandArgs(trailingOnly=TRUE)
Difftable=args[1]
# Parses DE-kupl contig table and generate annotation tables 
# D. Gautheret Dec 3, 2017

# ---- How to run this script? 
# Input files:
# - Fcontig_info: DE-kupl contig table (DiffContigsInfos.tsv)
# - Fevent_rules: describes criteria for event annotation (splice, polyA, unmapped, etc.) 
# - Fbase_functions: R functions for this script 
# Output files (in directory Outdir):
# - loci.tsv: "by locus" summary: retains 1 contig with highest FC for each locus (gene or intergenic region)
# - <events>.tsv: event tables (for events defined in Fevent_rules)

library(dplyr)

#------------ files and parameters ----------

Fbase_functions = "./contig-functions.r"
Fevent_rules = "./event-rules.4.tsv"     # Edit this table to change event rules or create new events
Fcontig_info = Difftable   # contig table produced by DE-kupl
Outdir=sub("DiffContigsInfos.tsv","Events",Difftable)  # directory for event output files

# set here coordinates of count columns
# varies with number of libraries.
# to be improved with a special prefix in column names

# max number of lines in output tables
# will report the maxoutputlines hits with highest log2FC
maxoutputlines = 10000

#-------- end of parameters --------------

if(!file.exists(Fbase_functions)){stop(paste0(Fbase_functions, ": file not found"))}
if(!file.exists(Fevent_rules)){stop(paste0(Fevent_rules, ": file not found"))}
if(!file.exists(Fcontig_info)){stop(paste0(Fcontig_info, ": file not found"))}
if(!file.exists(Outdir)){dir.create(path=Outdir)}
# adds "/" to Outdir if missing 
if ((substr(Outdir, nchar(Outdir), nchar(Outdir)))!="/"){Outdir=paste0(Outdir,"/")}

# read files
source (Fbase_functions)
rules <- get_rules (Fevent_rules)
ctgs <- readContigs(Fcontig_info) 
lastsample=which(names(ctgs)=='is_mapped')-1
counts=9:lastsample
# convert into Dplyer Tibble  
ctgt <- as_tibble(ctgs)

#round counts
ctgt[,colnames(ctgt[,counts])] = round(ctgt[,colnames(ctgt[,counts])],0)

# group by locus & write locus table
loci <- ctgt %>%  
  group_by(gene_id, upstream_gene_id, downstream_gene_id) %>%
  mutate (ncontigs=n()) %>%
  filter(abs(log2FC) == max(abs(log2FC))) %>%
  arrange(-abs(log2FC)) %>%
  select (gene_symbol, upstream_gene_id, downstream_gene_id, ncontigs, size, log2FC, 
          pvalue, du_pvalue, du_stat, meanA, meanB, chromosome, start, end, 
          tag, contig, counts)

if (nrow(loci)>maxoutputlines){
  loci=head(loci,maxoutputlines)
}

locusfile=paste0(Outdir,"/bestloci.tsv")
write.table(loci,file=locusfile,sep="\t",row.names=F, col.names=T, quote=F)

# sort contigs by increasing fold change
ctgt <- ctgt %>% arrange(-abs(log2FC))

for (i in (1:length(rules))){
  r=rules[i]
  eventname=names(rules)[i]
  eventfile=paste0(Outdir,eventname,".tsv")
  cat (paste0(eventname, " :"))
  ev <- ctgt %>%
    filter(eval(parse(text=r))) %>%
    select (gene_symbol, gene_id, chromosome, start, end, size, log2FC, pvalue, du_pvalue, du_stat, meanA, meanB, 
            counts, tag, contig,exonic,intronic,nb_mismatch,nb_hit,nb_snv)
  cat (paste0("  ",nrow(ev), " events.\n"))
  if (nrow(ev)>maxoutputlines){
    ev=head(ev,maxoutputlines)
  }
  write.table(ev,file=eventfile,sep="\t",row.names=F, col.names=T, quote=F)
}

