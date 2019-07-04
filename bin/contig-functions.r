# Fonctions for reading and parsing DE-kupl contig files

# read contig tsv file, and add/reformat some columns
readContigs<-function(Fcontig_info) {
  print ("reading contigs...")
  ctgs=read.table(Fcontig_info, header=T, sep="\t")
  print (paste0 (nrow(ctgs), " contigs read."))
  ctgs$contig = as.character(ctgs$contig) # need contig sequences as chars!
  ctgs <- ctgs[order(-abs(ctgs$log2FC)),]  # sort by FC
  ctgs$size=nchar(ctgs$contig)
  return(ctgs)
}

# reads a table with all event annotation rules
# creates filtering commands
get_rules <- function (fname) {
  annot=read.table(fname, header=T, sep="\t", check.names=F, as.is=T,quote="")
  nbev=nrow(annot)
  allrules=c()
  evnames=c()
  for (i in (1:nbev)){
    eventname=annot[i,1]
    evnames = c(evnames, eventname)
    rule=c()
    for (j in (2:(ncol(annot)-1))) {
      ru1=annot[i,j]
      if (ru1 != "") {
        field=colnames(annot)[j]
        str=paste0(field, ru1)
        rule=c(rule,str)
      }
    }  
    # special case: "other" rule is a full command
    ru1=annot[i,j+1]
    if (ru1 != "") {
      rule=c(rule,ru1)
    }
    com = paste0(rule,sep= " & ",collapse='')
    com = substr(com,1,nchar(com)-3)
    allrules= c(allrules, com)
  }
  names(allrules)=evnames
  return(allrules)
}

# Generates logical test for annotating a specific event 
#   and returns lines meeting test
# evname: name of event
# allrules: annotation rules table
# ctgs: all contigs

build_events <- function (evname, allrules, ctgs) {
  com= paste0("subset(ctgs,",allrules[evname],")",sep= " ",collapse='')
  print (com)
  temp=eval(parse(text=com))
  return(temp)
}
