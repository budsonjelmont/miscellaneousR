library(gdata)
library(stringr)

alt = "less" #one of "greater", "less", or "two.sided"

dat = read.xls("C:/Users/Judson/Documents/CSK CD45/CSK_CD45_LCK_pY_UCSF_10mowse_MediumScansiteStringency.xls")
#remove rows missing an IPI accession # (which should actually be UNIPROT)
dat = dat[-which(dat$UNIPROT.accession.number == ""),]
#remove rows where IPI accession is actually an IPI accession
dat = dat[-which(substr(dat$UNIPROT.accession.number,1,3) == "IPI"),]

####remove rows with peptides that have > 1 pp'ated site#####################################
dat = dat[-which(dat$Number.of.phosphorylated.sites > 1),]
#############################################################################################

#####Deredundination#########################################################################
#Expansion
# nr_dat = dat[0,]
# for (i in 1:nrow(dat)){
	# for (j in strsplit(as.character(dat$Return.delimited.phosphosites[i]),"\n")[[1]]){
		# newRow = dat[i,]
		# newRow$phosphosite.annotated = j
		# nr_dat = rbind(nr_dat, newRow)
	# }
# }
# #Contraction
# to.keep = c()
# for(i in 1:nrow(nr_dat)){
	# dupes = which(nr_dat$UNIPROT.accession.number == nr_dat$UNIPROT.accession.number[i] & nr_dat$phosphosite.annotated == nr_dat$phosphosite.annotated[i])
	# temp = nr_dat[dupes,]
	# to.keep = c(to.keep, dupes[which.max(median(c(temp$peakarea.manual.3.timepoint1,temp$peakarea.manual.3.timepoint2,temp$peakarea.manual.3.timepoint3,
	# temp$peakarea.manual.3.timepoint4,temp$peakarea.manual.3.timepoint5,temp$peakarea.manual.3.timepoint6,
	# temp$peakarea.manual.3.timepoint7,temp$peakarea.manual.3.timepoint8,temp$peakarea.manual.3.timepoint9,
	# temp$peakarea.manual.3.timepoint10),na.rm=TRUE))])
# }
# #Cleaning up
# nr_dat = nr_dat[to.keep,]
# nr_dat = nr_dat[which(substr(nr_dat[,4],1,1)=="Y"),]
# dat = nr_dat
#############################################################################################

#make unique list of all scansite annotation in the dataset
ss_motifs = apply(dat,1,function(x){print(strsplit(as.character(x[48]),"\n")[[1]])})
ssMotifsInData = unlist(ss_motifs)
ssMotifsInData = unique(ssMotifsInData)

#make list keyed on SS motif containing vector of all associated phosphosites in this dataset
gmt = list()
for(i in ssMotifsInData){
	for(j in 1:nrow(dat)){
		if(i %in% ss_motifs[[j]]){
			gmt[[i]] = c(gmt[[i]], paste(dat$UNIPROT.accession.number[j],gsub("\\*","",dat$phosphosite.annotated[j]),sep="_"))
		}
	}
}

#make list of all sites that are significantly up/down-regulated per dataset
significant = list()
for(i in 1:nrow(dat)){
	#CSK
	qvals = na.omit(c(dat$qvalues.for.SILAC.timepoint1[i], dat$qvalues.for.SILAC.timepoint2[i], dat$qvalues.for.SILAC.timepoint3[i]))
	if(any(qvals) < .05){
		ss = which(qvals < .05)
		ratios = na.omit(c(dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint1[i], dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint2[i],dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint3[i]))[ss] 
		if(all(ratios > 1) || all(ratios < 1)){
			print("CSK")
			print(ratios)
			significant[["CSK"]] = c(significant[["CSK"]], paste(dat$UNIPROT.accession.number[i],gsub("\\*","",dat$phosphosite.annotated[i]),sep="_"))
		}
	}
	#CD45
	qvals = na.omit(c(dat$qvalues.for.SILAC.timepoint4[i], dat$qvalues.for.SILAC.timepoint5[i], dat$qvalues.for.SILAC.timepoint6[i]))
	if(any(qvals < .05)){
		ss = which(qvals < .05)
		ratios = na.omit( c(dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint4[i], dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint5[i],dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint6[i]))[ss]
		if(all(ratios > 1) || all(ratios < 1)){
			significant[["CD45"]] = c(significant[["CD45"]], paste(dat$UNIPROT.accession.number[i],gsub("\\*","",dat$phosphosite.annotated[i]),sep="_"))
		}
	}
	#CSKCD45
	qvals = na.omit(c(dat$qvalues.for.SILAC.timepoint7[i], dat$qvalues.for.SILAC.timepoint8[i], dat$qvalues.for.SILAC.timepoint9[i]))
	if(any(qvals < .05)){
		ss = which(qvals < .05)
		ratios = na.omit(c(dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint7[i], dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint8[i],dat$SILAC.ratio.23.for.user.selected.SILAC.timepoint9[i]))[ss]
		if(all(ratios > 1) || all(ratios < 1)){
			significant[["CSKCD45"]] = c(significant[["CSKCD45"]], paste(dat$UNIPROT.accession.number[i],gsub("\\*","",dat$phosphosite.annotated[i]),sep="_"))
		}
	}
}
csk_sig = significant[["CSK"]]
cd45_sig = significant[["CD45"]]
cskcd45_sig = significant[["CSKCD45"]]

#set up table that will hold enrichment results
results = data.frame(Motif=ssMotifsInData,
	CSKpval=as.numeric(rep(NA,length(ssMotifsInData))),
	CD45pval=as.numeric(rep(NA,length(ssMotifsInData))),
	CSKCD45pval=as.numeric(rep(NA,length(ssMotifsInData))),
	CSKqval=as.numeric(rep(NA,length(ssMotifsInData))),
	CD45qval=as.numeric(rep(NA,length(ssMotifsInData))),
	CSKCD45qval=as.numeric(rep(NA,length(ssMotifsInData))),
	CSKratioEnriched=as.character(rep(NA,length(ssMotifsInData))),
	CD45ratioEnriched=as.character(rep(NA,length(ssMotifsInData))),
	CSKCD45ratioEnriched=as.character(rep(NA,length(ssMotifsInData))),
	CSKssSitesWithMotif=as.character(rep(NA,length(ssMotifsInData))),
	CD45ssSitesWithMotif=as.character(rep(NA,length(ssMotifsInData))),
	CSKCD45ssSitesWithMotif=as.character(rep(NA,length(ssMotifsInData))),
	stringsAsFactors=FALSE)
rownames(results) = ssMotifsInData
#iterate through list of ss motifs in data, and test for enrichment with Fisher's exact test
CSKresults = list()
CSKtables = list()
CD45results = list()
CD45tables = list()
CSKCD45results = list()
CSKCD45tables = list()
for (i in ssMotifsInData){
	#CSK
	AllWithMotif = length(gmt[[i]])
	SigWithMotif = length(intersect(gmt[[i]],significant[["CSK"]]))
	InsigWithMotif = AllWithMotif - SigWithMotif
	WithoutMotif = nrow(dat) - AllWithMotif
	SigWithoutMotif = length(significant[["CSK"]]) - SigWithMotif
	InsigWithoutMotif = WithoutMotif - length(significant[["CSK"]])
	
	cont.table = matrix(c(SigWithMotif, SigWithoutMotif, InsigWithMotif, InsigWithoutMotif),nrow=2,byrow=TRUE)
	test = fisher.test(cont.table, alternative = alt)
	# if(test$p.value < 0.05){
		# print(paste(i,"CSK:",sep=","))
		# print(cont.table)
		# print(test)
		CSKresults[[i]] = test
		CSKtables[[i]] = cont.table
		results[i,"CSKpval"] = test$p.value
		results[i,"CSKratioEnriched"] = paste(SigWithMotif,SigWithMotif+InsigWithMotif,sep="/")
		results[i,"CSKssSitesWithMotif"] = paste(intersect(gmt[[i]],significant[["CSK"]]),collapse="/")
	# }
	#CD45
	AllWithMotif = length(gmt[[i]])
	SigWithMotif = length(intersect(gmt[[i]],significant[["CD45"]]))
	InsigWithMotif = AllWithMotif - SigWithMotif
	WithoutMotif = nrow(dat) - AllWithMotif
	SigWithoutMotif = length(significant[["CD45"]]) - SigWithMotif
	InsigWithoutMotif = WithoutMotif - length(significant[["CD45"]])
	
	cont.table = matrix(c(SigWithMotif, SigWithoutMotif, InsigWithMotif, InsigWithoutMotif),nrow=2,byrow=TRUE)
	test = fisher.test(cont.table, alternative=alt)
	# if(test$p.value < 0.05){
		# print(paste(i,"CD45:",sep=","))
		# print(cont.table)
		# print(test)
		CD45results[[i]] = test
		CD45tables[[i]] = cont.table
		results[i,"CD45pval"] = test$p.value
		results[i,"CD45ratioEnriched"] = paste(SigWithMotif,SigWithMotif+InsigWithMotif,sep="/")
		results[i,"CD45ssSitesWithMotif"] = paste(intersect(gmt[[i]],significant[["CD45"]]),collapse="/")
	# }
	#CSK/CD45
	AllWithMotif = length(gmt[[i]])
	SigWithMotif = length(intersect(gmt[[i]],significant[["CSKCD45"]]))
	InsigWithMotif = AllWithMotif - SigWithMotif
	WithoutMotif = nrow(dat) - AllWithMotif
	SigWithoutMotif = length(significant[["CSKCD45"]]) - SigWithMotif
	InsigWithoutMotif = WithoutMotif - length(significant[["CSKCD45"]])
	
	cont.table = matrix(c(SigWithMotif, SigWithoutMotif, InsigWithMotif, InsigWithoutMotif),nrow=2,byrow=TRUE)
	test = fisher.test(cont.table, alternative=alt)
	# if(test$p.value < 0.05){
		# print(paste(i,"CSKCD45:",sep=","))
		# print(cont.table)Motif
		# print(test)
		CSKCD45results[[i]] = test
		CSKCD45tables[[i]] = cont.table
		results[i,"CSKCD45pval"] = test$p.value
		results[i,"CSKCD45pval"] = test$p.value
		results[i,"CSKCD45ratioEnriched"] = paste(SigWithMotif,SigWithMotif+InsigWithMotif,sep="/")
		results[i,"CSKCD45ssSitesWithMotif"] = paste(intersect(gmt[[i]],significant[["CSKCD45"]]),collapse="/")
	# }
}

results$CSKpval[which(results$CSKpval >= 1)]=0.999
results$CD45pval[which(results$CD45pval >= 1)]=0.999
results$CSKCD45pval[which(results$CSKCD45pval >= 1)]=0.999

##### Control the FDR #####
library(qvalue)
results$CSKqval = qvalue(results$CSKpval, fdr.level=0.01)$qvalues
results$CD45qval = qvalue(results$CD45pval, fdr.level=0.01)$qvalues
results$CSKCD45qval = qvalue(results$CSKCD45pval, fdr.level=0.01)$qvalues

write.table(results, "", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)