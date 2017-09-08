library(gdata)
library(stringr)
library(pROC)

ppnet = read.table("C:/Users/Judson/Documents/Phosphonet/KEGG_TCR_PPnetAll.tab", sep="\t", 
  colClasses = c("character","character","character","character","character","numeric","character","character","numeric"))
colnames(ppnet) = c("substrateUid", "substrateName", "psite", "kinase1Uid", "kinase1name", "k1naScore", "kinase2Uid", "kinase2name", "k2adjScore")

dupes = which(duplicated(ppnet))
ppnet = ppnet[-dupes,]

ppnet$UP_psite = as.character(paste(ppnet$substrateUid,ppnet$psite,sep="_"))

ppnetST = ppnet[union(which(substr(ppnet[,3],1,1)=="S"),which(substr(ppnet[,3],1,1)=="T")),]
ppnetY = ppnet[substr(ppnet[,3],1,1)=="Y",]

ppnetST$residue = "Ser/Thr"
ppnetY$residue = "Tyr"

uniqueSites = unique(ppnetST$UP_psite)

hprd2010 = read.table("C:/Users/Judson/Documents/Downloaded DBs/HPRD FLAT_FILES_072010/POST_TRANSLATIONAL_MODIFICATIONS.txt", sep="\t", header=FALSE, colClasses = "character")
keggTCRlist = read.table("C:/Users/Judson/Documents/Downloaded DBs/HPRD FLAT_FILES_072010/kegg TCR list.tab", sep="\t", header=FALSE, colClasses = "character")
idMap = read.table("C:/Users/Judson/Documents/Downloaded DBs/HPRD FLAT_FILES_072010/UNIPROT to HPRD map", sep="\t", header=FALSE, colClasses = "character")

colnames(idMap) = c("UNI", "HPRD")

#Remove entries on non-TCR proteins and non-phosphorylation interactions
hprd2010 = hprd2010[!is.na(match(hprd2010[,1],keggTCRlist[,1])),]
hprd2010 = hprd2010[which(hprd2010[,9]=="Phosphorylation"),]

#add colnames
colnames(hprd2010) = c("sub_id", "sub_name", "sub_id_iso", "gene", "site", "residue", 
  "kin_name", "kin_id", "type", "evid", "dunno")

#I don't have HPRD ids for PhosphoNet predictions, and the HPRD-supplied ID mapping file includes 
#deprecated UNIPROT IDs, so the relation between HPRD:UNIPROT IDs is one:many. 

#EDIT: Better idea, translate all UNIPROT IDs to HPRD using the UNIPROT to HPRD map file
ppnet_uids = unique(ppnet[,1])

ppnet$subhprd = idMap$HPRD[match(ppnet$substrateUid,idMap$UNI)]
ppnet$kin1hprd = idMap$HPRD[match(ppnet$kinase1Uid,idMap$UNI)]
ppnet$kin2hprd = idMap$HPRD[match(ppnet$kinase2Uid,idMap$UNI)]


ppnet$hprd_psite = paste(ppnet$subhprd, ppnet$psite, sep="_")

#add HPRD_psite column to hprd data frame
id = hprd2010$sub_id
site = hprd2010$site
res = hprd2010$residue

site = lapply(site, function(s){
	if(grepl(";",s)){
		pos = regexpr(";", s)
		s = substr(s,1,pos-1)
	}
	return(s)
})

hprd2010$hprd_psite = paste(id,"_",res,site,sep="")

#find all matching predictions between PPnet and HPRD for score threshold = 0, and use
#this to determine which isoforms should be discarded from HPRD
thresh = 0
hprd2010$outcome = NA
ppnet_threshed = ppnet[which(ppnet$k1naScore >= thresh),]
for(r in 1:nrow(hprd2010)){
	if(hprd2010$hprd_psite[r] %in% ppnet_threshed$hprd_psite){
		m = which(ppnet_threshed$hprd_psite == hprd2010$hprd_psite[r])
		if(hprd2010$kin_id[r] %in% na.omit(ppnet_threshed$kin1hprd[m])){
			hprd2010$outcome[r] = 1
		} else {
		#If HPRD has an unknown kinase annotated, then give full credit for any match in PPnet
		if(hprd2010$kin_id == "-"){
			hprd2010$outcome[r] = 1
			} else {
				hprd2010$outcome[r] = 0
			}
		}
	} else {
		hprd2010$outcome[r] = -1
	}
}

isoformDeredundinator = list()
for(id in unique(hprd2010$sub_id)){
	r = which(hprd2010$sub_id==id)
	iso = unique(hprd2010$sub_id_iso[r])
	for(i in iso){
		rows = which(hprd2010$sub_id_iso == i)
		outcomes = hprd2010$outcome[rows]
		isoformDeredundinator[[id]] = rbind(isoformDeredundinator[[id]],c(i,length(rows),sum(outcomes==-1)))
	}
}
#iterate through isoformDeredundinator and get the isoform that is best represented in PPnet
isoToKeep = c()
for(n in names(isoformDeredundinator)){
	isoToKeep = c(isoToKeep, isoformDeredundinator[[n]][which.min(as.numeric(isoformDeredundinator[[n]][,3])),1])
}
hprd2010_filt = hprd2010[which(hprd2010$sub_id_iso %in% isoToKeep),]
#Iterate over rows in PPnet table once to determine if each site is predicted in HPRD
ppnet_threshed = ppnet
ppnet_threshed = ppnet_threshed[-which(is.na(ppnet_threshed$kin1hprd)),]
ppnet_threshed$label = NA
ppnet_threshed$sitePredicted = NA
for(i in 1:nrow(ppnet_threshed)){
		#If PPnet site is annotated in HPRD
	if(ppnet_threshed$hprd_psite[i] %in% hprd2010_filt$hprd_psite){
		m = which(hprd2010_filt$hprd_psite == ppnet_threshed$hprd_psite[i])
		#If PPnet site is correctly annotated OR kinase is unknown in HPRD
		if(any(hprd2010_filt$kin_id[m] == "-") || any(hprd2010_filt$kin_id[m] == ppnet_threshed$kin1hprd[i])){
			ppnet_threshed$siteInHprd[i] = 1
		} else {
			ppnet_threshed$siteInHprd[i] = 0
		}
	} else {
		ppnet_threshed$siteInHprd[i] = 0
	}
}
#Iterate over rows in PPnet table and assign each prediction a label (TP, FP, TN, or FN)
results = list()
for(thresh in seq(100,600,100)){
	ppnet_threshed$sitePredicted = NA
	for(i in 1:nrow(ppnet_threshed)){
		if(ppnet_threshed$k1naScore[i] >= thresh){
				ppnet_threshed$sitePredicted[i] = 1
			} else {
				ppnet_threshed$sitePredicted[i] = 0
			}
	}
	results[[as.character(thresh)]] = ppnet_threshed
}
