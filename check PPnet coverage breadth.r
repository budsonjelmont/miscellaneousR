dat = read.table("C:/Users/Judson/Documents/Phosphonet/KEGG_TCR_PPnetAll.tab", sep="\t", 
  colClasses = c("character","character","character","character","character","numeric","character","character","numeric"))
colnames(dat) = c("substrateUid", "substrateName", "psite", "kinase1Uid", "kinase1name", "k1naScore", "kinase2Uid", "kinase2name", "k2adjScore")

dupes = which(duplicated(dat))
dat = dat[-dupes,]

datST = dat[union(which(substr(dat[,3],1,1)=="S"),which(substr(dat[,3],1,1)=="T")),]
datY = dat[substr(dat[,3],1,1)=="Y",]

datST$residue = "Ser/Thr"
datY$residue = "Tyr"

threshList = seq(0,750,by=25)
results = c()
for(thresh in threshList){
	#ST first
	datSTthreshed = datST[which(datST$k1naScore >= thresh),]
	allKinases = unique(datSTthreshed$kinase1name)
	count = 0
	for(kin in allKinases){
		if(sum(datSTthreshed$kinase1name >= 5)){
			count = count + 1
		}
	}
	results = rbind(results,c(thresh, as.numeric(as.character(count)), "Ser/Thr"))
	#now Y
	datYthreshed = datY[which(datY$k1naScore >= thresh),]
	allKinases = unique(datYthreshed$kinase1name)
	count = 0
	for(kin in allKinases){
		if(sum(datYthreshed$kinase1name >= 5)){
			count = count + 1
		}
	}
	results = rbind(results,c(thresh, as.numeric(as.character(count)), "Tyr"))
}

results = as.data.frame(results)
colnames(results) = c("Threshold","TotalKinases","Residues")
results$TotalKinases = as.numeric(as.character(results$TotalKinases))
results$Threshold = as.numeric(as.character(results$Threshold))

# resultsY = results[which(results$Residues == "Tyr"),]
# resultsST = results[which(results$Residues == "Ser/Thr"),]

#make the plots
p = ggplot(results, aes(x=results$Threshold, y=results$TotalKinases, group=results$Residues, colour=results$Residues)) + 
 geom_point() + geom_line() + scale_colour_discrete(name = "Motif type") + guides(fill=FALSE) +
 scale_x_continuous(name="PhosphoNet motif score") + scale_y_continuous(name="Number of kinases with at least 5 predictions")
plot(p)