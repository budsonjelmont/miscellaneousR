library(ggplot2)
library(gtable)
library(grid)

dat = read.table("C:/Users/Judson/Documents/Phosphonet/KEGG_TCR_PPnetAll.tab", sep="\t", 
  colClasses = c("character","character","character","character","character","numeric","character","character","numeric"))
colnames(dat) = c("substrateUid", "substrateName", "psite", "kinase1Uid", "kinase1name", "k1naScore", "kinase2Uid", "kinase2name", "k2adjScore")

dupes = which(duplicated(dat))
dat = dat[-dupes,]

dat$UP_psite = as.character(paste(dat$substrateUid,dat$psite,sep="_"))

datST = dat[union(which(substr(dat[,3],1,1)=="S"),which(substr(dat[,3],1,1)=="T")),]
datY = dat[substr(dat[,3],1,1)=="Y",]

datST$residue = "Ser/Thr"
datY$residue = "Tyr"

uniqueSites = unique(datST$UP_psite)

for(s in uniqueSites){
	scores = datST$k1naScore[which(datST$UP_psite==s)]
	p = qplot(scores, geom="histogram")
	ggsave(paste("SerThr/",s, " score distribution for 50 kinase predictions.png",sep=""), p)
}

uniqueSites = unique(datY$UP_psite)

for(s in uniqueSites){
	scores = datY$k1naScore[which(datY$UP_psite==s)]
	p = qplot(scores, geom="histogram")
	ggsave(paste("Tyr/",s, " score distribution for 50 kinase predictions.png",sep=""), p)
	#Make QQ plot
	jpeg()
	qqnorm(scores)
	dev.off()
}

####################################################################################################################################
eruption.lm = lm(eruptions ~ k1naScore, data=scores) 
eruption.stdres = rstandard(eruption.lm)

#We now create the normal probability plot with the qqnorm function, and add the qqline for further comparison.
qqnorm(eruption.stdres, 
     ylab="Standardized Residuals", 
     xlab="Normal Scores", 
     main="Old Faithful Eruptions") 
 qqline(eruption.stdres) 
####################################################################################################################################
skew <- function(iList){
    iUnList <- unlist(iList)
    m <- mean(iUnList, na.rm=TRUE)
    N <- length(iUnList)
    #for both kurtosis and skew, SD is calculated with N in denominator instead of N-1
    stdev <- sqrt(sum((iUnList - mean(iUnList))^2) / (N))
    skewness <- ((sum(iUnList - m)^3)/N)/(stdev^3)
    return(skewness)
}

dat$percentile = NA
dat$stDev = NA
dat$range = NA
dat$skew = NA

for(r in 1:nrow(dat)){
	datSite = dat[which(dat$UP_psite == dat$UP_psite[r]),]
	datSite = datSite[order(-datSite$k1naScore),]
	percentGT = ecdf(datSite$k1naScore)(dat$k1naScore[r])
	stdev = sd(datSite$k1naScore)
	dat$percentile[r] = percentGT
	dat$stDev[r] = stdev
	dat$range[r] = datSite$k1naScore[1] - datSite$k1naScore[length(datSite$k1naScore)]
	dat$skew[r] = skew(datSite$k1naScore)
}

dat$msXpXsd = dat$k1naScore * dat$percentile * dat$stDev
dat$msXp = dat$k1naScore * dat$percentile
dat$scoreLogNorm = dat$k1naScore * (dat$percentile + log(dat$stDev/mean(dat$stDev)))
dat$scoreLog2Norm = dat$k1naScore * (dat$percentile + log2(dat$stDev/mean(dat$stDev)))

ppnet = dat

####################################################################################################
####################################################################################################
#Check PhosphoNet against 2010 download of HPRD database############################################
####################################################################################################
####################################################################################################

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

#Iterate over score thresholds, get matching rows in ppnet table, and get the total # of hits from ppnet

ST_anyMatchScore_series = c()
Y_anyMatchScore_series = c()
ST_ExactMatchScore_series = c()
Y_ExactMatchScore_series = c()
#Now iterate through different thresholds and record yields
for(thresh in seq(0,1250,25)){
	hprd2010_filt$outcome = NA
	ppnet_threshed = ppnet[which(ppnet$scoreLogNorm >= thresh),]
	for(r in 1:nrow(hprd2010_filt)){
		if(hprd2010_filt$hprd_psite[r] %in% ppnet_threshed$hprd_psite){
			m = which(ppnet_threshed$hprd_psite == hprd2010_filt$hprd_psite[r])
			if(hprd2010_filt$kin_id[r] %in% na.omit(ppnet_threshed$kin1hprd[m])){
				hprd2010_filt$outcome[r] = 1
			} else {
			#If HPRD has an unknown kinase annotated, then give full credit for any match in PPnet
			if(hprd2010_filt$kin_id[r] == "-"){
				hprd2010_filt$outcome[r] = 1
				} else {
					hprd2010_filt$outcome[r] = 0
				}
			}
		} else {
			hprd2010_filt$outcome[r] = -1
		}
	}

	hprd2010_filt_ST = hprd2010_filt[union(which(hprd2010_filt$residue == "S"),which(hprd2010_filt$residue == "T")),]
	hprd2010_filt_Y = hprd2010_filt[which(hprd2010_filt$residue == "Y"),]

	ST_proportionAnyMatch = (sum(hprd2010_filt_ST$outcome == 0) + sum(hprd2010_filt_ST$outcome == 1))/nrow(hprd2010_filt_ST)
	Y_proportionAnyMatch = (sum(hprd2010_filt_Y$outcome == 0) + sum(hprd2010_filt_Y$outcome == 1))/nrow(hprd2010_filt_Y)
	ST_anyMatchScore_series = rbind(ST_anyMatchScore_series, c(thresh, ST_proportionAnyMatch,"Ser/Thr"))
	Y_anyMatchScore_series = rbind(Y_anyMatchScore_series, c(thresh, Y_proportionAnyMatch,"Tyr"))
	
	ST_proportionExactMatch = sum(hprd2010_filt_ST$outcome == 1)/nrow(hprd2010_filt_ST)
	Y_proportionExactMatch = sum(hprd2010_filt_Y$outcome == 1)/nrow(hprd2010_filt_Y)
	ST_ExactMatchScore_series = rbind(ST_ExactMatchScore_series, c(thresh, ST_proportionExactMatch,"Ser/Thr"))
	Y_ExactMatchScore_series = rbind(Y_ExactMatchScore_series, c(thresh, Y_proportionExactMatch,"Tyr"))
}

# plot(Y_ExactMatchScore_series,type="n",axes=F,col=3, ylab = "Proportion of HPRD 2010-annotated sites correctly predicted in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(Y_ExactMatchScore_series,col='green',lwd=2.5) 
# par(new=T)
# plot(ST_ExactMatchScore_series,type='n',col=2, ylab = "Proportion of HPRD 2010-annotated sites correctly predicted in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(ST_ExactMatchScore_series,col='blue',lwd=2.5) 
# par(new=F)
# legend('topright', # places a legend at the appropriate place
# c('Tyr','Ser/Thr'), # puts text in the legend
# lty=c(1,1), # gives the legend appropriate symbols (lines)
# lwd=c(2.5,2.5),col=c('green','blue'))

# plot(Y_anyMatchScore_series,type="n",axes=F,col=3, ylab = "Proportion of HPRD 2010-annotated sites with any prediction in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(Y_anyMatchScore_series,col='green',lwd=2.5) 
# par(new=T)
# plot(ST_anyMatchScore_series,type='n',col=2, ylab = "Proportion of HPRD 2010-annotated sites with any prediction in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(ST_anyMatchScore_series,col='blue',lwd=2.5) 
# par(new=F)
# legend('topright', # places a legend at the appropriate place
# c('Tyr','Ser/Thr'), # puts text in the legend
# lty=c(1,1), # gives the legend appropriate symbols (lines)
# lwd=c(2.5,2.5),col=c('green','blue'))

ST_anyMatchScore_series = data.frame(ST_anyMatchScore_series)
colnames(ST_anyMatchScore_series) = c("thresh","proportion","motif_type")
Y_anyMatchScore_series = data.frame(Y_anyMatchScore_series)
colnames(Y_anyMatchScore_series) = c("thresh","proportion","motif_type")
ST_ExactMatchScore_series = data.frame(ST_ExactMatchScore_series)
colnames(ST_ExactMatchScore_series) = c("thresh","proportion","motif_type")
Y_ExactMatchScore_series = data.frame(Y_ExactMatchScore_series)
colnames(Y_ExactMatchScore_series) = c("thresh", "proportion","motif_type")

ExactMatchScore_series = rbind(ST_ExactMatchScore_series, Y_ExactMatchScore_series)
ExactMatchScore_series = data.frame(ExactMatchScore_series)
colnames(ExactMatchScore_series) = c("thresh", "proportion", "motif_type")
ExactMatchScore_series[,1] = as.numeric(as.character(ExactMatchScore_series[,1]))
ExactMatchScore_series[,2] = as.numeric(as.character(ExactMatchScore_series[,2]))


ppnet$residue = NA 
ppnet$residue[union(which(substr(ppnet[,3],1,1)=="S"),which(substr(ppnet[,3],1,1)=="T"))] = "Ser/Thr"
ppnet$residue[which(substr(ppnet[,3],1,1)=="Y")] = "Tyr"


grid.newpage()

p1 = ggplot(ExactMatchScore_series, aes(x = thresh, y = proportion, group = motif_type,
  colour = motif_type)) + geom_line() + geom_point() + theme_classic()
p2 = ggplot(ppnet, aes(x=ppnet$k1naScore, group=ppnet$residue, colour=ppnet$residue, fill=ppnet$residue)) + 
  geom_density(alpha=0.2, position="identity") + theme_bw() %+replace% 
    theme(panel.background = element_rect(fill = NA), axis.ticks.y = element_blank())

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
    pp$l, pp$b, pp$l)

# axis tweaks
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
#ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
ax$grobs[[1]]$x <- ax$grobs[[1]]$x + unit(1, "npc") - unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

# draw it
grid.draw(g)

#Any match in PhosphoNet:

AnyMatchScore_series = rbind(ST_anyMatchScore_series, Y_anyMatchScore_series)
AnyMatchScore_series = data.frame(AnyMatchScore_series)
colnames(AnyMatchScore_series) = c("thresh", "proportion", "motif_type")
AnyMatchScore_series[,1] = as.numeric(as.character(AnyMatchScore_series[,1]))
AnyMatchScore_series[,2] = as.numeric(as.character(AnyMatchScore_series[,2]))

grid.newpage()

p1 = ggplot(AnyMatchScore_series, aes(x = thresh, y = proportion, group = motif_type,
  colour = motif_type)) + geom_line() + geom_point() + theme_classic()
p2 = ggplot(ppnet, aes(x=ppnet$k1naScore, group=ppnet$residue, colour=ppnet$residue, fill=ppnet$residue)) + 
  geom_density(alpha=0.2, position="identity") + theme_bw() %+replace% 
    theme(panel.background = element_rect(fill = NA), axis.ticks.y = element_blank())

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
    pp$l, pp$b, pp$l)

# axis tweaks
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
#ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
ax$grobs[[1]]$x <- ax$grobs[[1]]$x + unit(1, "npc") - unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

# draw it
grid.draw(g)

############################################################################################################################################################################
#Keep all hits at or above X percentile#####################################################################################################################################
############################################################################################################################################################################

#Iterate over score thresholds, get matching rows in ppnet table, and get the total # of hits from ppnet
ST_anyMatchScore_series = c()
Y_anyMatchScore_series = c()
ST_ExactMatchScore_series = c()
Y_ExactMatchScore_series = c()
#Now iterate through different thresholds and record yields
for(thresh in seq(0,1,0.1)){
	hprd2010_filt$outcome = NA
	ppnet_threshed = ppnet[which(ppnet$percentile >= thresh),]
	for(r in 1:nrow(hprd2010_filt)){
		if(hprd2010_filt$hprd_psite[r] %in% ppnet_threshed$hprd_psite){
			m = which(ppnet_threshed$hprd_psite == hprd2010_filt$hprd_psite[r])
			if(hprd2010_filt$kin_id[r] %in% na.omit(ppnet_threshed$kin1hprd[m])){
				hprd2010_filt$outcome[r] = 1
			} else {
			#If HPRD has an unknown kinase annotated, then give full credit for any match in PPnet
			if(hprd2010_filt$kin_id[r] == "-"){
				hprd2010_filt$outcome[r] = 1
				} else {
					hprd2010_filt$outcome[r] = 0
				}
			}
		} else {
			hprd2010_filt$outcome[r] = -1
		}
	}

	hprd2010_filt_ST = hprd2010_filt[union(which(hprd2010_filt$residue == "S"),which(hprd2010_filt$residue == "T")),]
	hprd2010_filt_Y = hprd2010_filt[which(hprd2010_filt$residue == "Y"),]

	ST_proportionAnyMatch = (sum(hprd2010_filt_ST$outcome == 0) + sum(hprd2010_filt_ST$outcome == 1))/nrow(hprd2010_filt_ST)
	Y_proportionAnyMatch = (sum(hprd2010_filt_Y$outcome == 0) + sum(hprd2010_filt_Y$outcome == 1))/nrow(hprd2010_filt_Y)
	ST_anyMatchScore_series = rbind(ST_anyMatchScore_series, c(thresh, ST_proportionAnyMatch,"Ser/Thr"))
	Y_anyMatchScore_series = rbind(Y_anyMatchScore_series, c(thresh, Y_proportionAnyMatch,"Tyr"))

	ST_proportionExactMatch = sum(hprd2010_filt_ST$outcome == 1)/nrow(hprd2010_filt_ST)
	Y_proportionExactMatch = sum(hprd2010_filt_Y$outcome == 1)/nrow(hprd2010_filt_Y)
	ST_ExactMatchScore_series = rbind(ST_ExactMatchScore_series, c(thresh, ST_proportionExactMatch,"Ser/Thr"))
	Y_ExactMatchScore_series = rbind(Y_ExactMatchScore_series, c(thresh, Y_proportionExactMatch,"Tyr"))
}

# plot(Y_ExactMatchScore_series,type="n",axes=F,col=3, ylab = "Proportion of HPRD 2010-annotated sites correctly predicted in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(Y_ExactMatchScore_series,col='green',lwd=2.5) 
# par(new=T)
# plot(ST_ExactMatchScore_series,type='n',col=2, ylab = "Proportion of HPRD 2010-annotated sites correctly predicted in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(ST_ExactMatchScore_series,col='blue',lwd=2.5) 
# par(new=F)
# legend('topright', # places a legend at the appropriate place
# c('Tyr','Ser/Thr'), # puts text in the legend
# lty=c(1,1), # gives the legend appropriate symbols (lines)
# lwd=c(2.5,2.5),col=c('green','blue'))

# plot(Y_anyMatchScore_series,type="n",axes=F,col=3, ylab = "Proportion of HPRD 2010-annotated sites with any prediction in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(Y_anyMatchScore_series,col='green',lwd=2.5) 
# par(new=T)
# plot(ST_anyMatchScore_series,type='n',col=2, ylab = "Proportion of HPRD 2010-annotated sites with any prediction in PhosphoNet", xlab = "PhosphoNet score threshold")
# lines(ST_anyMatchScore_series,col='blue',lwd=2.5) 
# par(new=F)
# legend('topright', # places a legend at the appropriate place
# c('Tyr','Ser/Thr'), # puts text in the legend
# lty=c(1,1), # gives the legend appropriate symbols (lines)
# lwd=c(2.5,2.5),col=c('green','blue'))

ST_anyMatchScore_series = data.frame(ST_anyMatchScore_series)
colnames(ST_anyMatchScore_series) = c("thresh","proportion","motif_type")
Y_anyMatchScore_series = data.frame(Y_anyMatchScore_series)
colnames(Y_anyMatchScore_series) = c("thresh","proportion","motif_type")
ST_ExactMatchScore_series = data.frame(ST_ExactMatchScore_series)
colnames(ST_ExactMatchScore_series) = c("thresh","proportion","motif_type")
Y_ExactMatchScore_series = data.frame(Y_ExactMatchScore_series)
colnames(Y_ExactMatchScore_series) = c("thresh", "proportion","motif_type")

ExactMatchScore_series = rbind(ST_ExactMatchScore_series, Y_ExactMatchScore_series)
ExactMatchScore_series = data.frame(ExactMatchScore_series)
colnames(ExactMatchScore_series) = c("thresh", "proportion", "motif_type")
ExactMatchScore_series[,1] = as.numeric(as.character(ExactMatchScore_series[,1]))
ExactMatchScore_series[,2] = as.numeric(as.character(ExactMatchScore_series[,2]))

ppnet$residue = NA 
ppnet$residue[union(which(substr(ppnet[,3],1,1)=="S"),which(substr(ppnet[,3],1,1)=="T"))] = "Ser/Thr"
ppnet$residue[which(substr(ppnet[,3],1,1)=="Y")] = "Tyr"

grid.newpage()

# p1 = ggplot(ppnet, aes(x=ppnet$k1naScore, group=ppnet$residue, colour=ppnet$residue, fill=ppnet$residue)) + 
  # geom_density(alpha=0.2, position="identity") + theme_classic()
# p2 = ggplot(ExactMatchScore_series, aes(x = thresh, y = proportion, group = motif_type,
  # colour = motif_type)) + geom_line() + geom_point() + theme_bw() %+replace% 
    # theme(panel.background = element_rect(fill = NA), axis.ticks.y = element_blank())

p1 = ggplot(ppnet, aes(x=ppnet$k1naScore, group=ppnet$residue, colour=ppnet$residue, fill=ppnet$residue)) + 
  geom_density(alpha=0.2, position="identity") + theme_classic()
p2 = ggplot(ExactMatchScore_series, aes(x = thresh, y = proportion, group = motif_type,
  colour = motif_type)) + geom_line() + geom_point() + geom_text(aes(label=thresh),hjust=0, vjust=0) + theme_bw() %+replace% 
    theme(panel.background = element_rect(fill = NA), axis.ticks.y = element_blank()) + scale_x_continuous(name = "Percentile threshold")

# extract gtable
g1 <- ggplot_gtable(ggplot_build(p1))
g2 <- ggplot_gtable(ggplot_build(p2))

# overlap the panel of 2nd plot on that of 1st plot
pp <- c(subset(g1$layout, name == "panel", se = t:r))
g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
    pp$l, pp$b, pp$l)

# axis tweaks
ia <- which(g2$layout$name == "axis-l")
ga <- g2$grobs[[ia]]
ax <- ga$children[[2]]
ax$widths <- rev(ax$widths)
ax$grobs <- rev(ax$grobs)
ax$grobs[[1]]$x <- ax$grobs[[1]]$x + unit(1, "npc") - unit(0.15, "cm")
g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

# draw it
grid.draw(g)
