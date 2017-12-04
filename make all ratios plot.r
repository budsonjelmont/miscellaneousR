#plot peak area ratios (10x vs 1x, 10x vs 3x, 3x vs 1x)  from accuracy dataset

p <- ggplot(rbind(RatioNF101, RatioNF103, RatioNF31), aes(x=peakRatios, fill=pass)) + geom_density(alpha=0.2, position="identity") +
	scale_x_continuous(breaks=seq(-10,10,1)) +
	labs(x=expression(log[e]~" peak area ratios"))
ggsave(filename=paste("All peak area ratios.jpg",sep=""), plot=p)

#plot all ratio comparisons together
p <- ggplot(rbind(RatioNF101, RatioNF103, RatioNF31), aes(x=peakRatios, fill=pass)) + geom_density(alpha=0.2, position="identity") +
	scale_x_continuous(breaks=seq(-10,10,1)) +
	labs(x=expression(log[e]~" peak area ratios")) +
	geom_vline(xintercept=log(10/1),linetype = "longdash") +
	geom_vline(xintercept=log(10/3),linetype = "longdash") +
	geom_vline(xintercept=log(3/1),linetype = "longdash") +
	geom_text(aes(log(10/1),0.53,label = signif(log(10/1),5), hjust=1.25, vjust = 0.25), size = 3.5) +
	geom_text(aes(log(10/3),1.25,label = signif(log(10/3),5), hjust=1.25, vjust = 0.25), size = 3.5) +
	geom_text(aes(log(3/1),0.8,label = signif(log(3/1),5), hjust=1.25, vjust = 0.25), size = 3.5)
ggsave(filename=paste("All peak area ratios_Vline.jpg",sep=""), plot=p)

#plot just a single ratio comparison
p <- ggplot(RatioNF31, aes(x=peakRatios, fill=pass)) + geom_density(alpha=0.2, position="identity") +
	scale_x_continuous(breaks=seq(-10,10,1)) +
	labs(x=expression(log[e]~" 3x/1x peak area ratios")) +
	geom_vline(xintercept=log(3/1),linetype = "longdash") +
	geom_text(aes(log(3/1),1.25,label = signif(log(3/1),5), hjust=1.25, vjust = 0.25), size = 3.5) +
	ggtitle("3X vs 1X Peak area ratios")
ggsave(filename=paste("3Xvs1X peak area ratios_Vline.jpg",sep=""), plot=p)