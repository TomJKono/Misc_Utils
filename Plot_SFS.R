#   An R script to plot a site frequency spectrum for various partitions of
#   SNPs. 
#   The data files are in the format
#   SNPNAME MAF
#   ...
#   separated by tabs
#   Fill in the correct values for your files, and your desired output file

#   Read the MAFs of the various SNP classes
SFS.1 <- read.table("FILE1", header=F)
SFS.2 <- read.table("FILE2", header=F)  
SFS.3 <- read.table("FILE3", header=F)
SFS.4 <- read.table("FILE4", header=F)

#   Bin them up!
bins = seq(0.0, 0.5, by=0.1)
SFS.1.bins <- cut(SFS.1[,2], breaks=bins, include.lowest=TRUE)
SFS.2.bins <- cut(SFS.2[,2], breaks=bins, include.lowest=TRUE)
SFS.3.bins <- cut(SFS.3[,2], breaks=bins, include.lowest=TRUE)
SFS.4.bins <- cut(SFS.4[,2], breaks=bins, include.lowest=TRUE)

#   Count up the frequencies
SFS.1.counts <- table(SFS.1.bins)/length(SFS.1[,2])
SFS.2.counts <- table(SFS.2.bins)/length(SFS.2[,2])
SFS.3.counts <- table(SFS.3.bins)/length(SFS.3[,2])
SFS.4.counts <- table(SFS.4.bins)/length(SFS.4[,2])

#   Stick them together, in the order that we want them to appear
SFS.alldata <- as.data.frame(cbind(SFS.1.counts, SFS.2.counts, SFS.3.counts, SFS.4.counts))

#   Open a handle to a PDF file for the plot
pdf(file="OUTPUT", width=8, height=6,family="Helvetica",pointsize=16)
#   Make the barplot, we have to do a separate x-axis later
barplot(t(SFS.alldata) ,ylim=c(0,0.35),beside=TRUE, axisnames=F,xlab="Allele Frequency",ylab="Proportion",col=c(grey(c(0.8, 0.6, 0.3, 0.0))))
#labels <- c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50")
labels <- c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5")

#   You'll have to play with these numbers to get the labels to line up with
#   the tick marks. As far as I know, there isn't a good way to figure this out
#   other than trial and error
at <- c(3, 8, 13, 18, 23)
axis(side=1, at=at, las=1,labels=labels,font=1,padj=1, cex.axis=0.8)
legend(inset=0,cex=0.55,"topright",c("CLASS 1", "CLASS 2", "CLASS 3", "CLASS 4"),fill=c(grey(c(0.8, 0.6, 0.3, 0.0))))
dev.off ()
