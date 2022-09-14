# rmf

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  print("Usage: Rscript scatterplot.R <youngVsOld table> <FBmirID> <FBmirName>");
  q();
}

inputFile <- args[1];
mirName <- args[2];

title = mirName;
xlabel<-"Zeitgeber Time";
ylabel<-"Expression (RPMMM)";
outPrefix = paste("/data/www/hendrixlab/html/data/circadian/microRNA_ATC/images_RPMM/",mirName,"_youngVsOld_replicates_44hrplot",sep="");

outputFile = paste(outPrefix,".pdf", sep="");
pdf(file=outputFile,useDingbats=F,family="Helvetica") 
par(mar=c(5, 4, 4, 0) + 0.1); #Extend margins. margin default=c(5, 4, 4, 2) + 0.1 # bottom, left, top and right
par(pin=c(5,5)); 

data = read.table(inputFile,col.name=c("col1","col2","col3","col4","col5"));

ZT <- c(0,4,8,12,16,20,24,28,32,36,40,44);

exprYoung <- c(data$col2,data$col3)
exprOld <- c(data$col4,data$col5)

dfYoung <- data.frame(ZT,exprYoung)
dfOld <- data.frame(ZT,exprOld)

Ymax<-max(data$col2, data$col3, na.rm=TRUE);
Omax<-max(data$col4, data$col5, na.rm=TRUE);

max=0;
if(Ymax > max) {
    max<-Ymax;
}
if(Omax > max) {
    max<-Omax;
}

maximum = 1.2*max;

xlim <- range(0,44);
ylim <- range(0,maximum);

#Not being used right now:
error.bars <- function(X,Y,Y0,Y1,w,col=1) {
X0 = X; X1 =X;
arrows(X0, Y0, X1, Y1, code=3,angle=90,length=w,col=col,lwd=1.5);
}

youngPlot <- plot(dfYoung,type = "l",col="red",lty=1,axes=F,ylim=ylim,xlim=xlim,xlab=xlabel,ylab=ylabel,main=title,bty='n',lwd=4,cex.lab=2.0,cex.main=3.5,font.main=3);
oldPlot <- lines(dfOld, type="l", pch=16, lty=2, col="blue",lwd=4);


tickLocations<-c(0,4,8,12,16,20,24,28,32,36,40,44);
axisLabels<-c(0,4,8,12,16,20,0,4,8,12,16,20);
axis(2,cex.axis=1.2);
axis(1,at=tickLocations,labels=axisLabels,cex.axis=.85);

legend('topleft',c('Young','Old'),col=c("red","blue"),lty=c(1,2),border=NA,bty = "n",lwd=4,cex=1.6);



supress <- dev.off(); 

