setwd(".")

#install.packages("PopGenome")
library("PopGenome")

GENOME.class <- readData("FASTA")
GENOME.class@n.sites
GENOME.class@region.names
get.sum.data(GENOME.class)
GENOME.class@region.data
GENOME.class  <- neutrality.stats(GENOME.class)
get.neutrality(GENOME.class)
get.neutrality(GENOME.class)[[1]]

GENOME.class@Tajima.D
GENOME.class  <- F_ST.stats(GENOME.class)
GENOME.class  <- diversity.stats(GENOME.class)
get.diversity(GENOME.class)[[1]]
GENOME.class.slide  <- sliding.window.transform(GENOME.class,width=50,jump=10,start.pos = 50,end.pos = 2700,type=2,whole.data = T)
GENOME.class@n.sites
GENOME.class.slide@region.names
GENOME.class.slide<-diversity.stats(GENOME.class.slide)

get.diversity(GENOME.class.slide)[[1]]
pos<-c()
n<-length(get.diversity(GENOME.class.slide)[[1]][,1])
j=0

for (i in 1:n){
  pos<-c(pos,as.integer((2*j+1+50)/2))
  j<-j+10
}

pdf("Figure3c.pdf",width=6,height = 5)
  plot(pos,get.diversity(GENOME.class.slide)[[1]][,1],xlab="Position(bp); win=50, step=10",ylab = "Nucleotide Diversity",cex=0.3,cex.lab=1,col="Blue",pch=19)
  lines(x = c(1,399),y=c(-0.2,-0.2),col="blue",lwd=2)
  lines(x = c(511,1374),y=c(-0.2,-0.2),col="blue",lwd=2)
  lines(x = c(1741,2794),y=c(-0.2,-0.2),col="blue",lwd=2)
  dev.off()
  write.table(get.diversity(GENOME.class.slide)[1], "Figure3c.txt",sep="\t",quote = F)
