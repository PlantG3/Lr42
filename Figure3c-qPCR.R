setwd(".")
library(dplyr)
#install.packages("multcompView")
library(multcompView)

datapath <- "./Lr42/qRT_PCR/native_transgene/"
d <- read.delim(paste0(datapath, "/qRTPCR10012021.txt"), stringsAsFactors=F)
### expression formula
#100 x 2^(Cp-actin - Cp-GOI)
# GOI: gene of interest
d$SampleID <- d$Sample
d$Sample <- gsub("-[123]$", "", d$SampleID)
d <- d[d$Sample != "NTC" & d$Sample != "TA2450", ]


d$Cq[is.na(d$Cq)] <- 40 # NA to 40

# actin
actin0 <- d[d$Target == "Actin", c("Cq", "SampleID")]
actin <- actin0 %>% group_by(SampleID) %>% summarise_at(vars("Cq"), mean)
colnames(actin) <- c("SampleID", "actinCq")

# genes of interest
goi0 <- d[d$Target == "Lr42_qrt_FR6", c("Target", "Cq", "Sample", "SampleID")]
goi <- goi0 %>% group_by(Target, Sample, SampleID) %>% summarise_at(vars("Cq"), mean)
goi2 <- merge(goi, actin, by = "SampleID")
goi2$Expression <- 100 * 2 ^ (goi2$actinCq - goi2$Cq)
goi2
# output
write.table(goi2, "nativeProm.qRTPCR.plotting.expdata.txt", quote=F, row.names=F, sep="\t")

# average biological reps:
exp.mean <- tapply(goi2$Expression, paste(goi2$Targe, goi2$Sample), mean, na.rm=T)
exp.sd <- tapply(goi2$Expression, paste(goi2$Targe, goi2$Sample), sd, na.rm=T)
label.names <- names(exp.mean)

names(exp.mean) <- label.names
names(exp.sd) <- label.names

label.names.new <- c("Bobwhite", "G18W019-4", "933-1", "968-3")
label.renames <- c("Bobwhite", "Ubi::Lr42", "Lr42P::Lr42 1", "Lr42P::Lr42 2")

#######################################################################################
### lr42
#######################################################################################
lr42.exp.mean <- exp.mean[grep("Lr42", names(exp.mean))]
lr42.exp.sd <- exp.sd[grep("Lr42", names(exp.sd))]

names(lr42.exp.mean) <- gsub("Lr42_qrt_FR6 ", "", names(lr42.exp.mean))
names(lr42.exp.sd) <- gsub("Lr42_qrt_FR6 ", "", names(lr42.exp.sd))
lr42.exp.mean <- lr42.exp.mean[label.names.new]
lr42.exp.sd <- lr42.exp.sd[label.names.new]

names(lr42.exp.mean) <- label.renames
names(lr42.exp.sd) <- label.renames 

#######################################################################################
### lr42
#######################################################################################

bottom_ceil <- 12
top_floor <- 45
top_ceil <- 62

total_height <- top_ceil - (top_floor - bottom_ceil)
gap <- total_height / 25
total_height <- total_height + gap

nbar <- length(lr42.exp.mean)
bar_width <- 1
bar_dist <- bar_width / 5
xmax <- bar_width * nbar + bar_dist * (nbar - 1)
bar_cols <- c("gray50", "lightblue4", "lightblue", "lightblue")

pdf("transgene.qRTPCR.gapped.pdf", width=3.5, height=5)
par(mar=c(4, 4.5, 3, 1))
plot(NULL, NULL, ylim=c(0, total_height), xlim=c(0, xmax), axes=F, bty="n",
     main=substitute(paste(italic(Lr42), " expression")),
     xlab="", ylab="Relative expression")

# bottom figure
lr42.exp.mean.bottom <- pmin(lr42.exp.mean, bottom_ceil)
bar <- 1
xstart <- 0
for (exp in lr42.exp.mean.bottom) {
  rect(xstart, 0, xstart + bar_width, exp, col=bar_cols[bar], border=NA)
  
  group_name <- label.names.new[bar]
  individual_exp <- goi2$Expression[goi2$Sample == group_name]
  individual_exp_bottom <- individual_exp[individual_exp <= bottom_ceil]
  points(jitter(rep(xstart,length(individual_exp_bottom)) + bar_width/2, amount=0.3),
         individual_exp_bottom, col="gray30", pch=19)
  
  if (lr42.exp.mean[bar] - lr42.exp.sd[bar] < bottom_ceil) {
    arrows(xstart + bar_width/2, lr42.exp.mean[bar],
           xstart + bar_width/2, lr42.exp.mean[bar] + lr42.exp.sd[bar],
           angle = 90, length = 0.1, code = 2)
    
  }
  # xlabels
  text(xstart + bar_width/2, -0.8, pos=2, labels=names(lr42.exp.mean.bottom)[bar], srt=30, xpd=T)
  # next
  xstart <- xstart + bar_width + bar_dist
  bar <- bar + 1
}
axis(2, at=c(0, bottom_ceil), las=2, lwd=1.2)

# top figure
bar <- 1
xstart <- 0
for (exp in lr42.exp.mean) {
  if (exp>bottom_ceil) {
    bar_bottom <- bottom_ceil + gap
    bar_top <- exp - top_floor + bottom_ceil + gap
    rect(xstart, bar_bottom , xstart + bar_width, bar_top, col=bar_cols[bar], border=NA)
    
    # individual dots
    group_name <- label.names.new[bar]
    individual_exp <- goi2$Expression[goi2$Sample == group_name]
    
    individual_exp_top <- individual_exp[individual_exp >= top_floor]
    points(jitter(rep(xstart,length(individual_exp_top)) + bar_width/2, amount=0.3),
           individual_exp_top - top_floor + bottom_ceil + gap,
           col="gray30", pch=19)
    
    individual_exp_top <- individual_exp[individual_exp < top_floor & individual_exp > bottom_ceil]
    if (length(individual_exp_top)>0) {
      points(jitter(rep(xstart, length(individual_exp_top)) + bar_width/2, amount=0.5),
             rep(bottom_ceil + gap/2, length(individual_exp_top)),
             col="gray30", pch=19)
    }
    
    
    if (exp - lr42.exp.sd[bar] > bottom_ceil) {
      arrows(xstart + bar_width/2, bar_top,
             xstart + bar_width/2, bar_top + lr42.exp.sd[bar],
             angle = 90, length = 0.1, code = 2)
    }
  }
  # next
  xstart <- xstart + bar_width + bar_dist
  bar <- bar + 1
}
axis(2, at=c(bottom_ceil + gap, total_height), labels=c(top_floor, top_ceil), las=2, lwd=1.2)

dev.off()
