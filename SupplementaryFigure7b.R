setwd(".")
# install if not a package has not installed yet
list.of.packages <- c("dplyr", "multcompView")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# loading package
library(dplyr)
library(multcompView)

# data path
datapath <- "./Lr42/qRT_PCR/VIGS/"
d <- read.delim(paste0(datapath, "/qRTPCR_VIGS.txt"), stringsAsFactors=F)
### expression formula
#100 x 2^(Cp-rl - Cp-GOI)
# GOI: gene of interest
d$SampleID <- d$Sample
d$Sample <- gsub("_[123]$", "", d$SampleID)
d$Cq[is.na(d$Cq)] <- 40 # NA to 40

# rl
rl0 <- d[d$Target == "RL1", c("Cq", "Sample")]
rl <- rl0 %>% group_by(Sample) %>% summarise_at(vars("Cq"), mean)
colnames(rl) <- c("Sample", "rlCq")

# genes of interest and average technical reps
goi0 <- d[d$Target == "Lr42", c("Target", "Cq", "Sample", "SampleID")]
goi <- goi0 %>% group_by(Target, Sample) %>% summarise_at(vars("Cq"), mean)
goi2 <- merge(goi, rl, by = "Sample")
goi2$Expression <- 100 * 2 ^ (goi2$rlCq - goi2$Cq)
goi2

# output
write.table(goi2, "VIGS.qRTPCR.plotting.expdata.txt", quote=F, row.names=F, sep="\t")

# average biological reps:
goi2$Group <- gsub("-[123]", "", goi2$Sample)
exp.mean <- tapply(goi2$Expression, goi2$Group, mean, na.rm=T)
exp.sd <- tapply(goi2$Expression, goi2$Group, sd, na.rm=T)
label.names <- names(exp.mean)

names(exp.mean) <- label.names
names(exp.sd) <- label.names
exp.mean
exp.sd

label.names.new <- c("LIC", "Lr42")
label.renames <- c("LIC", "Lr42")
#######################################################################################
### lr42
#######################################################################################
lr42.exp.mean <- exp.mean[label.names.new]
lr42.exp.sd <- exp.sd[label.names.new]

names(lr42.exp.mean) <- label.renames
names(lr42.exp.sd) <- label.renames 

lr42.ymax.bar  <- max(lr42.exp.mean + lr42.exp.sd)
lr42.ymax <- lr42.ymax.bar * 1.05
#lr42.ymax <- 12

### plot

### plot
pdf("VIGS.qRTPCR.pdf", width = 3, height = 3)

par(mar=c(3,3,3,1))
barcenters <- barplot(lr42.exp.mean, ylab = "Relative expression", las = 1,
                      names.arg = label.renames, ylim = c(0, lr42.ymax),
                      main = substitute(paste(italic(Lr42), " expression")),
                      col = c("lightblue4", "lightblue", "lightblue", "gray50"))
arrows(barcenters, lr42.exp.mean - lr42.exp.sd, barcenters,
       lr42.exp.mean + lr42.exp.sd, angle = 90,
       length = 0.1, code = 3)

# ttest
lic <- goi2[grep("LIC", goi2$Sample), "Expression"]
lr42 <- goi2[grep("Lr42", goi2$Sample), "Expression"]
ttest <- t.test(lic, lr42)
pval <- round(ttest[[3]],3)

# add expressions of individuals
for (i in 1:length(lr42.exp.mean)) {
        en <- names(lr42.exp.mean)[i]
        eexp <- goi2[grep(en, goi2$Sample), "Expression"]
        barpos <- barcenters[i]
        dotxpos <- jitter(rep(barpos, length(eexp)), factor=2)
        points(dotxpos, eexp, pch=19, col="gray30", cex=2, xpd=T)
}

legend("topright", legend=paste0("p=", pval), bty='n')

dev.off()



# What is the effect of the treatment on the value ?
lr42exp <- goi2$Expression[goi2$Target == "Lr42_qrt_FR6"]
lr42samples <- goi2$Sample[goi2$Target == "Lr42_qrt_FR6"]
lr42samples <- gsub("\\-", "_", lr42samples)
lr42model <- lm(lr42exp ~ lr42samples)
lr42aov <- aov(lr42model)
summary(lr42aov)

# Tukey test to study each pair of treatment :
lr42.tukey <- TukeyHSD(x=lr42aov, 'lr42samples', conf.level=0.95)
lr42levels <- lr42.tukey[["lr42samples"]][,4] %>% multcompLetters
lr42levels <- lr42levels[['Letters']]
lr42labels <- lr42levels
names(lr42levels) <- gsub("_", "-", names(lr42labels))
lr42labels <- lr42levels[label.names.new]
text(barcenters, lr42.exp.mean + lr42.exp.sd, label=lr42labels, pos=3, xpd=T)

lr42.tukey 

dev.off()

