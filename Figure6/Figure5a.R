setwd("/data1/home/liu3zhen/lr42/CIMMIT/")
options(stringsAsFactors = F)

counts0 <- read.delim("./data/Lr42_cimmyt.txt")
libsize <- read.delim("./data/cimmyt2017-2018_ReadsPerSample.txt")
cimmit <- read.csv("./data/Lr42postulation_Sandesh_SM.csv")
head(cimmit)
nrow(cimmit)
plr42lines <- cimmit$Line[!is.na(cimmit$LR42Postulation) & cimmit$LR42Postulation == 1]

counts <- merge(counts0, libsize[, 1:2], by.x="Lines", by.y = "FullSampleName")
nrow(counts)
head(counts)

nonzerocounts <- function(x) { sum(x > 0) }

LineTagC <- apply(counts[, -c(1, ncol(counts))], 1, nonzerocounts)
TagLineC <- apply(counts[, -c(1, ncol(counts))], 2, nonzerocounts)
hist(LineTagC)
hist(TagLineC)
nrow(counts[LineTagC > 20, ])

###########################################################################################
par(mfrow = c(1, 2))
log10size <- log10(counts$goodBarcodedReads + 1)
lineTagC.jitter <- jitter(LineTagC, factor = 2.5)

plot(log10size, lineTagC.jitter,
     xlab = "log10(Number of reads)",
     ylab = "Number of target tags",
     cex = 0.4, lwd = 0.2,
     col = "dodgerblue3", main = "CIMMIT lines")

plot(log10size, lineTagC.jitter,
     xlab = "log10(Number of reads)",
     ylab = "Number of target tags",
     cex = 0.4, lwd = 0.2,
     col = "dodgerblue3",
     xlim = c(4,7), main = "CIMMIT lines + Lr42 introgression highlighted")
points(log10size[counts$Lines %in% plr42lines], lineTagC.jitter[counts$Lines %in% plr42lines],
       cex = 0.4, lwd = 0.3, col = "red")
###########################################################################################

tags.tauschii <- read.delim("./data/TA2450.tags.txt")
tags.tauschii.sub <- tags.tauschii[tags.tauschii$Tag %in% colnames(counts), ]

plot(tags.tauschii.sub$num_lines, TagLineC[tags.tauschii.sub$Tag])


more.specific.tags <- intersect(names(TagLineC[TagLineC < 5000]), tags.tauschii.sub[tags.tauschii.sub$num_lines< 100, "Tag"])
length(more.specific.tags)


TagLineC[more.specific.tags]

LineTagCsub <- apply(counts[, more.specific.tags], 1, nonzerocounts)
hist(LineTagCsub)
hist(LineTagCsub, ylim = c(0, 1000))
sum(LineTagCsub >= 6)




LineTagC <- apply(counts[, more.specific.tags], 1, nonzerocounts)
TagLineC <- apply(counts[, more.specific.tags], 2, nonzerocounts)
hist(LineTagC)
hist(TagLineC)
nrow(counts[LineTagC > 20, ])

###########################################################################################
par(mfrow = c(1, 2))
log10size <- log10(counts$goodBarcodedReads + 1)
lineTagC.jitter <- jitter(LineTagC, factor = 2.5)

plot(log10size, lineTagC.jitter,
     xlab = "log10(Number of reads)",
     ylab = "Number of target tags",
     cex = 0.4, lwd = 0.2,
     col = "dodgerblue3", main = "CIMMIT lines")

plot(log10size, lineTagC.jitter,
     xlab = "log10(Number of reads)",
     ylab = "Number of target tags",
     cex = 0.4, lwd = 0.2,
     col = "dodgerblue3",
     xlim = c(4,7), main = "CIMMIT lines + Lr42 introgression highlighted")
points(log10size[counts$Lines %in% plr42lines], lineTagC.jitter[counts$Lines %in% plr42lines],
       cex = 0.4, lwd = 0.3, col = "red")
###########################################################################################
