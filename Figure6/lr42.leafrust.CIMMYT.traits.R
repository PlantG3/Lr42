setwd("/bulk/liu3zhen/research/projects/Lr42/Lr42_fromLab83/lr42/CIMMYT_traits")

source("~/scripts/visualization/bararrow.R")

td <- read.delim("data_from_Juliana.txt", stringsAsFactors = F) # trait data

colnames(td)
td$Leafrust_seedling_MBJSP
td$Lr42

# seedling
seedlinglr42p <- td$Leafrust_seedling_MBJSP[!is.na(td$Leafrust_seedling_MBJSP) & td$Lr42 == 1]
sum(!is.na(td$Leafrust_seedling_MBJSP) & td$Lr42 == 1)
seedlinglr42a <- td$Leafrust_seedling_MBJSP[!is.na(td$Leafrust_seedling_MBJSP) & td$Lr42 == 0]
sum(!is.na(td$Leafrust_seedling_MBJSP) & td$Lr42 == 0)
seedling.num.diff <- length(seedlinglr42a) - length(seedlinglr42p)
seedling.plotdf <- data.frame(Plus = c(seedlinglr42p, rep(NA, seedling.num.diff)), Minus = seedlinglr42a)
p1 <- t.test(seedlinglr42p, seedlinglr42a)[[3]]
p1 <- format(p1, digits = 2)

# APR
adultlr42p <- td$Leafrust_Obregon[!is.na(td$Leafrust_Obregon) & td$Lr42 == 1]
sum(!is.na(td$Leafrust_Obregon) & td$Lr42 == 1)
adultlr42a <- td$Leafrust_Obregon[!is.na(td$Leafrust_Obregon) & td$Lr42 == 0]
sum(!is.na(td$Leafrust_Obregon) & td$Lr42 == 0)
adult.num.diff <- length(adultlr42p) - length(adultlr42a)


adultlr42p2 <- td$Leafrust_ElBatan[!is.na(td$Leafrust_ElBatan) & td$Lr42 == 1]
sum(!is.na(td$Leafrust_ElBatan) & td$Lr42 == 1)
adultlr42a2 <- td$Leafrust_ElBatan[!is.na(td$Leafrust_ElBatan) & td$Lr42 == 0]
sum(!is.na(td$Leafrust_ElBatan) & td$Lr42 == 0)
adult.num.diff2 <- length(adultlr42p2) - length(adultlr42a2)
adult.plotdf2 <- data.frame(Plus = adultlr42p2, Minus = c(adultlr42a2, rep(NA, adult.num.diff2)))

max.entries <- max(length(adultlr42p), length(adultlr42a), length(adultlr42p2), length(adultlr42a2))
adult.plotdf <- data.frame(Plus1 = c(adultlr42p, rep(NA, max.entries - length(adultlr42p))),
                           Minus1 = c(adultlr42a, rep(NA, max.entries - length(adultlr42a))),
                           Plus2 = c(adultlr42p2, rep(NA, max.entries - length(adultlr42p2))),
                           Minus2 = c(adultlr42a2, rep(NA, max.entries - length(adultlr42a2)))
                           )

# ANOVA
adult.elbatan <- td[!is.na(td$Leafrust_ElBatan), c("Lr42", "Leafrust_ElBatan")]
colnames(adult.elbatan) <- c("Lr42", "Severity")
adult.elbatan$Location <- "ElBatan"
adult.obregon <- td[!is.na(td$Leafrust_Obregon), c("Lr42", "Leafrust_Obregon")]
colnames(adult.obregon) <- c("Lr42", "Severity")
adult.obregon$Location <- "Obregon"

adult.data <- rbind(adult.elbatan, adult.obregon)
head(adult.data)
adult.data$Lr42 <- as.factor(adult.data$Lr42)
aov.res <- aov(Severity ~ Lr42 + Location, data = adult.data)
aov.summary <- summary(aov.res)

p2 <- aov.summary[[1]][1, 5]
p2 <- format(p2, digits = 2)


# plot
#pdf("lr42.leafrust.CIMMYT.traits.1.pdf", width=2.5, height=3)
par(mar = c(3, 4, 4, 1), mfrow = c(1, 1))
bararrow(dataf = seedling.plotdf, plotcols = 1:2, pmain = "Seedling Leaf Rust",
         pxlab = "", pylab = "Disease severity",
         plab.cex = 1.1, pmain.cex = 1, barcolors = c("dark green", "blue"),
         barnames = c("Lr42+", "Lr42-"), barnames.cex = 1, barnames.las = 1,
         arrow.colors = NULL, arrow.lwd = 1, arrow.datatype = "se",
         arrow.length = 0.2, arrow.code = 3,
         ymin = 0, ymax = NULL, yaxis.cex = 1)
legend("topleft", legend = paste0("p=", p1), bty = "n")
#dev.off()

#pdf("lr42.leafrust.CIMMYT.traits.1.pdf", width=4, height=3)
bararrow(dataf = adult.plotdf, plotcols = 1:4, pmain = "Adult Leaf Rust", pxlab = "", pylab = "Disease severity",
         plab.cex = 1.1, pmain.cex = 1,
         barcolors = c("dark green", "blue", "dark green", "blue"),
         barnames = c("Lr42+", "Lr42-", "Lr42+", "Lr42-"), barnames.cex = 1, barnames.las = 1,
         arrow.colors = NULL, arrow.lwd = 1, arrow.datatype = "se",
         arrow.length = 0.15, arrow.code = 3,
         ymin = 0, ymax = NULL, yaxis.cex = 1)
legend("topleft", legend = paste0("p=", p2), bty = "n")
#dev.off()

