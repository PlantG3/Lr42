setwd(".")

## plot function
oneplot_pdf<-function (input, yrange = NULL, xlab.text="Chromosome",
                       ylab.text = "Probability of linkage",
                       main.text="",axisline.width=1.5, chr.set,
                       chr.size, order.by.chrsize=F, label.rm=NULL, cexaxis=0.5,
                       saveplot=FALSE, plot.width=6, plot.heigh=4.5,
                       xaxis.draw = F,
                       plot.path=".", plot.filename="default.pdf") {
        
        # input should have three columns, "Chr", "Pos", "Prob"
        # size contains length information for every chromosome or contigs
        # size has two columns: chr and size
        plot.path <- gsub("/$", "", plot.path)
        colnames(input) <- c("Chr", "Pos", "Prob")
        input$Chr <- as.character(input$Chr)
        input <- input[input$Chr %in% chr.set, ]
        
        # chromosome length:
        colnames(chr.size) <- c("Chr", "Size")
        chr.size$Chr <- as.character(chr.size$Chr)
        chr.size <- chr.size[chr.size$Chr %in% chr.set, ]
        if (order.by.chrsize) {
                chr.size <- chr.size[order(chr.size$Size, decreasing=T), ]
        }
        # to judge the number is odd or even:
        odd <- function(x) {
                if ((round(x/2,0) - i/2)==0) {
                        y <- 0
                } else {
                        y <- 1
                }
                return(y)
        }
        
        # plotting
        accum <- 0
        all.col <- NULL
        all.chr <- NULL
        centers <- NULL
        gap <- sum(chr.size$Size)/100
        if (is.null(yrange)) {
                ymax <- max(input$Prob)
                yrange <- c(0, ymax)
        } else {
                ymax <- max(yrange)
        }
        if (saveplot) {
                pdf(paste(plot.path, plot.filename, sep="/"), width=plot.width, heigh=plot.heigh,pointsize = 8)
                #png(paste(plot.path, plot.filename, sep="/"), pointsize = 4, res = 2000, units = "in", 
                #   width=plot.width, heigh=plot.heigh)
        }
        
        plot(NULL, NULL, ylim=yrange,
             xlim=c(0, gap*nrow(chr.size)+sum(chr.size$Size)),
             xaxt="n",xlab=xlab.text, ylab=ylab.text,
             main=main.text)
        #box(col = 'gray44')
        all.accum <- NULL
        for (i in 1:(nrow(chr.size))) {
                all.accum <- c(all.accum, accum)
                pre.accum <- accum
                chr <- chr.size[i, "Chr"]
                len <- chr.size[i, "Size"]
                if (odd(i)) {
                        plot.col = 'blue'
                } else {
                        plot.col = "dark green"
                }
                pos <- input[input$Chr==chr, "Pos"]
                prob <- input[input$Chr==chr, "Prob"]
                if (xaxis.draw) {
                        lines(c(accum, accum+len), c(-(ymax/50), -(ymax/50)), col=plot.col, lwd=axisline.width, lend=1)
                }
                points(accum+pos, prob, pch=16, cex=0.5, col=plot.col)
                accum <- accum + len + gap
                center.point <- (pre.accum + accum - gap)/2
                all.col <- c(all.col, plot.col)
                all.chr <- c(all.chr, chr)
                centers <- c(centers, center.point)
        }
        if (!is.null(label.rm)) {
                for (each.label.rm in label.rm) {
                        all.chr <- gsub(each.label.rm, "", all.chr)
                }
        }
        print(all.chr)
        axis(side=1, at=centers, labels=all.chr, tick=F, cex.axis=cexaxis)
        if (saveplot) dev.off()
        
        names(all.accum) <- chr.set
        return(all.accum)
}





### plot 1
outfile <- "Sourcedata_Figure1c_BSR.txt"
probs <- read.delim(outfile, stringsAsFactors = F)

maxlength <- tapply(probs$POS, probs$CHR, max)
chrsize <- data.frame(Chr = names(maxlength), Length = maxlength)
chrsize$Length <- as.numeric(as.character(chrsize$Length))
chrset<-chrsize$Chr[-8]
oneplot_pdf(input = probs[, c("CHR", "POS", "ppp")],
        ylab.text = "Probability of complete linkage",
        main.text = "TA2450xTA2433 BSR-seq",
        chr.set = chrset, cexaxis = 1,
        chr.size = chrsize,
        order.by.chrsize = F, label.rm = "Chr",
        saveplot = T, plot.width = 3, plot.heigh = 3,
        plot.path = ".", plot.filename = "Figure1c.pdf")






### pop2

outfile2 <- "Sourcedata_SupplementaryFigure1c_BSR.txt"

probs2 <- read.delim(outfile2, stringsAsFactors = F)

### plot 2
oneplot_pdf(input = probs2[, c("CHR", "POS", "ppp")],
        ylab.text = "Probability of complete linkage",
        main.text = "TA2450xTA10132 BSR-seq",
        chr.set = chrset,
        chr.size = chrsize, cexaxis = 1,
        order.by.chrsize = F, label.rm = "Chr",
        saveplot = T, plot.width = 3, plot.heigh = 3,
        plot.path = ".", plot.filename = "SupplementaryFigure1c.pdf")

