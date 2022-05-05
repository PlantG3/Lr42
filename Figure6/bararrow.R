#'@param dataf data frame with the numeric data in certain columns; each column contains data for a group
#'@param plotcols columns for plotting
bararrow <- function(dataf, plotcols, pmain = "", pxlab = "", pylab = "",
                     points=T, points.jitter.factor=3, points.col="gray80", points.cex=0.6,
                     pmain.cex = 1.5, plab.cex = 1.2,
                     barcolors = NULL, barnames = NULL, barnames.cex = 1, barnames.las = 1,
                     arrow.colors = NULL, arrow.lwd = 1, arrow.datatype = "se",
                     arrow.length = 0.1, arrow.code = 3,
                     ymin = 0, ymax = NULL, yaxis.cex = 1) {
  
  se <- function(x) {
  # standand errors
    sd(x, na.rm = T) / sqrt(length(x) - sum(is.na(x)))
  }
  
  # plot data
  pdata <- dataf[ , plotcols]
  pmeans <- apply(pdata, 2, mean, na.rm = T)
  
  # add individual data to a list
  ind.data <- vector("list", ncol(pdata))
  for (i in 1:ncol(pdata)) {
  	ind.data[[i]] <- pdata[, i]
  	print(head(ind.data[[i]]))
  }

  if (arrow.datatype == "sd") {
    parrow.data <- apply(pdata, 2, sd, na.rm = T)
  } else {
    parrow.data <- apply(pdata, 2, se)
  }
  
  if (is.null(barnames)) {
    barnames <- colnames(pdata)
  }
  
  if (is.null(barcolors)) {
    barcolors <- "grey"
  }
  
  if (is.null(arrow.colors)) {
    arrow.colors <- "grey"
  }
  
  if (is.null(ymax)) {
    ymax <- max(pmeans + parrow.data)
  }
  
  cat("ymax=", ymax, "\n")
  
  # bars
  barcenters <- barplot(pmeans, xlab = pxlab, ylab = pylab,
                        cex.lab = plab.cex, cex.main = pmain.cex,
                        main = pmain, ylim = c(ymin, ymax),
                        names.arg = barnames,
                        las = barnames.las,
                        col = barcolors,
                        cex = barnames.cex,
                        cex.axis = yaxis.cex)
  
  # points
  if (points) {
  	for (i in 1:ncol(pdata)) {
  	  data.points <- ind.data[[i]]
  	  npoints <- length(data.points)
		  points(jitter(rep(barcenters[i], npoints), factor=points.jitter.factor),
		         data.points, col=points.col, cex=points.cex, xpd=T)
  	}
	}

  # arrows
  arrows(barcenters, pmeans - parrow.data, barcenters, pmeans + parrow.data,
         angle = 90, length = arrow.length, code = 3, lwd = arrow.lwd)
  
  invisible(barcenters)
}
