# tcga

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Convert aid to R aid
#'
#' This function converts aid into R aid
#' @param x A vector of AIDs
#' @keywords aid r.aid
#' @export
#' @examples
#' aid2r.aid()
aid2r.aid = function (x) {
  x <- gsub("-", ".", x)
  first.chr <- suppressWarnings(as.numeric(substr(x,1,1)))
  x[which(first.chr != "NA")] <- paste("X", x[which(first.chr != "NA")], sep="")
  return(x)
}

#' Convert R aid to aid
#'
#' This function converts R AIDs into AIDs
#' @param x A vector of R AIDs.
#' @keywords aid r.aid
#' @export
#' @examples
#' r.aid2aid()
r.aid2aid = function (x) {
  x <- gsub(".", "-", x, fixed=TRUE)
  x <- sub("^X", "", x)
  return(x)
}

#' Get AOV pvalue
#'
#' This function returns the -log10(pvalue) of AOV
#' @param x Continuous data vector
#' @param g Group labels
#' @keywords aov
#' @export
#' @examples
#' func.anova(x=c(3,4,5,10,12,15,100,120,130), g=c("A","A","A","B","B","B","C","C","C"))
func.anova = function(x, g) {
  -log(unlist(summary(aov(as.numeric(x)~as.factor(g))))[9],base=10)
}

#' Refactor
#'
#' This function refactors data
#' @param data
#' @param startdata
#' @param annotation
#' @param colclasses
#' @keywords refactor
#' @export
#' @examples
#' refactor()
refactor = function (data, startdata, annotation, colclasses) {
  d   = matrix(nrow=nrow(data), ncol=0)
  grp = c()
  for (fac in levels(annotation[[colclasses]])) {
    index = which(colnames(data) %in% annotation$name[which(annotation[[colclasses]]==fac)])
    grp = c(grp, rep(fac, length(index)))
    d   = cbind(d, data[,index])
  }

  colnames(d) = r.aid2aid(colnames(d))

  # remove rows that contain all zeros
  zi = c()
  for(i in 1:nrow(d)) { zi[i] <- all(d[i,]==0) }
  d <- d[which(zi==FALSE),]

  # calculate Pvalues
  Pval <- apply(d, 1, func.anova, g=grp)

  # get gene annotation
  gene <- data[which(zi==FALSE),1:startdata-1]

  d <- cbind(gene, Pval, d)

  return (list("data" = d, "groups" = grp, "zero.index" = zi))
}

#' Sample types to Factors
#'
#' This function converts numerical TCGA Sample types to a factor by replacing with "T" or "N"
#' @param v Vector of Sample types
#' @keywords sample type
#' @export
#' @examples
#' sampleTypeToFactor(c(1,2,3,6,10,11))
sampleTypeToFactor = function (v) {
  t = data.frame(num=c(1,2,3,6,10,11), type=c("T","T","T","T","N","N"))
  f = t[match(v, t$num),]$type
  return (f)
}

