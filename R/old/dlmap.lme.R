`dlmap.lme` <-
function(object, phename, fixed=NULL, seed=1, maxit=60, multtest=c("holm", "bon"), alpha=.05, filestem="dl")
{
  if (!require(nlme)) 
	stop("nlme package must be installed to use this function. Please install it from CRAN before proceeding")

  fixed.forma <- NULL
  set.seed(seed)
  type <- attr(object, "type")
 
  if (missing(multtest)) multtest <- "holm"

  if (length(phename)!=1)
  	stop("Can only model one phenotypic trait")

  if (is.null(fixed))
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))

  if (!is.null(fixed)&(length(as.character(fixed))!=2))
  {
	print("Fixed formula has incorrect form: reverting to default of ~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }
  
  if (length(grep("\\|", as.character(fixed))>0))
	stop("Error: Fixed formula cannot contain grouping levels")

  if ((length(fixed)==2)&(as.character(fixed)[1]=="~"))
	fixed.forma <- as.formula(paste(phename, as.character(fixed)[1], as.character(fixed)[2], sep=""))

  idname <- object$idname

  object$envModel <- list()
  object$envModel$fixed <- fixed.forma
  object$dfMerged$grp1 <- 1
  nphe <- object$nphe

  object$alpha <- alpha
  object$multtest <- multtest

  #####
  # Recentering
  #####
  dfMdim <- ncol(object$dfMerged)-1
  n.mrk <- unlist(lapply(object$map, length))

  if (length(which(n.mrk==1))>0)
	message("Warning: Linkage groups ", which(n.mrk==1), " contain only one marker")

   object$dfMerged[, (nphe+1):dfMdim] <- t(t(object$dfMerged[,(nphe+1):dfMdim]) - apply(object$dfMerged[,(nphe+1):dfMdim], 2, function(x) return(mean(x,na.rm=TRUE))))

  message("Beginning detection stage...")
  det.out <- detect.lme(object, filestem)

  # List of no. detected QTL from output of detection step
  detqtl <- vector(length=length(object$map))

  for (i in 1:length(object$map))
    detqtl[i] <- length(grep(paste("C", i, "M", sep=""), det.out$qtls$pos))

  if (type=="f2") detqtl <- detqtl/2
  
  message("Detection stage complete. ", sum(detqtl), " total QTL detected")
  loc.out <- det.out
  if ((object$loc)&(sum(detqtl)>0))
  {
    message("Beginning localization stage...")
    loc.out <- localize.lme(object, QTLperChr=detqtl)
  }
  message("Fitting final model with all QTL")

  final.fixed <- fixed.forma
  if (sum(detqtl)>0) final.fixed <- paste(as.character(fixed.forma)[2], "~", as.character(fixed.forma)[3], "+", paste(loc.out$qtls$pos, collapse="+"), sep="")

  final.fixed <- as.formula(final.fixed)
  mod.fin <- lm(final.fixed, data=object$dfMerged, na.action=na.omit)

  output <- list()
  output$input <- object
  output$no.qtl <- detqtl
  output$final.model <- mod.fin
  output$final.model$call$fixed <- c(summary(mod.fin)$terms[[1]], summary(mod.fin)$terms[[2]], summary(mod.fin)$terms[[3]])
 
  if (sum(detqtl)>0)
  {
    table <- list()

    table[[1]] <- rep(names(object$map), detqtl)
    if (type=="f2") table[[1]] <- rep(table[[1]], each=2)

    mappos <- unlist(object$mapp)
     if (type=="f2")
        names(mappos) <- names(object$dfMerged)[(nphe+1):dfMdim][seq(1, dfMdim-nphe, 2)] else
     names(mappos) <- names(object$dfMerged)[(nphe+1):dfMdim]
     tmpord <- match(loc.out$qtls$pos, names(mappos))

    loc.out$qtls$pos <- loc.out$qtls$pos[order(tmpord)]

    if (type=="f2") pos <- unique(substr(loc.out$qtls$pos, 1, nchar(loc.out$qtls$pos)-1)) else 
  	pos <- loc.out$qtls$pos
    if (type=="f2") 	posD <- paste(pos, "D", sep="") else posD <- pos

    nmmap <- unlist(lapply(object$mapp, names))
    table[[2]] <- round(mappos[sort(tmpord)], 2)
    if (type=="f2") table[[2]] <- rep(table[[2]], each=2)

    table[[3]] <- table[[4]] <- NULL

    table[[5]] <- round(mod.fin$coefficients[match(loc.out$qtls$pos, names(mod.fin$coefficients))], 4)
    table[[6]] <- round(summary(mod.fin)$coefficients[match(loc.out$qtls$pos, rownames(summary(mod.fin)$coefficients)), 2], 4)

    table[[7]] <- round(table[[5]]/table[[6]], 2)
    table[[8]] <- round(sapply(table[[7]], function(x) 2*(1-pnorm(abs(x)))), 4)

    # Wald profile for each chromosome containing QTL
    output$profile <- loc.out$profile

    if (type!="other")
    {
	mark <- grep("M", names(object$dfMerged)[(nphe+1):dfMdim])+nphe
	if (type=="f2")
	  mark <- mark[seq(1, length(mark), 2)]

	endmrkL <- substr(pos, nchar(pos)-1, nchar(pos))=="M1"
	endmrkR <- substr(names(mappos)[match(posD, names(mappos))+1], nchar(pos)-1, nchar(pos))=="M1"
	endmrkR[posD==names(mappos)[length(mappos)]] <- TRUE

	lm <- vector(length=length(posD))
	rm <- vector(length=length(posD))
	lm[!endmrkL] <- sapply(match(posD[!endmrkL], names(object$dfMerged)), function(x) return(max(mark[mark<x])))
	rm[!endmrkR] <- sapply(match(posD[!endmrkR], names(object$dfMerged)), function(x) return(min(mark[mark>x])))
	
	lm[endmrkL] <- match(posD[endmrkL], names(object$dfMerged))
	rm[endmrkR] <- match(posD[endmrkR], names(object$dfMerged))

	table[[3]] <- nmmap[match(names(object$dfMerged)[lm], names(mappos))]
	table[[4]] <- nmmap[match(names(object$dfMerged)[rm], names(mappos))]
	
	if (type=="f2")
	{
	  table[[3]] <- rep(table[[3]], each=2)
	  table[[4]] <- rep(table[[4]], each=2)
	}
     }

     names(table) <- c("Chr", "Pos", "Left Marker", "Right Marker", "Effect", "SD", "Z-value", "p-value")
     output$Summary <- as.data.frame(do.call("cbind", table))
     #rownames(output$Summary) <- names(table$SD)
     rownames(output$Summary) <- NULL

   } # End of check for detected QTL
   
   if (sum(detqtl)==0)
 	print("No QTL detected in data. See log file for more details of testing.")
   class(output) <- c("dlmap", class(object))

   return(output)
}

