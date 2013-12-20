`dlmap.asreml` <-
function(object, phename, baseModel, fixed=NULL, random=NULL, rcov=NULL, sparse=NULL, pedigree, seed=1, n.perm=0, multtest=c("holm", "bon"), alpha=.05, filestem="dl", ...)
{
  idname <- object$idname
  type <- attr(object, "type")

  if (!require(asreml))
	stop("ASReml must be installed to use this function. Please use dlmap.lme if you do not have a license")

  if (missing(object)) stop("Input object is required to perform dlmap analysis")
  if (!inherits(object, "dlcross"))  stop("Input object should have class dlcross, see create.dlcross")

  trace <- paste(filestem, ".trace", sep="")
  ftrace <- file(trace, "w")
  sink(trace, type = "output", append = FALSE)
  on.exit(sink(type = "output"))
  on.exit(close(ftrace), add = TRUE)

  set.seed(seed)

  # Check objects
  if (missing(multtest))  multtest <- "holm"

  if (missing(phename))
	stop("Phenotypic trait name is a required argument and is missing")

  if (length(phename)!=1)
	stop("Error: Can only analyze one trait")
 
  fixed.forma <- NULL
  if (is.null(fixed))
  {
	message("Warning: No fixed effects object for spatial model; will
		fit the default model: ", phename,"~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }

  if ((!is.null(fixed))&(length(as.character(fixed))!=2))
  {
	message("Fixed formula has incorrect form: reverting to default of ~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }

  if ((length(fixed)==2)&(as.character(fixed)[1]=="~"))
	fixed.forma <- as.formula(paste(phename, as.character(fixed)[1], as.character(fixed)[2], sep=""))

  if (is.na(match(phename, colnames(object$dfMerged))))
	stop("The response ", phename, " is not included in the data; please check names of variables")
 
  if (!missing(pedigree))
  if (colnames(pedigree)[1]!=idname)
	stop("Error: The identifier for the pedigree does not match the unique identifier in the genotype file")

  Ainv <- NULL
  if (!missing(pedigree))
  {
	ped.Ainv <- asreml.Ainverse(pedigree)$ginv
	
	Ainv <- list()
	Ainv[[idname]] <- ped.Ainv
  } 

  object$nperm <- n.perm
  object$alpha <- alpha
  object$multtest <- multtest

  if (missing(baseModel))
  object$envModel <- list(fixed=fixed.forma, random=random, sparse=sparse, rcov=rcov, ginverse=Ainv)
  else object$envModel <- as.list(baseModel$call)[2:length(baseModel$call)]

  nphe <- object$nphe
  object$dfMerged[[idname]] <- as.factor(object$dfMerged[[idname]])

  n.mrk <- unlist(lapply(object$map, length))
  dimdM <- ncol(object$dfMerged)

   if (length(which(n.mrk==1))>0)
	message("Warning: Linkage groups ", which(n.mrk==1), " contain only one marker")

   object$dfMerged[, (1+nphe):dimdM] <- t(t(object$dfMerged[,(1+nphe):dimdM]) - apply(object$dfMerged[,(1+nphe):dimdM], 2, function(x) return(mean(x,na.rm=TRUE))))

  if (n.perm>0)
   object$dfMrk[, 2:ncol(object$dfMrk)] <- t(t(object$dfMrk[,2:ncol(object$dfMrk)]) - apply(object$dfMrk[,2:ncol(object$dfMrk)], 2, function(x) return(mean(x,na.rm=TRUE))))

  message("Beginning detection stage...")
  # note that which function is used will depend on the type of cross
  # e.g. 2 alleles, 3 alleles, or association
  det.out <- detect.asreml(object, filestem=filestem, ...)
  detqtl <- vector(length=length(object$map)) 

  for (i in 1:length(object$map))
	detqtl[i] <- length(grep(paste("C", i, "M", sep=""), det.out$qtls$pos))

  if (type=="f2") detqtl <- detqtl/2
   
  message("Detection stage complete. ", sum(detqtl), " total QTL detected")

  loc.out <- det.out
  if ((sum(detqtl)>0) & (object$loc==TRUE))
  {
    message("Beginning localization stage...")
    loc.out <- localize.asreml(object, QTLperChr=detqtl, ...)
  }

  message("Fitting final model with all QTL")

#  loc.out$qtls$pos <- sort(loc.out$qtls$pos)
  envModel <- object$envModel
  final.fixed <- paste(as.character(envModel$fixed[2]), "~", as.character(envModel$fixed[3]), sep="")
  if (sum(detqtl)>0)
     final.fixed <- paste(final.fixed, "+",paste(loc.out$qtls$pos, collapse="+"), sep="")

  final.fixed <- as.formula(final.fixed)
  final <- list(fixed=final.fixed, random=envModel$random, sparse=envModel$sparse, rcov=envModel$rcov, ginverse=Ainv, data=object$dfMerged, ...)
  final <- final[!sapply(final, is.null)]
  mod.fin <- do.call("asreml", final)

   output <- list()
   output$input <- object 
   output$no.qtl <- detqtl
   output$final.model <- mod.fin
 
   if (sum(detqtl)>0)
   {
     table <- list()

     table[[1]] <- rep(names(object$map), detqtl)
     if (type=="f2") table[[1]] <- rep(table[[1]], each=2)

     mappos <- unlist(object$mapp)
     if (type=="f2")
	names(mappos) <- names(object$dfMerged)[(nphe+1):ncol(object$dfMerged)][seq(1, ncol(object$dfMerged)-nphe, 2)] else 
     names(mappos) <- names(object$dfMerged)[(nphe+1):ncol(object$dfMerged)]
     tmpord <- match(loc.out$qtls$pos, names(mappos))
 
     loc.out$qtls$pos <- loc.out$qtls$pos[order(tmpord)]

     if (type=="f2") pos <- unique(substr(loc.out$qtls$pos, 1, nchar(loc.out$qtls$pos)-1)) else 
     pos <- loc.out$qtls$pos
     if (type=="f2") 	posD <- paste(pos, "D", sep="") else posD <- pos

     nmmap <- unlist(lapply(object$mapp, names))
     if (type!="other")
     table[[2]] <- round(mappos[sort(tmpord)], 2)
     if (type=="f2") table[[2]] <- rep(table[[2]], each=2)

     table[[3]] <- loc.out$qtls$pos
#     table[[3]] <- table[[4]] <- NULL
     table[[4]] <- NULL

     table[[5]] <- round(mod.fin$coefficients$fixed[match(loc.out$qtls$pos, names(mod.fin$coefficients$fixed))], 3)
     table[[6]] <- round((mod.fin$vcoeff$fixed[match(loc.out$qtls$pos, names(mod.fin$coefficients$fixed))]*mod.fin$sigma2)^.5, 3)

     table[[7]] <- round(table[[5]]/table[[6]], 2)
     table[[8]] <- round(sapply(table[[7]], function(x) 2*(1-pnorm(abs(x)))), 4)
     # Wald profile for each chromosome containing QTL
     output$profile <- loc.out$profile

     # how to deal with flanking markers? still output for anything but assoc
     if (type != "other")
     {
	mark <- grep("M", names(object$dfMerged)[(nphe+1):dimdM])+nphe
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
   if (is.null(table$Pos)) names(table)[3] <- "Marker"
   table <- table[!sapply(table, is.null)]

   output$Summary <- as.data.frame(do.call("cbind", table))
   rownames(output$Summary) <- NULL

   } # End of check for detected QTL
  
   if (sum(detqtl)==0)
 	message("No QTL detected in data. See log file for more details of testing.")
 
   class(output) <- c("dlmap", class(object))

   output
}

