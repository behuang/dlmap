`map.loc.asreml` <-
function(input, s.chr, chrSet, prevLoc=NULL, ...)
{
  dfMrk <- input$dfMrk
  dfMerged <- input$dfMerged
  envModel <- input$envModel
  nphe <- input$nphe  
  map <- input$mapp[[s.chr]]
  mrk <- grep(paste("C", s.chr, "M", sep=""), names(dfMerged))
  chr <- sort(c(mrk, grep(paste("C", s.chr, "P", sep=""), names(dfMerged))))
  type <- attr(input, "type")
  n.chr <- length(input$map)

  f.pos <- vector()
  f.mrk <- vector()
  if (length(prevLoc)>0) {
    f.pos <- prevLoc$pos
    f.mrk <- prevLoc$mrk
 
    fp <- match(f.pos, names(input$dfMerged))
    fm <- match(f.mrk, names(input$dfMerged))

    if (type=="f2") {
      	f.pos <- paste(unique(substr(f.pos, 1, nchar(f.pos)-1)), "D", sep="")
      	f.mrk <- unique(substr(f.mrk, 1, nchar(f.mrk)-1))
    }
  }

  formula <- envModel

  ### redefine dfMerged
  ## pull out the phenotypic data and fixed effect markers/pos
  if (ncol(input$dfMrk) > n.chr*nrow(dfMrk)) {
   dfm1 <- input$dfMerged[,c(1:nphe, fp, fm)]
   ## we're going to merge everything else onto this. 
   ## create separate groups of indices
   index <- list()
   mat <- list()
   for (kk in 1:n.chr)
   {
    index[[kk]] <- setdiff(grep(paste("C", kk, "M", sep=""),colnames(input$dfMrk)[2:ncol(input$dfMrk)]), match(prevLoc$mrk, colnames(input$dfMrk)[2:ncol(input$dfMrk)])) + 1
    mat[[kk]] <- input$dfMrk[,index[[kk]]]
    mat[[kk]] <- mroot(mat[[kk]] %*% t(mat[[kk]]))
    ncolm[kk] <- ncol(mat[[kk]])
   }
   cumind <- c(0, cumsum(ncolm))
  
   dfm2 <- as.data.frame(do.call("cbind", mat))
   dfm2 <- cbind(input$dfMrk[,1], dfm2)

   dfMerged2 <- merge(dfm1, dfm2, by=names(dfm2)[1], all.x=TRUE, sort=FALSE)

   for (kk in 1:n.chr) 
   formula$group[[paste("g_", kk, "chr", sep="")]] <- ncol(dfm1) + (cumind[kk]+1):cumind[kk+1]
  } else {
   for (kk in 1:n.chr) 
     formula$group[[paste("g_", kk, "chr", sep="")]] <- setdiff(grep(paste("C", kk, "M", sep=""),colnames(dfMerged)[(nphe+1):ncol(dfMerged)]), match(c(prevLoc$mrk, prevLoc$pos), colnames(dfMerged)[(nphe+1):ncol(dfMerged)])) + nphe
   dfMerged2 <- dfMerged
  }
  groups <- formula$group

 ##################################
  if (type=="f2") mrk <- mrk[seq(1, length(mrk), 2)]

  results <- list()
  wald <- rep(0, length(map))

  # Initialize convergence flag for output
  results$converge <- TRUE

  int <- vector()
  if (length(f.pos)>0) {
    fmrk <- sapply(match(f.pos, names(dfMerged)), function(x) {
	if (x==min(mrk)) return(c(min(mrk), min(mrk[mrk>x]))) else if (x==max(mrk)) return(c(max(mrk[mrk<x]), max(mrk))) else return(c(max(mrk[mrk<x]), min(mrk[mrk>x])))})

  # because we want to exclude both additive and dominant effects
  if (type=="f2") fmrk[2,] <- fmrk[2,]+1

  int <- eval(parse(text=paste("c(", paste(apply(fmrk, 2, function(x) return(paste(x[1], ":", x[2], sep=""))), collapse=","), ")", sep="")))
  }

  for (jj in 1:length(map))
  {
	if (type=="f2") pos <- c(2*jj-1, 2*jj) else pos <- jj
	
	# Check that position is not in the same interval as any previously mapped QTL
	#  need to check up on write up to make sure which fixed model elements are being fit on a given iteration - are we including fixed effects for QTL on other chromosomes? 

      if (all(!(chr[pos] %in% int)))
      {
	### create data frame with position and transformed random effects
	dfMerged3 <- dfMerged2
	if (length(intersect(names(dfMerged)[chr[pos]], names(dfMerged2)))==0){
	dfMerged3 <- cbind(dfMerged[, chr[pos]], dfMerged2)
	for (kk in 1:n.chr) formula$group[[paste("g_",kk,"chr",sep="")]] <- groups[[paste("g_",kk,"chr",sep="")]]+1 
	names(dfMerged3)[1] <- names(dfMerged)[chr[pos]] }

	chrnam <- paste("idv(grp(g_", setdiff(chrSet,s.chr), "chr))", sep="")
	formula$random <- paste("~", paste(chrnam, collapse="+"))

	# Include spatial/environmental random effects
	if (!is.null(envModel$random))
	  formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")

	formula$random <- as.formula(formula$random)

	formula$fixed <- paste(as.character(envModel$fixed)[2], "~", as.character(envModel$fixed[3]), sep="")

	if (length(f.pos) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(prevLoc$pos, collapse="+"), sep="")

	# not going to work for f2, need to correct this. 
	if (type=="f2")	
	formula$fixed <- paste(formula$fixed, "+", paste(paste(names(dfMerged)[chr[pos]], collapse="+"), sep="")) else formula$fixed <- paste(formula$fixed, "+", names(dfMerged)[chr[pos]], sep="")

	formula$fixed <- as.formula(formula$fixed)

	formula$data <- dfMerged3
	formula$Cfixed <- TRUE
	formula <- c(formula, ...)
	formula <- formula[!duplicated(formula)]
	formula <- formula[!sapply(formula, is.null)]

	if (length(chrSet)>1) model <- do.call("asreml", formula)

	if (length(chrSet)==1)
	{
	formula1 <- formula
	formula1$random <- envModel$random
	formula1 <- formula1[!sapply(formula1, is.null)]
	model <- do.call("asreml", formula1)
	}

	if (model$converge==FALSE) 	results$converge <- FALSE

	if (model$coefficients$fixed[names(model$coefficients$fixed) %in% names(dfMerged)[chr[pos]]] != 0)
	wald[jj] <- wald.test.asreml(model, list(list(which(model$coefficients$fixed[names(model$coefficients$fixed) %in% names(dfMerged3)[chr[pos]]]!=0), "zero")))$zres$zwald
      } # end of check for distinct intervals
  }

  results$wald <- wald
  return(results)
}

