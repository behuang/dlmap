`map.det.lme` <-
function(input, s.chr, chrSet, prevLoc=NULL)
{
  dfMerged <- input$dfMerged
  fixed <- input$envModel$fixed
  n.chr <- length(input$map)
  nphe <- input$nphe

  type <- attr(input, "type")
  results <- list()
  formula <- list()
  chrRE <- vector()
  wald <- rep(0, length(input$map[[s.chr]]))

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  chrRE[kk] <- paste("pdIdent(~", paste(setdiff(names(dfMerged)[grep(paste("C", kk, "M", sep=""), names(dfMerged))], prevLoc), collapse="+"), "-1)", sep="")

  # Loop over positions on the selected chromosome
  mrkloop <- unlist(lapply(strsplit(names(dfMerged)[setdiff(grep(paste("C", s.chr, "M", sep=""), names(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)]))+nphe], "M"), function(x) return(x[2])))

  if (type=="f2")
	mrkloop <- unique(substr(mrkloop, 1, nchar(mrkloop)-1))
  mrkloop <- as.numeric(mrkloop)

  for (jj in 1:length(mrkloop))
  {
    # Only include random effects for chromosomes not being scanned
    if (length(chrSet)>2)	
	formula$random <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet, s.chr)], collapse=","), "))", sep="")

    if (length(chrSet)==2)
	formula$random <- chrRE[setdiff(chrSet,s.chr)]

    formula$fixed <- paste(as.character(fixed)[2], "~", as.character(fixed)[3], sep="")

	# If there are markers already mapped on other chromosomes, add in as fixed effects 
    if (length(prevLoc) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(prevLoc, collapse="+"), sep="")

    	effectnames <- paste("C", s.chr, "M", mrkloop[jj], sep="")
    	if (type=="f2") effectnames <- paste(effectnames, c("A", "D"), sep="")
    	
	formula$fixed <- paste(formula$fixed, "+", paste(effectnames, collapse="+"), sep="")

	formula$fixgrp <- paste(formula$fixed, "| grp1", sep="")
	formula$fixed <- as.formula(formula$fixed)
	formula$fixgrp <- as.formula(formula$fixgrp)

	gd <- groupedData(formula$fixgrp, data=dfMerged)

	# Fit model - different forms depending on relevant terms
	if (length(chrSet)>1)
	model <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	if (length(chrSet)==1)
	model <- lme(fixed=formula$fixed, random=~1|grp1, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	# Get output - Wald statistic for position fixed effect
	wald[jj] <- anova(model, Terms=effectnames)[1,3]
  }
  results$wald <- wald

  return(results)
}

