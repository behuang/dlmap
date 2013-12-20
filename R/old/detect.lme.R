`detect.lme` <-
function(input, filestem)
{
  chr.names <- names(input$map)
  type <- attr(input, "type")
  n.chr <- length(input$map)

  locations <- list()
  results <- list()

  logfile = paste(filestem, ".det.log", sep="")

  write("Note: p-values and threshold are analytical", logfile)
  write("****************************************************", logfile)

  # Fit initial detection step with all chromosomes
  # Output contains p-values for each chromosome variance component
  test <- test.lme(input, chrSet=1:n.chr)

  # Output observed statistics to logfile
  write("Iteration 1: ", logfile, append=TRUE)
  write(paste(c("",chr.names), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh,4), sep=""), logfile, append=TRUE)

  # Flag for significant detection step - one stopping rule is when this is 0
  found <- (min(test$adj.pval) < input$alpha)

  # Subset of significant chromosomes - only use these in localization steps
  chrSet <- which(test$adj.pval < input$alpha)

  write("Significant chromosomes to be used for scanning/testing:", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)

  # Keep track of iteration number
  itnum = 1

  while (found==1) 
  {
    itnum=itnum+1

    # list of locations from previous iteration - keeps track of markers
    # which have already been fitted on each chromosome
    loc.temp <- locations

    # Loop over all chromsomes which were significant in the preceding step
    # Scan each one using a 1/2/3D scan, then insert all mapped QTLs in model as fixed effects
    markers <- vector(length=length(chrSet))
    for (jj in 1:length(chrSet))
    {
	# Count the number of QTLs which have already been mapped on the selected chromosome
     	no.qtls <- length(grep(paste("C", chrSet[jj], "M", sep=""), locations))
	if (type=="f2") no.qtls <- no.qtls/2

     	map.results <- map.det.lme(input, chrSet[jj], chrSet, loc.temp)
	
	# Find position of maximum wald score and marker markers
	sel.mrk <- which.max(map.results$wald)

	markers[jj] <- sel.mrk

	if (type=="f2")
	 locations <- c(locations, paste("C", chrSet[jj], "M", sel.mrk, c("D", "A"), sep="")) else
	 locations <- c(locations, paste("C", chrSet[jj], "M", sel.mrk, sep=""))

	if (itnum==2)
	{
 	  results$profile[[names(input$map)[chrSet[jj]]]] <- rbind(as.vector(input$map[[chrSet[jj]]]), map.results$wald)
	  rownames(results$profile[[names(input$map)[chrSet[jj]]]]) <- c("Position", "Wald")
	}
    } # end of loop over significant chromosomes

    locations <- unlist(locations)

    # Output selected markers to logfile
    write(paste(c("Mrk:", markers), collapse="\t"), logfile, append=TRUE)

    # Do not allow more than one QTL per interval on a chromosome
    ex <- vector()
    chrSet.old <- chrSet

    for (ii in 1:length(chrSet))
 	if ((no.qtls+2)>=length(input$map[[chrSet[ii]]]))
	  ex <- c(ex, ii)

    chrSet <- setdiff(chrSet, chrSet[ex])
  
    # Make note in logfile if chromosome removed for this reason
    if (length(ex)>0)
    {
	write("*******************************************************", logfile, append=TRUE)
	write("Chromosomes full - One QTL located per interval already", logfile, append=TRUE)
	write(paste("Chromosome(s): ", paste(c("", chr.names[chrSet.old[ex]]), collapse="\t"), " removed from testing set", sep=""), logfile, append=TRUE)
	write("*******************************************************", logfile, append=TRUE)
    }

    found=0
    # now reperform detection - rescan over all chromosomes for detection
    if (length(chrSet)>0)
    {
	chrSet.old <- chrSet
	# Test previous subset of chromosomes for significance 
  	test <- test.lme(input, chrSet, prevLoc=locations)

	# Check whether QTLs are detected on any chromosomes
  	found <- (min(test$adj.pval) < input$alpha)

	# Find the subset of chromosomes which is significant
  	chrSet <- chrSet[which(test$adj.pval < input$alpha)]
    }


  # Output other iterations to logfile
  write("***************************************************************", logfile, append=TRUE)
  write(paste("Iteration ", itnum,":", sep=""), logfile, append=TRUE)

  write("Chromosomes from previous iteration: ", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet.old]), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh,4), sep=""), logfile, append=TRUE)

  if (length(chrSet)>0)
  {
   write("Significant chromosomes for next round of testing/scanning:", logfile, append=TRUE)
   write(paste(c("",chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)
  }
 }  # End of iterative search

  results$qtls$pos <- locations
  return(results)
}

