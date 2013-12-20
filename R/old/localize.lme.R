`localize.lme` <-
function(input, QTLperChr)
{
  locations <- list()
  results <- list()

  type <- attr(input, "type") 
  dfMerged <- input$dfMerged
  chrSet <- which(QTLperChr>0)

    # Loop over all chromsomes which were significant in the preceding detection step
    for (kk in 1:length(chrSet))
    {
	# Count the number of QTLs which have already been mapped on the selected chromosome
     	no.qtls <- QTLperChr[chrSet[kk]][[1]]
	
	mrk <- grep(paste("C", chrSet[kk], "M", sep=""), names(dfMerged))
	chr <- sort(c(mrk, grep(paste("C", chrSet[kk], "P", sep=""), names(dfMerged))))

	if (type=="f2") {
	  mrk <- mrk[seq(1, length(mrk), 2)]
	  chr <- chr[seq(1, length(chr), 2)]
	}

	loc1 <- list()
	# 1D scan
	for (mm in 1:no.qtls)
	{
	  # 1D scan for QTL
     	  map.results <- map.loc.lme(input, s.chr=chrSet[kk], chrSet=chrSet, prevLoc=loc1)
	
	  # Find position of maximum wald score and marker markers
	  # or alternately use the likelihood for localization
	  sel.pos <- which.max(map.results$wald)

	  if (sel.pos==1) mark.l <- mrk[1] else 
		mark.l <- max(mrk[mrk<chr[sel.pos]])
	  if (sel.pos==length(chr)) mark.r <- mrk[length(mrk)] else
	  	mark.r <- min(mrk[mrk>chr[sel.pos]])

	  if (type=="f2")
	    loc1$mrk <- c(loc1$mrk, names(dfMerged)[c(mark.l:(mark.l+1), mark.r:(mark.r+1))]) else
	  loc1$mrk <- c(loc1$mrk, names(dfMerged)[c(mark.l, mark.r)])

	  if (type=="f2")
	    loc1$pos <- c(loc1$pos, names(dfMerged)[chr[sel.pos]:(chr[sel.pos]+1)]) else
	  loc1$pos <- c(loc1$pos, names(dfMerged)[chr[sel.pos]])

	  if (mm==1)
	  {
		results$profile[[names(input$map)[chrSet[kk]]]] <- t(rbind(input$mapp[[chrSet[kk]]], map.results$wald))
		colnames(results$profile[[names(input$map)[chrSet[kk]]]]) <- c("Position", "Wald")
	  }
  	}

	locations$mrk <- c(locations$mrk, loc1$mrk)
 	locations$pos <- c(locations$pos, loc1$pos)

    } # end of loop over significant chromosomes
 
  results$qtls <- locations
  return(results)
}

