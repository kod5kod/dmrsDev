methylSmoother<- function (BSseq, minM = 70, minB = 1000, sampIndx, indexes){
    ########################
    ##This function smoothes the samples methylation along each sample of BSseq 
    ##Args:
    ##  BSseq: BSseq DMRs object
    ##  minM: A smoothed BSseq object 
    ##  minB: a vector specifying the sample indexes of group 1 of the 
    ##        BSseq object (e.g. c(1,3,5,7,9,11) )
    ##  sampIndx: a vector specifying the sample indexes of group 2 
    ##            (e.g. c(2,4,6,8,10,12) )
    ##  indexes: e.g. 1:5 or 10:20
    ##Returns:
    ##  A 
    ##Requires: library('GenomicRanges'), library('bsseq')
    ## BASED ON BSseq - BSsmooth
    ########################
    ##Getting the coverage readouts per  specific sample, over  specific indexes
    Cov <- unname(as.array(getCoverage(BSseq, type = "Cov")[indexes, sampIndx])) 
    ##Getting the methylation readouts per specific sample, overspecific indexes
    M <- unname(as.array(getCoverage(BSseq, type = "M")[indexes, sampIndx])) 
    pos <- start(BSseq)[indexes] # getting the position of all hte readouts
    ##Add checks that positives is ___stopifnot(all(diff(pos) > 0))
    positives <- which(Cov != 0) # getting the non-zero methylation loci
    winRatio <- minM / length(positives) # locfit lp parameter 
    if (length(positives) <= minM){
        cat("  The number of positive coverage loci is: ", length(positives),"
            , while the minimum methylation window is: ",minM)
        stop("the number of positive coverage loci is lower than the minimum 
              methulation loci windows. Please either decrease teh min 
              methylation loci, or increase the input.\n")
    }
    ##Creating a dataframe with the relevant positions and the methylation, 
    ## making sure there that when M=Cov we take M or Cov - 0.01 to avoid issues
    methylData <- data.frame(pos = pos[positives], 
                             M = pmin(pmax(M[positives], 0.01), 
                                      Cov[positives] - 0.01), 
                             Cov = Cov[positives])
    # Applying locfit with local polynomial term
    fit <- locfit(M ~ lp(pos, nn = winRatio, h = minB),
                  data = methylData, weights = Cov, family = "binomial", 
                  maxk = 10000)
    pp <- preplot(fit, where = "data", band = "local",
                  newdata = data.frame(pos = pos)) 
    se.coef <- pp$se.fit

    plot(pos[positives], pp$coef[positives], type = "p", cex = 0.1)
    lines(pos[positives], M[positives], type = "p", col = "red", cex = 0.1)
    return(list(coef = pp$fit, se.coef = se.coef,
                trans = pp$trans, minB = minB, minM = minM))

}
