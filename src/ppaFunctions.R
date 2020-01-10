
getBFFromPval = function(df) {
  df[,'Z'] <- pToZStat(df$P_Val)
  apply(df[,c('Z','F','N')], 1, function(d) calcLogBF_z(d['Z'],d['F'],d['N']))
}

pToZStat = function(pVals) {
  sqrt(qchisq(p=pVals, df=1, lower.tail=F))
}

calcLogBF_z = function(Z, f, N) {
  WW <- 0.1
  V <- approx_v(f, N)
  r <- WW / (V + WW)
  toreturn <- log( sqrt(1-r) ) + (Z*Z*r / 2)
  toreturn
}

approx_v = function(f, N) {
  1 / (2*f*(1-f) * N)
}

getPPAsFromPval = function(df, segmentCol, pCol = "P_Val", defaultMAF = 0.2, defaultN = 10000)
{
  df$P_Val = df[, pCol, drop=T]
  if (!has_name(df, "logBF")) {
    if (!has_name(df, "F")) {
      df$F = defaultMAF
      warning(sprintf("No column 'F' in dataframe, using default MAF of %f", defaultMAF))
    }
    df$F[is.na(df$F)] <- defaultMAF
    if (!has_name(df, "N")) {
      df$N = defaultN
      warning(sprintf("No column 'N' in dataframe, using default N of %d", defaultN))
    }
    df$N[is.na(df$N)] <- defaultN
    df$logBF = getBFFromPval(df)
  }
  getAllSegmentPPAs(df, segmentCol)
}

calcLogBF_beta = function(beta, se) {
  W <- 0.1
  V = se^2
  Z2 = beta^2 / V
  log( sqrt((V + W) / V) ) - (Z2 * W / (2 * (V + W)))
}

getBFFromBetaSE = function(df, betaCol, seCol) {
  #calcLogBF_beta gives a negative result for evidence against the null hypothesis,
  # but this is the opposite sign from what we want later.
  -apply(df[,c(betaCol,seCol)], 1, function(d) calcLogBF_beta(d[1],d[2]))
}

getPPAsFromBetaSE = function(df, segmentCol, betaCol, seCol)
{
  df$logBF = getBFFromBetaSE(df, betaCol, seCol)
  getAllSegmentPPAs(df, segmentCol)
}

getAllSegmentPPAs = function(df, segmentCol) {
  segmentVals = df[, segmentCol, drop=T]
  segIndices <- getUniqueIndices(segmentVals)
  if (length(segIndices) > length(unique(segmentVals))) {
    warning("getPPAs: ERROR - data.frame must be ordered by the segmentCol column.")
    ppas = NULL
  } else {
    ppas <- vector(mode="double")
    # For each gene/segment...
    for (i in 1:length(segIndices)) {
      startIndex <- segIndices[[i]][[2]]
      endIndex <- segIndices[[i]][[3]]
      ppas[startIndex:endIndex] <- getNaivePPA(df$logBF[startIndex:endIndex])
    }
  }
  ppas
}

getNaivePPA = function(vecLogBF)
{
  logsegbfNaive <- -1000
  for (i in 1:length(vecLogBF)) {
    if (!is.na(vecLogBF[i])) {
      logsegbfNaive <- sumlog(logsegbfNaive, vecLogBF[i])
    }
  }
  vecPPA <- vecLogBF - logsegbfNaive
  exp(vecPPA)
}

sumlog = function(logx, logy)
{
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}

# First determine the indices of the genes in the table
getUniqueIndices = function(vec)
{
  if (class(vec) == "factor") {
    # For some reason factors are INCREDIBLY slow if used in the code below
    vec <- as.character(vec)
  }
  indicesList <- list()
  if (length(vec) == 0) {
    return(indicesList)
  }
  
  lastVal = vec[1]
  lastIndex = 1
  for (i in 1:length(vec)) {
    if (vec[i] != lastVal) {
      indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i-1))))
      lastVal <- vec[i]
      lastIndex <- i
    }
  }
  indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i))))
  return(indicesList)
}
