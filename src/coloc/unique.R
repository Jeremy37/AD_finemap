#!/usr/bin/env Rscript
# Use this script to remove duplicates in a table based on the values in one
# or more columns being duplicated. Use as follows, with 1-based column indexes:
#     unique.R --args file=input.txt cols=2,3 > output.txt
# or using column names:
#     unique.R --args file=input.txt header=T cols=2,3 > output.txt
# or using stdin:
#     cat file.txt | unique.R --args cols=2,3 > output.txt
#
myargs <- NULL

main = function()
{
  myargs <<- getArgs()
  #print(myargs)
  
  if (is.null(myargs$cols)) {
    stop("Missing parameter --cols")
  }
  if (is.integer(myargs$cols)) {
    cols <- myargs$cols
  } else {
    cols <- as.integer(unlist(strsplit(myargs$cols, split=",")))
  }
  
  sep <- "\t"
  #if (!is.null(myargs$sep)) {
  #  sep <- myargs$sep
  #}
  hasHeader <- F
  if (!is.null(myargs$header)) {
    hasHeader <- as.logical(myargs$header)
  }
  
  f <- NULL
  if (is.null(myargs$file)) {
    f <- file("stdin")
  } else {
    f <- file(myargs$file, open="r")
  }
  options(stringsAsFactors=F)
  df <- read.table(f, header=hasHeader, sep=sep, row.names=NULL, strip.white=F, quote="", fill=T, na.strings=c("NA"))
  df.subset <- df[!duplicated(df[,cols]),]
  write.table(df.subset, file="", sep=sep, quote=F, row.names=F, na="", col.names=hasHeader)
}



###########################################################################
##' commandArgs parsing
##' 
##' return a named list of command line arguments
##'
##' Usage:
##' call the R script thus
##'   ./myfile.R --args myarg=something
##' or
##'   R CMD BATCH --args myarg=something myfile.R
##'
##' Then in R do
##'   myargs <- getArgs()
##' and myargs will be a named list
##' > str(myargs)
##' List of 2
##' $ file : chr "myfile.R"
##' $ myarg: chr "something"
##'
##' @title getArgs
##' @param verbose print verbage to screen 
##' @param defaults a named list of defaults, optional
##' @return a named list
##' @author Chris Wallace
getArgs = function(verbose=FALSE, defaults=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}


###########################################################################

main()
