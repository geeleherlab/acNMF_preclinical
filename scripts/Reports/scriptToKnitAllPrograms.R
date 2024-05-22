library("knitr")

# create new versions of this script that input the new program number, spin them and create new output files...

for(i in 1:41) # numPrograms, can't put variable if I'm rm() and gc()'ing everything below.
{
    print(paste("Spinning Program ", i, sep=""))
    
    # Read in the original script....
    fileNameIn <- "gepSummaryKnitr.R"
    baseScript <- readChar(fileNameIn, file.info(fileNameIn)$size)
    
    # this will set the program number, which will create the report for that Gene expression program.
    thisProgSting <- paste("#+ setTheProgNumAutomatically, echo=FALSE, include=FALSE\nthisProgNumber <- ", i, "\n\n", sep="")
    
    # Append this line of code to our main script, write out a new script with a unique file name, then spin it.
    theNewScriptsText <- paste(thisProgSting, baseScript, sep="")
    
    # write the new file
    outFileName <- paste("gepSummaryKnitr_", i , ".R", sep="")
    fileConn <- file(outFileName)
    writeLines(theNewScriptsText, fileConn)
    close(fileConn)
    
    spin(outFileName)
    
    # remove all objects and do garbage collection to get back memory before we go at it agin.
    rm(list=ls())
    gc()
}

# This code will write the HTML links to the files that are created here
for(i in 1:41)
{
    filename <- paste0("gepSummaryKnitr_", i , ".html")
    cat(paste0("<br><a href='", filename, "'>Program ", i, "</a>\n"))
    
    #Fix the show/hide button
    #Replaces html code for single quotes (&#39;) to single quotes in order to restore functionality to the button
    x <- readLines(filename)
    y <- gsub("&#39;", "'", x)
    cat(y, file = filename, sep = "\n")
}


