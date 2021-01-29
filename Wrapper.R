
#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("tidyr", "dplyr", "plyr", "reshape2", "seqinr", "ggplot2", "coefplot", 
             "forcats", "tibble", "stringr", "purrr", "gridExtra", "pracma", "hablar")  

invisible(lapply(packages, library, character.only = TRUE)) # add imported packages to library

#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# START CODE FOR RUNNING IN RSTUDIO (comment out if running from TurnoveR)
#------------------------------------------------------------------------------------


filepath <<- "C:/Users/alimarsh/Documents/Turnover/Data/report.csv"
tool.dir <<- "C:\\Users\\alimarsh\\Documents\\Turnover\\ToolHolder"
diet.enrichment <- as.numeric ("0.999999") # Leucine percent enrichment in diet
min.avg.turnover.score <<- as.numeric ("0")
min.isotope.dot.product <<- as.numeric ("0")
folder.name <- "Data"
Detection.Qvalue.threshold <- as.numeric ("1")

setwd("C:/Users/alimarsh/Documents/Turnover/Data")

#------------------------------------------------------------------------------------
# END CODE FOR RUNNING IN RSTUDIO
#------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------
# START CODE FOR RUNNING FROM TURNOVER (comment out if running from RSTUDIO)
#------------------------------------------------------------------------------------
# 
# #------------------------------------------------------------------------------------
# # LOAD ARGUMENTS FROM SKYLINE #
# 
# arguments <- commandArgs(trailingOnly=TRUE)
# cat(length (arguments))
# if ( length (arguments) != 7)
#   # expected arguments not present -- error
#   stop ("USAGE: R --slave --no-save --args '<textbox><textbox><textbox><textbox>'")
# for (i in 1:7) {
#   arg <- arguments [i]
#   # remove leading and trailing blanks
#   arg <- gsub ("^ *", "", arg)
#   arg <- gsub (" *$", "", arg)
#   # remove any embedded quotation marks
#   arg <- gsub ("['\'\"]", "", arg)
#   #report file is brought in as an argument, this is specified in TestArgsCollector.properties
#   
#   #TODO put [arg] back for all
#   if (i==1) filepath <<- arg
#   if (i==2) tool.dir <<- arg
#   if (i==3) diet.enrichment <- as.numeric (arg) # Leucine percent enrichment in diet
#   if (i==4) min.avg.turnover.score <<- as.numeric (arg)
#   if (i==5) min.isotope.dot.product <<- as.numeric (arg)
#   if (i==6) folder.name <- arg
#   if (i==7) Detection.Qvalue.threshold <- as.numeric (arg)
# }
# 
# dir.create(file.path(getwd(), folder.name), showWarnings = FALSE) # Create folder for script output
# setwd(file.path(getwd(), folder.name))
# 


#------------------------------------------------------------------------------------
# END CODE FOR RUNNING FROM TURNOVER
#------------------------------------------------------------------------------------


diet.enrichment <- ifelse(diet.enrichment==1, 0.999999, diet.enrichment)

#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------
# SET CONSTANT ARGUMENTS #
min.abundance <- 0.0001 #as.numeric (arg) # minimum abundance
resolution <- 0.1 #as.numeric (arg) # resolution for distinguishing peaks
p.tolerance <- 0.05 #as.numeric (arg) # tolerance for combining masses in observed data

#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------
# LOAD DATA #

df.input <- read.csv(file = filepath, stringsAsFactors = F)


if (is.null(df.input$Replicate.Name)) {
  # rename df columns with periods between words
  df.input <- rename_with(df.input, .fn = function(vector){
    return(c("Protein", "Replicate.Name", "Protein.Description", "Protein.Accession", "Protein.Gene", "Peptide", 
             "File.Name", "Timepoint", "Condition", "Precursor.Charge", "Precursor.Mz", "Molecule.Formula", 
             "Precursor.Neutral.Mass", "Modified.Sequence", "Is.Decoy", "Detection.Q.Value", "Total.Area.Ms1", 
             "Isotope.Dot.Product", "Product.Mz", "Product.Charge", "Fragment.Ion", "Isotope.Dist.Index", 
             "Isotope.Dist.Rank", "Isotope.Dist.Proportion", "Fragment.Ion.Type", "Area"))
    
  })
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# RUN STEP 1 #
source(paste(tool.dir, "Step1.R", sep="/"))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# RUN STEP 2 #
# single leucine data set (1 leucine)
data.s.holder <- read.csv(paste(getwd(),"/Step1_Data_Output_Skyline_singleleucine_peps_date.csv", sep=""), stringsAsFactors = F)

# multiple leucine data set (2,3,4 leucines)
data.m.holder <- read.csv(paste(getwd(),"/Step1_Data_Output_Skyline_multileucine_peps_date.csv", sep=""), stringsAsFactors = F)
data.s <- data.s.holder
data.m <- data.m.holder
source(paste(tool.dir, "Step2.R", sep="/"))

#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------
# RUN STEP 3 #
data.s <- data.s.holder
data.m <- data.m.holder
# medians of x-intercepts by cohort from step 3
df.x.int.medians <- read.csv(paste(getwd(),"/Table_step2_xintercepts.csv", sep=""), stringsAsFactors = F)

source(paste(tool.dir, "Step3.R", sep="/"))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# RUN STEP 4 #
data.s <- data.s.holder
data.m <- data.m.holder
# medians of x-intercepts by cohort from step 3
df.x.int.medians <- read.csv(paste(getwd(),"/Table_step2_xintercepts.csv", sep=""), stringsAsFactors = F)

source(paste(tool.dir, "Step4.R", sep="/"))
#------------------------------------------------------------------------------------



cat("\n---------------------------------------------------------------------------------------")
cat(" ALL COMPLETED ")
cat("---------------------------------------------------------------------------------------\n\n")
cat("Output at: ")
cat(getwd())






