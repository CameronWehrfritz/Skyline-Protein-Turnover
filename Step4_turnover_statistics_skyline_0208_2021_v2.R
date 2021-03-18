#Written by Cameron Wehrfritz
#and Natan Basisty, PhD
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#March, 2020
#updated: March 13, 2021

# PROTEIN TURNOVER ANALYSIS
# STEP 4:
# Linear model of log(Percent.Newly.Synthesized) by timepoints and conditions and their interaction
#
# OUTPUT:
# i. Data table of statistics

###########################
### Begin Script Step 4 ###
###########################

#------------------------------------------------------------------------------------
# START CODE FOR RUNNING IN RSTUDIO (comment out if running from TurnoveR)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# set working directory
setwd("//bigrock/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #
packages = c("tidyr", "dplyr", "reshape2", "seqinr", "ggplot2", "coefplot", "plyr")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LOAD DATA

# single leucine data set (1 leucine)
data.s <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts/Additional_Output_Data/Step1_Data_Output_Skyline_singleleucine_peps.csv", stringsAsFactors = F) # PC

# multiple leucine data set (2,3,4 leucines)
data.m <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts/Additional_Output_Data/Step1_Data_Output_Skyline_multileucine_peps.csv", stringsAsFactors = F) # PC

# medians of x-intercepts by condition from Step 2 
df.x.int.medians <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0814_Skyline_Turnover_Tool/Turnover_R_scripts/Additional_Output_Data/Turnover_step2_xintercepts.csv", stringsAsFactors = F) # PC
#------------------------------------------------------------------------------------

# Reference condition
Reference.Condition <- "OCon" # this should be assigned by the user

#------------------------------------------------------------------------------------
# END CODE FOR RUNNING IN RSTUDIO
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Combine Single Leucine and Multiple Leucine data sets together for modeling

df <- data.m %>%
  bind_rows(data.s) %>% # combine data; retains all columns, fills missing columns in with NA
  filter(Perc.New.Synth>0) # filter out data with negative percent newly synthesized value
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PREP FOR MODEL 

# conditions
conditions <- unique(df$Condition)

# make a conditions vector for looping through all comparisons
conditions.loop <- conditions[!conditions==Reference.Condition] # keep all conditions except for the reference conditions

# proteins
prots <- unique(df$Protein.Accession)

# genes
genes <- unique(df$Protein.Gene)

# time points
time <- sort(unique(df$Timepoint))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# create modified time 
# by subtracting the median x-intercept time (shifting left toward the origin)
# unless x-intercepts are negative, then modified.time is simply the same as time
for(i in conditions){
  if(all(df.x.int.medians %>% pull(Median.x.intercept)>0, df.x.int.medians %>% pull(Median.x.intercept)<min(time))){ # check if all median x-intercepts are positive and less than minimum timepoint
    df$Modified.Time[df$Condition==i] <- df %>% filter(Condition==i) %>% pull(Timepoint) - df.x.int.medians %>% filter(Condition==i) %>% pull(Median.x.intercept) # modify timepoints by translating left by respective median x-intercept
  } else { 
    df$Modified.Time[df$Condition==i] <- df %>% filter(Condition==i) %>% pull(Timepoint) # else do not modify
  }
} # end for
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LINEAR MODELLING

# initialize data frame to write out results from linear model 
# first figure out column names since the number of columns is built into col size
col.names <- c("Protein.Accession", "Gene", "Comparison", "No.Peptides", "No.Points", "Interaction" , "Std.Error",  "t.value", "Unadj.P", "Qvalue", "DF",
               "Slope.Numerator", "Slope.Denominator", "Half.Life.Numerator", "Half.Life.Denominator")
df.model.output <- data.frame(matrix(nrow = length(conditions.loop)*length(prots), ncol = length(col.names)))
names(df.model.output) <- col.names

# Loop
row.index <- 1 # counter
for(i in prots){
  
  # subset data for protein 
  data.protein.loop <- subset(df, Protein.Accession == i) 
  
  # subset data from reference condition - do this prior to the conditions loop below
  data.ref <- subset(data.protein.loop, Condition==Reference.Condition) # ref refers to the user defined reference condition
  
  # loop through conditions - in order to generate each comparison to the reference condition
  for(j in conditions.loop){
    
    # subset data for variable condition and reference condition; for use in combined linear model
    data.condition.loop <- subset(data.protein.loop, Condition == j | Condition == Reference.Condition )
    
    # subset data for variable condition
    data.var <- subset(data.protein.loop, Condition==j) # var refers to the variable condition, which will be compared against the user defined reference condition
    
    # create name of comparison
    comparison <- paste(j, "/", Reference.Condition, sep="") # variable condition vs reference condition
  
    # write out comparison name
    df.model.output[row.index, colnames(df.model.output)=="Comparison"] <- comparison
    
    # write out Protein.Accession
    df.model.output[row.index, colnames(df.model.output)=="Protein.Accession"] <- i 
    
    # write out Gene name
    df.model.output[row.index, colnames(df.model.output)=="Gene"] <- data.protein.loop %>% pull(Protein.Gene) %>% unique() 
    
    # write out number of unique peptides
    df.model.output[row.index, colnames(df.model.output)=="No.Peptides"] <- data.condition.loop %>% pull(Modified.Peptide.Seq) %>% unique() %>% length()
    
    # write out number of data points
    df.model.output[row.index, colnames(df.model.output)=="No.Points"] <- nrow(data.condition.loop)

    # LINEAR MODEL #
    if( nrow(data.condition.loop) >= 1.5*length(time) ){ # quick assessment: if the number of data points is greater than length of time points then the model should hopefully converge
      tryCatch(
        expr={ 
          # Model
          model <- lm(log(1-Perc.New.Synth) ~ 0 + Condition*Modified.Time, data = data.condition.loop %>% filter(Perc.New.Synth<1)) # this model matches the model used in step 3
          
          # write out statistics from combined model
          df.model.output[row.index, c(6:8)] <- summary(model)$coef[4, 1:3] %>% round(., digits=4) # model statistics: estimate, standard error, t value # fourth row should be the interaction term of condition and time variables
          df.model.output[row.index, 9] <- summary(model)$coef[4, 4] # p-value from combined linear model; not rounded, so we can sort by this variable # fourth row should be the interaction term of condition and time variables
          df.model.output[row.index, colnames(df.model.output)=="DF"] <- summary(model)$df[2] # degrees of freedom

          # model reference condition against its timepoints
          model.ref <- lm( log(1-Perc.New.Synth) ~ 0 + Modified.Time, data = data.ref %>% filter(Perc.New.Synth<1)) # this model matches the model used in step 3
          # model variable condition against its timepoints
          model.var <- lm( log(1-Perc.New.Synth) ~ 0 + Modified.Time, data = data.var %>% filter(Perc.New.Synth<1)) # this model matches the model used in step 3
          
          # write out slopes and half-lives from individual linear models
          # variable condition
          df.model.output[row.index, "Slope.Numerator"] <- -summary(model.var)$coef[1] %>% round(., digits=4) # slope of linear model
          df.model.output[row.index, "Half.Life.Numerator"] <- -log(2)/summary(model.var)$coef[1] %>% round(., digits=4) # half life
          # reference condition
          df.model.output[row.index, "Slope.Denominator"] <- -summary(model.ref)$coef[1] %>% round(., digits=4) # slope of linear model
          df.model.output[row.index, "Half.Life.Denominator"] <- -log(2)/summary(model.ref)$coef[1] %>% round(., digits=4) # half life
        },
        error=function(e){
          message("Caught an error!")
          print(e)
        },
        warning=function(w){
          message("Caught a warning!")
          print(w)
        },
        finally={
          message("All done, quitting.")
        }
      ) # end tryCatch
      
    } else{ # otherwise there may not be enough data points to run the model, then write out NA and continue looping
      df.model.output[row.index, c(6:9)] <- NA
      df.model.output[row.index, colnames(df.model.output)=="DF"] <- NA
    } # end else
    # increment row.index counter before iterating through condition loop
    row.index <- row.index + 1
  } # end for; condition level
} # end for; protein level
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# final calculations
df.model.output <- df.model.output %>%
  mutate(Log2.Ratio.Half.Life=log2(Half.Life.Numerator/Half.Life.Denominator)) %>% # calculate log2.ratio of half lives
  mutate(Adj.Pvalue=p.adjust(Unadj.P, method="BH")) # calculate adjusted pvalues using Benjamini & Hochberg (1995) "BH" method
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# clean model output
df.model.output <- df.model.output %>%
  select(-Qvalue) %>%  # Qvalue package is not working, remove Qvalue variable
  arrange(Unadj.P) # arrange best Pvalue top down
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write out model output
write.csv(df.model.output, file = "Half-lives_statistics.csv", row.names = FALSE)
#------------------------------------------------------------------------------------


# END 
