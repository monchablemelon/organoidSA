options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()
pacman::p_load(pacman, gridExtra, reshape2,tidyverse, cluster, ggpubr, rstatix, factoextra, pROC, dplyr, GGally, glmnet, ggplot2, rjava, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, FSelector, rpart, caret, rpart.plot, xlsx, data.tree, caTools, car)

#####################################################################################
##  Before running this program, please make sure of the following:                ##
##  Remove all N/As, blanks, extraneous data, clean the table, and any weird fonts ##
##  Enter relevant data below.                                                     ##
##  Ensure excel sheet follows the format described in the README.                 ##
######################## Please set data content here!###############################

file_name                     <- "Example.xlsx"                               # < These settings can  
target_file_location          <- "~/Desktop/File1/File2"                      # < be ignored if manual_import = TRUE
manual_import                 = FALSE                                         # < default = FALSE, for automatic file upload. 
area_sheetname                <- "area"                                       # < area sheet name
dark_sheetname                <- "dark"                                       # < darkness sheet name
ecc_sheetname                 <- "ecc"                                        # < eccentricity sheet name 
x_variable_normalisation      <- "none"  #types = mm, zs, ,rs, fc, minfc, absfc (minmax, z-score, robust standardisation, fold-change, minimum-value fold-change, **sets minimum  ONLY FOR dark/time/SA)
# If roc is to be viewed (change_in_dark) view via none (dont use mm/zs/rs) - fc is susceptible to fluctuation from initial value as it is a singular x value, not a slope (slopes are less affected as slope gradient doesn't change). 
# DO NOT use x_variable_normalisation if possible. 
true_area_model               <- "TRUE" #If TRUE, model as ellipsoid, if FALSE, model as sphere).
graph_title                   <- "Graph Title Here"
y_axis_lab                    <- "AVG ΔDarkness/um^2*h"

light_abs_adjustment          <- "TRUE" #If true, adjusts for non-linear light absorption for darkness measurements
light_abs_k                   <- 0.52 #default light transmission (default = 0.52)
#### troubleshooting ######
abs_dark_troubleshooting     <- "FALSE" #(default = FALSE, TRUE to view cumulative darkness)
change_in_dark               <- "TRUE" #(default = TRUE, TRUE to view change in darkness/unit time / unit area. Overrides abs_dark_troubleshooting and Surfacegraphing)
Surfacegraphing              <- "FALSE" #(default = FALSE, TRUE to view graph of surface area instead of darkness. default false for darkness measurements)
#view(SA_graph) to view SA converted version of Area
#View(true_dark_transposed)) to view relevant data table
###########################################################################################
################# CAUTION: Do not edit content below unless necessary #####################
###########################################################################################
templist                      <- list("mm","zs","rs","mn","fc", "minfc")
#import conditions
if (manual_import == FALSE) {
  setwd(target_file_location)
  mydatarawarea <- read.xlsx(file_name, sheetName = area_sheetname)
  mydatarawdark <- read.xlsx(file_name, sheetName = dark_sheetname)
  mydatarawecc <- read.xlsx(file_name, sheetName = ecc_sheetname)}

#basic pre-processing function 
#Remove first column + replace elapsed with id

all_df <- list(mydatarawarea, mydatarawdark, mydatarawecc)
for (i in 1:length(all_df)) {
  # Remove the first two columns + 1 row
  all_df[[i]] <- all_df[[i]][-c(1,2) ]
  all_df[[i]] <- all_df[[i]][-c(1), ]
  # Rename column if "Elapsed" is found
  #all_df[[i]]$NA.[ all_df[[i]]$NA. == "Elapsed"] <- "id"
}

tall_dfs<- list()
for (i in 1:length(all_df)) {
  cat("Transposing DataFrame", i, ":\n")
  temp <- t(all_df[[i]])
  tall_dfs[[i]] <- temp
  tall_dfs[[i]] <- as.numeric(tall_dfs[[i]])
}

num_rows <- nrow(all_df[[1]])
num_cols <- ncol(all_df[[1]])

#Light absorption adjustment 
if (light_abs_adjustment == TRUE){
  tall_dfs[[2]] <- -log((1-tall_dfs[[2]]/100)/light_abs_k)
}


# Reshape the vector of doubles into a matrix
# Adjust byrow if needed
# Convert the matrix back to a dataframe
#true_SA <- 2*(tall_dfs[[1]])*(1+((1-(tall_dfs[[3]])^2)/(tall_dfs[[3]]))*atanh(tall_dfs[[3]]))
true_SA <- 2*pi*(sqrt((tall_dfs[[1]])/(pi*sqrt(1-(tall_dfs[[3]])^2))))^2*(1+((1-(tall_dfs[[3]])^2)/(tall_dfs[[3]]))*atanh(tall_dfs[[3]]))
#true_SA <- (2*(tall_dfs[[1]]))*(1+(asin(tall_dfs[[3]]))*(((tall_dfs[[1]]/pi))/(1-(tall_dfs[[3]])^2))^0.5/(((tall_dfs[[1]]/pi)^0.5)*(tall_dfs[[3]])))
### ^ prolate formula
#test to overestimate area 
#true_SA <- 2*(tall_dfs[[1]])/(1-(tall_dfs[[3]]))*(1+((1-(tall_dfs[[3]])^2)/(tall_dfs[[3]]))*atanh(tall_dfs[[3]]))
SA_graph <- true_SA
true_dark <- tall_dfs[[1]]*tall_dfs[[2]]/true_SA 

###used to switch to spherical model (i.e., dark=dark, eccentricity independent)
if (true_area_model == FALSE){
  true_dark <- tall_dfs[[2]]
}
if (true_area_model == FALSE & Surfacegraphing == TRUE){
  true_dark <- tall_dfs[[1]]
}
SA_graph <- matrix(SA_graph, nrow = num_rows , ncol = num_cols, byrow = TRUE) 
SA_graph <- as.data.frame(SA_graph)
true_dark <- matrix(true_dark, nrow = num_cols , ncol = num_rows, byrow = TRUE) 
true_dark <- as.data.frame(true_dark)
if (Surfacegraphing == TRUE & true_area_model == TRUE){
  true_dark <- SA_graph
}


#for troubleshooting absolute darkness
if (abs_dark_troubleshooting == TRUE | change_in_dark ==TRUE){
  abs_dark<- tall_dfs[[1]]*tall_dfs[[2]]
  
  ###used to switch to spherical model (i.e., dark=dark, eccentricity independent)
  if (true_area_model == FALSE){
    abs_dark <- tall_dfs[[2]]
  }
  
  true_dark <- matrix(abs_dark, nrow = num_cols , ncol = num_rows, byrow = TRUE) 
  true_dark <- as.data.frame(abs_dark)
} else {print ("Processing true_darkness")}

#for change in true_darkness/change in unit time (amount of shedding/unit time)
if (change_in_dark == TRUE){
  zero_list <- rep(0, num_cols)
  df_shifter <- c(zero_list,abs_dark)
  f_true_dark <- c(abs_dark, zero_list)
  dif_true_dark <- f_true_dark - df_shifter
  dif_true_dark <- tail(dif_true_dark, length(dif_true_dark)-num_cols) 
  dif_true_dark <- dif_true_dark/true_SA #adjust this so divide by average SA between timepoints?
  dif_true_dark <- head(dif_true_dark, length(dif_true_dark)-num_cols)
  
  true_dark <- matrix(dif_true_dark, nrow = num_cols , ncol = num_rows-1, byrow = TRUE) 
  true_dark <- as.data.frame(true_dark)
}

################################################################################################
#normalisation here 
################################################################################################
#tempshift back to transpose
if (change_in_dark == TRUE){
  true_dark <- as.data.frame(matrix(as.numeric(t(true_dark)), nrow=num_rows-1, ncol = (num_cols), byrow = TRUE))} else {
    true_dark <- as.data.frame(matrix(as.numeric(t(true_dark)), nrow=num_rows, ncol = num_cols, byrow = TRUE))
  }

mm <- function(x) {(x - min(x)) / (max(x) - min(x))}
minfc <- function(x) {x / min(x)}
rs<- function(x){(x- median(x)) /(quantile(x,probs = .75)-quantile(x,probs = .25))}
fc <- function(x) {x / true_dark[i,1]}
start <- true_dark[1,]


normalise <- function(true_dark) {
  for (i in 1:ncol(true_dark)){
    if (x_variable_normalisation == "mm") {
      print("min-max")
      true_dark[i] <<- as.data.frame(lapply(true_dark[i], mm))}
    else if (x_variable_normalisation == "zs") {
      true_dark <<-  as.data.frame(scale(true_dark, center = TRUE, scale = TRUE))}
    else if (x_variable_normalisation == "rs") {
      true_dark <<- as.data.frame(lapply(true_dark, rs))}
    else if (x_variable_normalisation == "minfc") {
      print("minimum-fold change")
      true_dark <<- as.data.frame(lapply(true_dark, minfc))}
    else if (x_variable_normalisation == "fc") {
      print("fold change")
      for (i2 in 1:nrow(true_dark)) {
        true_dark[i2,i] <<- (as.numeric(true_dark[i2,i]))/start[i]
        print(true_dark[i2,i])
      }}
    else if (!(x_variable_normalisation %in% templist)){
      print("No Normalisation")}
  }}
normalise(true_dark)

#tempshift back to transpose
true_dark <- as.data.frame(t(true_dark))
################################################################################################ 
#adjusting row/column labels 
tempname <- mydatarawarea[-c(1)]
tempname <- tempname[-c(1),]
if (change_in_dark == TRUE){
  tempname <- tempname[-c(1),]}
#colnames(true_dark_norm) <- tempname[,1]
true_dark <- rbind(tempname[,1], true_dark)

tempname <- mydatarawarea[-c(1,2)]
temp <- cbind("time", tempname[1,])
rownames(true_dark) <- temp

true_dark_transposed <- true_dark
##################################################################################################
#transpose for graphing - if needtrui comparison analysis use above for anova processing
#tempshift back to transpose
true_dark_transposed <- as.data.frame(t(true_dark_transposed))

###renaming to shorten column names
colnames(true_dark_transposed) <- str_replace(colnames(true_dark_transposed), "(?i)terminal ileum", "TI")
colnames(true_dark_transposed) <- str_replace(colnames(true_dark_transposed), "(?i)sigmoid column", "SC")
colnames(true_dark_transposed) <- str_replace(colnames(true_dark_transposed), "(?i)duodenum", "Duo")

# Convert the list of row means to a data frame



#for change in true_darkness/change in unit time (amount of shedding/unit time)
if (change_in_dark != TRUE){
  df_long <- pivot_longer(true_dark_transposed, cols = -time, names_to = "wells", values_to = "dark")
  df_long$time <- as.numeric(df_long$time)
  df_long$dark <- as.numeric(df_long$dark)
  darkmin <- 0.98*as.numeric(min(df_long$dark)) 
  if (x_variable_normalisation=="mm"){
    darkmax <- 1} else {
      darkmax <- 1.02*as.numeric(max(df_long$dark))}
  darkmin_r <- round(darkmin, digits = 1)
  darkmax_r <- round(darkmax, digits = 1)
  
  ##Averaging lines - -to find average line of each group 
  models <- lapply(split(df_long, df_long$wells), function(group_data) lm(dark ~ time, data = group_data))
  model_coefficients <- lapply(models, coef)
  val_intercept <- sapply(model_coefficients, "[[", "(Intercept)")
  val_slope <- sapply(model_coefficients, "[[", "time")
  avg_val_intercept <- mean(val_intercept)
  avg_val_slope <- mean(val_slope)
  
  ## actual plot
  p <- ggplot(df_long, aes(x = time, y = dark, color = wells, group=wells)) + 
    geom_line() + 
    labs(x = "Time (h)", y = y_axis_lab, title = graph_title) +
    geom_abline(intercept = avg_val_intercept, slope = avg_val_slope, color = "black", linetype="dashed") +  # Add line with specified intercept and slope
    scale_y_continuous(
      #  limits = c(darkmin_r, darkmax_r),  # Set the limits of the y-axis
      #  breaks = seq(darkmin_r, darkmax_r, by = dif_lim),  # Set the breaks (tick marks) on the y-axis
      #  labels = seq(darkmin_r, darkmax_r, by = dif_lim), # Set the labels for the breaks
      name = "Arbitrary Darkness level") +  # Set the axis label
    geom_smooth(method = "lm", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme(panel.background = element_rect(fill = "#dff5f3")) # Change background color
  print(p)
  paste("Normalisation:", x_variable_normalisation)
  paste("Darkness change/unit time:", change_in_dark)
} else if (change_in_dark == TRUE) {
  tdtcol <- ncol(true_dark_transposed)
  averageroc <- list()
  for (i in 2:tdtcol){
    averageroc[i]<- mean(as.numeric(as.list(true_dark_transposed[[i]])))
  }
  averageroc <- averageroc[-c(1)] 
  averageroc <- as.data.frame(averageroc)
  averageroc <-   rbind(colnames(true_dark_transposed[-c(1)]), averageroc)
  rownames(averageroc) <- list("sample", "dark")
  averageroc <- t(averageroc)
  averageroc <- as.data.frame(averageroc)
  averageroc$dark <- as.numeric(averageroc$dark)
  
 
  tdtcol <- tdtcol-1
  if (x_variable_normalisation == "absfc") {
    for (i in 1:tdtcol){
      print("absfc - setting lowest roc to 1")
      print(averageroc$dark[i]/min(averageroc$dark))
      averagerocmod$dark[i] <-(averageroc$dark[i]/min(averageroc$dark))
    }
    averageroc <- averagerocmod
  }
  
  
  #plot bar graph
  p <- ggplot(averageroc, aes(x = sample, y = dark, fill = sample)) +
    geom_bar(stat = "identity") +  # Create bars based on the Value column
    labs(x = "Sample", y = y_axis_lab, title = graph_title)+
    #scale_fill_manual(values = colors, name = "sample") - manual colour fill
    scale_fill_viridis_d(name = "sample")+ # Add color legend
    theme(axis.text.x = element_blank()) +# Remove x-axis labels
    theme(plot.title = element_text(face = "bold"))
  #geom_errorbar(aes(ymin = averageroc$dark - stdevaverageroc, ymax= averageroc$dark + stdevaverageroc), width = 0.4, position = position_dodge(0.9)) 
  
  p
}

#for standard deviations 
stdevaverageroc <- abs(averageroc$dark)^0.5








