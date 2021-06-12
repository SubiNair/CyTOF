### Processor for CyTOF Application
#Includes gating and logging

library(flowCore)
library(tools)
library(bestNormalize)
library(gateR)
library(logr)



#Command parsing
cli <- commandArgs(trailingOnly = TRUE)
args <- strsplit(cli, "=", fixed = TRUE)


#https://stackoverflow.com/questions/28017141/select-the-last-n-columns-of-data-frame-in-r
move_to_start <- function(x, to_move) {
  x[, c(to_move, setdiff(colnames(x), to_move))]
} 

#This indicates how many conditions there are
#by default it is set to 1 but if C2 is detected then it will be changed to 2
num_c <- 1

# process all of the arguments here:

fcspath <- args[[2]][2]
correlation_val <- args[[3]][2]
alpha_val <- as.numeric(args[[4]][2])
csvpath <- args[[5]][2] #actually just the name of the csv
markerinput <- as.character(args[[6]][2])
numerator_val <- args[[7]][2]
arcsinh_transform <- args[[8]][2]
runname <- args[[9]][2] #name of the run that the user selected

#create the logfile
log_name <- paste(runname, "log", sep = '.')
lf <- log_open(file.path(fcspath, log_name))
log_print(paste("Log file for", runname, sep = ': '), console = FALSE)

log_print("Parameters", console = FALSE,  hide_notes = TRUE)
log_print(paste('Directory', fcspath, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Run name', runname, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Metadata Filename', csvpath, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Correlation', correlation_val, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Alpha', alpha_val, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Arcsinh Transformation', arcsinh_transform, sep = ': '), console = FALSE, hide_notes = TRUE)
log_print(paste('Extracting cells from case clusters', numerator_val, sep = ': '), console = FALSE, hide_notes = TRUE)

#Split up those markers that were selected based on the string we used
selected_markers <- unlist(strsplit(markerinput, "xyz"))


#Read in the files and establish a blacklist of .fcs columns we will not use
filecsv <- suppressWarnings(read.csv(paste(fcspath, csvpath, sep = "/"), header = TRUE))
metadata <- data.frame(lapply(filecsv, factor))
fs <- read.flowSet(path = fcspath, pattern = ".fcs")
blacklist <- c("Time", "Cell_length", "beadDist")

#Build out a dataframe of the markers and channels so we have a reference table
channels <- colnames(fs[[1]])
channels <- channels[!(channels %in% blacklist)]

m_dim <- length(channels)

markers <- as.data.frame(matrix(nrow = 2, ncol = m_dim))
for (ch in 1:m_dim) {
  
  if(!is.na(match(getChannelMarker(fs[[1]],channels[ch])$desc, markers[1,]))) {
    
    firstindex <- match(getChannelMarker(fs[[1]],channels[ch])$desc, markers[1,])
    suffix <- gsub("[()]", "", channels[ch])
    orig_suffix <- gsub("[()]", "", channels[firstindex])
    
    markers[1, ch] <- paste(getChannelMarker(fs[[1]],channels[ch])$desc, suffix, sep = "-")
    markers[1, firstindex] <- paste(getChannelMarker(fs[[1]],channels[ch])$desc, orig_suffix, sep = "-")
    markers[2, ch] <- channels[ch]
    
  }
  else {
    markers[1, ch] <- getChannelMarker(fs[[1]],channels[ch])$desc
    markers[2, ch] <- channels[ch]
  }
  #gsub("[()]", "", x)
}

log_print('Markers for Gating', console = FALSE, hide_notes = TRUE)

#extract the channel names needed for the analysis
vars_channels <- selected_markers
for(m in 1:length(selected_markers)) {
  vars_channels[m] <- markers[2, match(selected_markers[m], markers[1,])]
  log_print(paste('Marker-', m, ': ' ,selected_markers[m], sep = ''), console = FALSE, hide_notes = TRUE)
}

vars_channels <- make.names(vars_channels)

#full_data will be the data_frame for the gating function
#manual_colnames are the columns that we will be adding
#If a second condition is in the file, "C2" will be appended later
full_data <- NULL
manual_colnames <- c("ID", "C1")

### Building out the input data frame

for(frame in 1:length(fs)) {
  
  
  #We go through each of the files in the flowset and create a little dataframe with the necessary info
  #We then bind all of those dataframes together into full_data
  #fs_frame holds the data from each run of the loop and is cleared out each run
  
  #All of the file data
  #We then remove the columns that are not a part of the analysis -- listed in the blacklist
  fs_frame <- as.data.frame(exprs(fs[[frame]]))
  fs_frame <- fs_frame[,-which(names(fs_frame) %in% blacklist)]
  
  # if(log_ten == TRUE) {
  #   fs_frame <- log(fs_frame)
  # }
  
  #find the number of cells (rows) for the file
  fs_row <- nrow(fs_frame)
  
  #Repeat the filename for the number of rows (id_col), and bind it to the frame with the header Dataset
  #the identifier() function is from flowCore and extracts the filename from the current file
  #file_path_sans_ext is from tools and removes the .fcs extension
  id_col <- rep(identifier(fs[[frame]]), fs_row)
  fs_frame <- cbind(fs_frame, Dataset = id_col)
  file_val <- file_path_sans_ext(identifier(fs[[frame]]))
  
  #Using the metadata file we can now extract the conditions and their names
  #There is no need to rename them, they can retain their original names
  
  #C1 is needed but C2 is optional
  #We deal with C1 first
  condition1val <- metadata$C1[metadata$Filename==file_val]
  condition1col <- rep(condition1val, fs_row)
  fs_frame <- cbind(fs_frame, C1 = condition1col)
  
  #And C2 if it is present
  if('C2' %in% colnames(metadata)) {
    num_c <- 2
    c2_control <- metadata$C2[which(metadata$C2.Control==TRUE)]
    c2_case <- metadata$C2[which(metadata$C2.Control==FALSE)]
    condition2val <- metadata$C2[metadata$Filename==file_val]
    condition2col <- rep(condition2val, fs_row)
    fs_frame <- cbind(fs_frame, C2 = condition2col)
    #fs_frame$C2 <- as.factor(fs_frame$C2)
    if(!('C2' %in% manual_colnames)) {
      manual_colnames <- append(manual_colnames, "C2")
    }
  }
  
  full_data <- rbind(full_data, fs_frame)
  
}

#Add a column for the row numbers called "ID"
row_numbers <- seq.int((nrow(full_data)))
full_data <- cbind(full_data, ID=row_numbers, row.names = NULL)
full_data <- move_to_start(full_data, manual_colnames)

#should be changed to the case value instead
c1_control <- metadata$C1[which(metadata$C1.Control==TRUE)]
c1_case <- metadata$C1[which(metadata$C1.Control==FALSE)]


#arcsinh transformation
if(arcsinh_transform == TRUE) {
  
  manual_colnames <- append(manual_colnames, "Dataset")
  
  for(col in colnames(full_data)) {
    if(!col %in% manual_colnames) {
      full_data[,col] <- predict(arcsinh_x(full_data[,col]))
    }
    
  }
  
}

#single vs double log 

if(numerator_val == FALSE) {
  if('C2' %in% colnames(metadata)) {
    c2val = as.character(c2_control)
    c1val = as.character(c1_control)
  } else {c1val = as.character(c1_control)}
} else {
  if('C2' %in% colnames(metadata)) {
    c2val = as.character(c2_case)
    c1val = as.character(c1_case)
  } else {c1val = as.character(c1_case)}
}


if (length(vars_channels) %% 2 != 0 ) {
  stop(length(vars_channels) %% 2)
}

# add path specify picture saving path


#creating the path for the images to be saved - will also be used to save RDS
picture_path <- paste(fcspath, "/", sep = '')

#We can now remove the .fcs files in the folder since the object for gating has been created 
g <- dir(fcspath, pattern = '.fcs')
file.remove(file.path(fcspath, g))


#back to app for incrementing progress

if('C2' %in% colnames(metadata)) {
  final_obj <- gating(dat = full_data, n_condition = num_c, vars=vars_channels, plot_gate=TRUE, save_gate = TRUE, path_gate = picture_path, alpha = alpha_val,  p_correct = correlation_val, numerator = numerator_val, c1n=c1val, c2n=c2val)
} else {
  final_obj <- gating(dat = full_data, n_condition = num_c, vars=vars_channels, plot_gate=TRUE, save_gate = TRUE, path_gate = picture_path, alpha = alpha_val,  p_correct = correlation_val, numerator = numerator_val, c1n=c1val)
}

alt_markers <- markers
alt_markers[2,] <- make.names(markers[2,])

#add the markers from this run onto the object
final_obj[[length(final_obj) + 1]] <- alt_markers
names(final_obj)[length(final_obj)] <- "markerTable"


#prepare object path and save
rds_obj <- final_obj

log_print("The following numbers indicate the number of cells after each gate.", console = FALSE,  hide_notes = TRUE)
log_print(paste('n', rds_obj$n, sep = ': '), console = FALSE)

rds_name <- paste(runname, "RDS", sep = ".")

rds_obj[[length(rds_obj) + 1]] <- rds_name
names(rds_obj)[length(rds_obj)] <- "run"

rds_obj[[length(rds_obj) + 1]] <- selected_markers
names(rds_obj)[length(rds_obj)] <- "gatingMarkers"

rds_path <- paste(picture_path, rds_name, sep="/")
saveRDS(rds_obj, rds_path)

