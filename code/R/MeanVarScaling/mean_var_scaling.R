
options(warn = 1) # to issue warnings immediately rather than at end of loop

# Load a few utility functions
source('../lib/read_data.R')
source('../lib/study_info.R')

# Read in the table of dataset details
dataset_details <- read.csv(
  '../../../data/dataset_details.csv'
)

# Run on both types of datasets?
types <- c("One_Predator_One_Prey", "One_Predator_Two_Prey")
# types <- c("One_Predator_One_Prey")

# Initiate counts of datasets that were fit
cnts <- cntsof <- c(0, 0)

# Use no weights (i.e. equal weights) or weight treatments
# (i.e. mean-variance estimates) by each treatment's number of replicates
# (only relevant when replication varies across treatments)
weightbyrepl <- TRUE

for(tp in seq_along(types)){
  type <- types[tp]
  
  # Create empty container in which to aggregate all summarized datasets for
  # subsequent (hierarchical) analysis
  AggData <- dim(0)
  
  # Master list of datasets
  datasets <- list.files(paste0('../',type,'/Dataset_Code'), 
                         pattern=".R$", full.names=TRUE, include.dirs=FALSE)
  
  # Manually remove some atypical datasets
  # (Note: Elliot_2006 and Iyer_1996 datasets are okay despite subsequent warnings)
  rem <- c(
    grep('Krylov', datasets), # no replication for more than 1 pred density
    grep('Eveleigh_1982_ad', datasets), # high proportion of zero-variance treatments
    grep('Eveleigh_1982_ap', datasets), # high proportion of zero-variance treatments
    grep('Wong_2005', datasets) # insufficient treatment variation in both datasets
    )
  if(length(rem) > 0){  datasets <- datasets[-rem] }
  
  cntsof[tp] <- length(datasets)
  
  # datasets <- datasets[35]
  
  # Fit everything on a dataset by dataset basis
  for(i in seq_along(datasets)){
  
    # Create a short nickname for the dataset
    datasetsName <- sub('*.R$','', sub(paste0('*../',type,'/Dataset_Code/'),'', datasets[i]))
    
    # Grab info about how to find a dataset
    source(datasets[i])
    
    # Use the above information to read the data into a variable d
    d <- read.data(datadir, filename, type, columns)
    
    # Check if data has actually be read in, only then should we assess it
    if(is.null(d)){
      # Print out which dataset WILL NOT be assessed
      cat(paste0("Skipping ", "(#",i,"): ", datasetsName, "\n"))
    }else{
      # Print out which dataset WILL be assessed
      cat(paste0("Assessing ", "(#",i,"): ", datasetsName, "\n"))
      
      # Grab info about experimental design, etc
      this.study <- study.info(
        dataset_details,
        datadir,
        type
      )
      
      # Some datasets (e.g., Omkar 2014) have replicates of differing time duration.
      # for these we standardize the counts before summarizing
      # as doing so won't change the relationship between mean and variance 
      if(length(unique(d$Time)) > 1){
        if(type=='One_Predator_One_Prey'){
          d$Nconsumed <- d$Nconsumed / d$Time
        }
        if(type=='One_Predator_Two_Prey'){
          d$Nconsumed1 <- d$Nconsumed1 / d$Time
          d$Nconsumed2 <- d$Nconsumed2 / d$Time
        }
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Summarize raw datasets by treatment
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(!this.study$bootstrap){ # !bootstrap indicates raw datasets
        if(type=='One_Predator_One_Prey'){
          dsum <- aggregate(list(Nconsumed = d$Nconsumed), 
                            by = list(Npredator = d$Npredator,
                                      Nprey = d$Nprey),
                            FUN = function(x){ c(
                              mean = mean(x),
                              var = var(x),
                              se = sd(x)/length(x),
                              n = length(x)
                            )})
        }
        if(type=='One_Predator_Two_Prey'){
          dsum <- merge(
            aggregate(list(Nconsumed1 = d$Nconsumed1),
                      by = list(Npredator = d$Npredator,
                                Nprey1 = d$Nprey1,
                                Nprey2 = d$Nprey2),
                      FUN = function(x){ c(
                        mean = mean(x),
                        var = var(x),
                        se = sd(x)/length(x),
                        n = length(x)
                        )}),
            aggregate(list(Nconsumed2 = d$Nconsumed2),
                      by = list(Npredator = d$Npredator,
                                Nprey1 = d$Nprey1,
                                Nprey2 = d$Nprey2),
                      FUN = function(x){ c(
                        mean = mean(x),
                        var = var(x),
                        se = sd(x)/length(x)
                      )}))
        }
        
        dsum <- as.data.frame(as.matrix(dsum))
        colnames(dsum) <- sub('Nconsumed.n', 'n', colnames(dsum))
        colnames(dsum) <- sub('Nconsumed1.n', 'n', colnames(dsum))
        
      }else{ # Estimate variance for summarized datasets
        if(type=='One_Predator_One_Prey'){
          dsum <- d
          dsum$Nconsumed.var <- (dsum$Nconsumed.se * dsum$n)^2
        }
        if(type=='One_Predator_Two_Prey'){
          dsum <- d
          dsum$Nconsumed1.var <- (dsum$Nconsumed1.se * dsum$n)^2
          dsum$Nconsumed2.var <- (dsum$Nconsumed2.se * dsum$n)^2
        }
        
        dsum <- dsum[, !(names(dsum) %in% 'Time')] # drop 'Time' column
      }
      
      # Remove unreplicated treatments, those with zero feeding, and those with 
      # no variance in feeding (the latter two because we'll be fitting models
      # on the log-scale).  Also clears data for use based on remaining replication.
      # Only fit a dataset if there are at least 
      nMin = 4 # replicates per treatment
      #     for each of at least 
      lMin = 4 # lMin treatment levels
      #     and if estimates of uncertainty are actually provided (e.g., Katz 1982 doesn't)
      if(type=='One_Predator_One_Prey'){ # this creates dat
        source('lib/clean_one_pred_one_prey.R')
      }
      if(type=='One_Predator_Two_Prey'){ # this creates dat1 and dat2
        source('lib/clean_one_pred_two_prey.R')
      }
      
      # Fit dataset-specific models & save summarized data to aggregated data file
      if(!cleared){ # variable 'cleared' produced by cleaning scripts above
        cat(paste0("          Did not use (#",i,"):", datasetsName, "\n"))
      }else{
        
        cat(paste0(" Fitting ", "(#",i,"):", datasetsName, "\n"))
  
        if(type=='One_Predator_One_Prey'){
          source('lib/models_one_pred_one_prey.R')
          sdat <- dat
          nTrtmts <- nrow(dat)
        }
        if(type=='One_Predator_Two_Prey'){
          source('lib/models_one_pred_two_prey.R')
          sdat <- list(Prey1 = list(dat1),
                       Prey2 = list(dat2))
          nTrtmts <- list(Prey1 = nrow(dat1),
                          Prey2 = nrow(dat2))
        }
        cnts[tp] <- cnts[tp] + 1
        
        # Include summarized data in the all-summarized-data container.
        # Some rows may have been removed before model-fitting, so here we keep
        # the derived "dat" data.frame instead of the original "dsum" data.frame.
        dat <- data.frame(dataname = datasetsName,
                          replacement = this.study$replacement,
                          dat)
        AggData <- rbind(AggData, dat)
        
        # Save data-specific summarized data, study info, and model fits
         ffr.fit <- list(
           study.info = c(
             datasetName = datasetsName,
             datadir = datadir,
             data = list(sdat),
             nTrtmts = nTrtmts,
             this.study
             ),
           fits = fits,
           # profile = profile,
           ttest = ttest)
 
        saveRDS(ffr.fit, 
                file = paste0('../../../results/R/MeanVarScaling/IndividualFits/', 
                              type, '/', datasetsName, '.Rdata'))
      }
    } # end of try dataset loop
  } # end of dataset loop
  
  AggData <- data.frame(AggData)
  write.csv(AggData,
            file = paste0('../../../results/R/MeanVarScaling/AggregatedData/AggData_',
                          type, '.csv'),
            row.names = FALSE)
  
} # end of dataset type loop

cat(paste("Successfully fit", cnts, types,"of", cntsof, "datasets. \n"))

###############################################################################
###############################################################################
###############################################################################
