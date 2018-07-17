#####################################################################################
# Check to ensure all usable data sets are in fact in use (i.e. have scripts written).
#####################################################################################
require(RCurl)
options(stringsAsFactors=FALSE)

# specify where the data files are located
dropboxdir <- '../../../dropbox_data/Data' # Stouffer
dropboxdir <- '~/Dropbox/Research/Projects/GenFuncResp/Data' # Novak

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find the folders that should have usable datasets in them
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
googledoc<-read.csv(
  text=getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vQcykYqM8Pkmgrlmp9S2jorZZEOlZ14a0AINRuDc2Y_29f6dTR9ojhOOBV2rcattJO5LXA5ATVn_nK6/pub?gid=0&single=true&output=csv"),
  header=T,
  skip=1,
  sep=",",
  na.strings = c("*", "NA")
)

folders2use <- googledoc[googledoc$Variable_Predator_Abundance == TRUE,]$Dataset_Folder
folders2use <- folders2use[!is.na(folders2use) & folders2use!='pdf'  & folders2use!='Unused']

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find the filenames of the datasets in these folders
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datasets <- Sys.glob(paste0(dropboxdir,'/*/*.csv'))  # Use instead of list.files in order to only go one subfolder in.
datasetfolders <- sapply(strsplit(dirname(datasets),'/'),tail,1)
datasets2use <- basename(datasets[datasetfolders %in% folders2use])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What datasets do we already have scripts for
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scriptfiles <- list.files('./Dataset_Code',full.names=FALSE)
scriptfiles.inprep <- scriptfiles[grepl('zzz', scriptfiles)]
  print(scriptfiles.inprep)
scriptfiles <- scriptfiles[!grepl('zzz', scriptfiles)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Which scripts must still be written?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inhand <- tolower(sub('.csv','',datasets2use))
written <- tolower(sub('.R','',scriptfiles))

done <- datasets2use[inhand %in% written]
needed <- datasets2use[!inhand %in% written] # May include filenames with typos!
needfix <- scriptfiles[!written %in% inhand]

length(done)
  print(done)
length(needfix)
  print(needfix)
length(needed)
  print(needed)

################################################################################
################################################################################
################################################################################
