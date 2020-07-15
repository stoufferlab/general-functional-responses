# reads in a dataset
# first checks if there is a file within the repository
# second reads the file from dropbox when no such file exists
read.data <- function(datadir, filename, datatype, dropboxdir=""){
	fname <- paste("../../../data", datatype, datadir, filename,sep="/")
  	if(file.exists(fname)){
  		rawdata <- read.csv(fname)
  	}else{
  		rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))
  	}
  	return(rawdata)
}