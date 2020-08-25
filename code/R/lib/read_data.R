# reads in a dataset
# first checks if there is a file within the repository
# second reads the file from dropbox when no such file exists
# if no option works then returns NULL
read.data <- function(datadir, filename, datatype, columns=NULL){
	fname <- paste("../../../data", datatype, datadir, filename,sep="/")
	if(file.exists(fname)){
		rawdata <- read.csv(fname,row.names=1)
	}else{
		# specify where the data files are located
		dropboxdir <- switch(
			Sys.getenv("LOGNAME"),
			stouffer = '~/Dropbox/Projects/GenFuncResp/Data',
			marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data',
			'gobbledygook'
		)
		fname <- paste(dropboxdir,datadir,filename,sep="/")
		if(file.exists(fname)){
			rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))
			rawdata <- rawdata[,columns[,2]]
			colnames(rawdata) <- columns[,1]
		}else{
			rawdata <- NULL
		}
	}
	return(rawdata)
}