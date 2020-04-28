# Identify the (number of) unique experimental designs.
# Also calculate the maximum number (percent) of prey that are consumed
# at the highest prey density treatment (distinguishing between
# replacement and non-replacement studies).

rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
  Sys.getenv("LOGNAME"),
  stouffer = '../../../dropbox_data/Data',
  marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
source('../lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Places to store the designs of all analysed datasets
all.designs <- list(0)
design.dims <- dim(0)

# Place to store max number of prey eaten
prey.eaten <- list(0)

j<-1
for(i in seq_len(length(datasets))){
  source(datasets[i])
  
  # grab info from the google doc
  this.study <- study.info(datadir)
  
  # dataset
  datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))
  datasetsName <- gsub('_',' ', datasetsName)
  
  # used in our analyses?
  if(grepl("H|R", this.study$runswith)){

    # study treatment design
    design <- data.frame(Npredator=d$Npredator,Nprey=d$Nprey)#,Time=d$Time)
    design <- unique(design)
    design <- design[do.call(order, design), ]
    
    all.designs[[j]] <- list(dataname = datasetsName, 
                             design = design)
    design.dims <- rbind(design.dims, dim(design))
    
    # maximum prey treatment and number of prey eaten at it
    c.Nprey <- grep('Nprey',colnames(d))
    c.Ncons <- grep('Nconsumed',colnames(d))
    max.prey.trtmts <- which(d[,c.Nprey]==max(d[,c.Nprey]))
    max.Nprey <- d[max.prey.trtmts[1],c.Nprey]
    max.Ncons <- max(d[max.prey.trtmts, c.Ncons])
    max.perCons <- ifelse(is.integer(max.Nprey), max.Ncons/max.Nprey, NA)
    
    prey.eaten[[j]] <- list(dataname = datasetsName,
                            replacement = this.study$replacement,
                            max.Nprey = max.Nprey,
                            max.Ncons = max.Ncons,
                            max.perCons = max.perCons)
    j <- j+1
  }
}

#~~~~~~~~~~~~~~~~~~~~~~
# Percent of prey eaten
#~~~~~~~~~~~~~~~~~~~~~~
prey.eaten = as.data.frame(t(sapply(prey.eaten, function(x) x[1:max(lengths(prey.eaten))])))

par(mfrow=c(1,2), yaxs='i', xaxs='i')
  hist(as.numeric(subset(prey.eaten, replacement==TRUE)$max.perCons), 
       breaks=20, main='Repl = TRUE', col='grey',
       xlab='Fraction of max(N) prey eaten in max(N) treatment')
  hist(as.numeric(subset(prey.eaten, replacement==FALSE)$max.perCons), 
       breaks=20, main='Repl = FALSE', col='grey',
       xlab='Fraction of max(N) prey eaten in max(N) treatment')
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Experimental treatment designs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimum possible number of designs to deal with based on number of treatments alone
nrow(design.dims)
nrow(unique(design.dims))

# Identify unique designs
matequal <- function(x, y){
  dim(x) == dim(y) && all(x == y)
}

unique.designs <- list(0)
unique.designs[[1]] <- all.designs[[1]]
rem.designs <- 2:length(all.designs)

for(i in rem.designs){
  di <- all.designs[[i]]$design
  for(j in 1:length(unique.designs)){
    match<-FALSE
    dj <- unique.designs[[j]]$design
    if(matequal(di,dj)){
      # yes
      match<-TRUE
      unique.designs[[j]]$dataname <- c(unique.designs[[j]]$dataname,
                                            all.designs[[i]]$dataname)
      rem.designs <- rem.designs[!rem.designs%in%i]
    }}
  if(!match){
      # no and no
      unique.designs[[length(unique.designs)+1]] <- all.designs[[i]]
      rem.designs <- rem.designs[!rem.designs%in%i]
    }
}

length(unique.designs)

# There are 68 unique designs (of 81) when considering NPred,NPrey,Time
# There are 60 unique designs (of 81) when considering NPred,NPrey

# Postscript:  Could have flattened all matrices into a single array and applied duplicated().



