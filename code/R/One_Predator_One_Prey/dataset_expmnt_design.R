# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Various aspects of experimental design
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# - Identify the (number of) unique experimental designs.
# - Calculate the maximum number (percent) of prey that are consumed
#     at the highest prey density treatment (distinguishing between
#     replacement and non-replacement studies).

rm(list = ls())
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(Sys.getenv("LOGNAME"),
                     stouffer = '~/Dropbox/Projects/GenFuncResp/Data',
                     marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
source('../lib/study_info.R')
source('../lib/read_data.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in the table of dataset details
dataset_details <- read.csv('../../../data/dataset_details.csv')

# master list of datasets
datasets <- list.files(
  './Dataset_Code',
  pattern = ".R$",
  full.names = TRUE,
  include.dirs = FALSE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Places to store the designs of all analysed datasets
all.designs <- list(0)
design.dims <- dim(0)

# Place to store max number of prey eaten
prey.eaten <- list(0)

j <- 1
for (i in seq_len(length(datasets))) {
  source(datasets[i])
  
  # read the data
  d <-
    read.data(datadir, filename, "One_Predator_One_Prey", columns)
  
  # grab info about experimental design, etc
  this.study <- study.info(dataset_details,
                           datadir,
                           "One_Predator_One_Prey")
  
  # dataset
  datasetsName <-
    sub('*.R$', '', sub('*./Dataset_Code/', '', datasets[i]))
  datasetsName <- gsub('_', ' ', datasetsName)
  
  # study treatment design
  design <-
    data.frame(Npredator = d$Npredator, Nprey = d$Nprey)#,Time=d$Time)
  design <- unique(design)
  design <- design[do.call(order, design),]
  
  replacement <- this.study$replacement
  
  all.designs[[j]] <- list(
    dataname = datasetsName,
    design = design,
    replacement = this.study$replacement
  )
  design.dims <- rbind(design.dims, dim(design))
  
  # maximum prey treatment and number of prey eaten at it
  c.Nprey <- grep('Nprey', colnames(d))
  c.Ncons <- grep('Nconsumed', colnames(d))
  max.prey.trtmts <- which(d[, c.Nprey] == max(d[, c.Nprey]))
  max.Nprey <- d[max.prey.trtmts[1], c.Nprey]
  max.Ncons <- max(d[max.prey.trtmts, c.Ncons])
  max.perCons <-
    ifelse(is.integer(max.Nprey), max.Ncons / max.Nprey, NA)
  
  prey.eaten[[j]] <- list(
    dataname = datasetsName,
    replacement = this.study$replacement,
    max.Nprey = max.Nprey,
    max.Ncons = max.Ncons,
    max.perCons = max.perCons
  )
  j <- j + 1
  
}

#~~~~~~~~~~~~~~~~~~~~~~
# Percent of prey eaten
#~~~~~~~~~~~~~~~~~~~~~~
prey.eaten = as.data.frame(t(sapply(prey.eaten,
                                    function(x)
                                      x[1:max(lengths(prey.eaten))])))

par(mfrow = c(1, 2),
    yaxs = 'i',
    xaxs = 'i')
hist(
  as.numeric(subset(prey.eaten, replacement == TRUE)$max.perCons),
  breaks = 20,
  main = 'Repl = TRUE',
  col = 'grey',
  xlab = 'Fraction of max(N) prey eaten in max(N) treatment'
)
hist(
  as.numeric(subset(prey.eaten, replacement == FALSE)$max.perCons),
  breaks = 20,
  main = 'Repl = FALSE',
  col = 'grey',
  xlab = 'Fraction of max(N) prey eaten in max(N) treatment'
)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Experimental treatment designs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimum possible number of designs to deal with based on number of treatments alone
nrow(design.dims)
nrow(unique(design.dims))

# Identify unique designs
matequal <- function(x, y) {
  dim(x) == dim(y) && all(x == y)
}

unique.designs <- list(0)
unique.designs[[1]] <- all.designs[[1]]
rem.designs <- 2:length(all.designs)

for (i in rem.designs) {
  di <- all.designs[[i]]$design
  for (j in 1:length(unique.designs)) {
    match <- FALSE
    dj <- unique.designs[[j]]$design
    if (matequal(di, dj)) {
      # yes
      match <- TRUE
      unique.designs[[j]]$dataname <-
        c(unique.designs[[j]]$dataname,
          all.designs[[i]]$dataname)
      rem.designs <- rem.designs[!rem.designs %in% i]
    }
  }
  if (!match) {
    # no and no
    unique.designs[[length(unique.designs) + 1]] <-
      all.designs[[i]]
    rem.designs <- rem.designs[!rem.designs %in% i]
  }
}

length(unique.designs)

# There are 68 unique designs (of 81) when considering NPred,NPrey,Time
# There are 60 unique designs (of 81) when considering NPred,NPrey

# Postscript:  Could have flattened all matrices into a single array and applied duplicated().

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot experimental designs
# ~~~~~~~~~~~~~~~~~~~~~~~~~
library(sfsmisc) # for eaxis()

alpha <- function(col, alpha) {
  col <- as.vector(col2rgb(col) / 255)
  rgb(col[1], col[2], col[3], alpha)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
focal.designs <- unique.designs

# # Restrict to replacement studies
# repl <- sapply(focal.designs, function(x) {
#   x$replacement == TRUE
# })
# focal.designs <- focal.designs[repl]
# 
# # Restrict to integer-valued studies
# repl <- sapply(focal.designs, function(x) {
#   any(is.integer(c(x$design$Nprey, x$design$Npredator))) == TRUE
# })
# focal.designs <- focal.designs[repl]

maxAbunds <- data.frame(t(sapply(focal.designs, function(x) {
  c(
    Prey = max(x$design$Nprey),
    Pred = max(x$design$Npredator),
    Ratio = max(x$design$Nprey / x$design$Npredator)
  )
})))

cairo_pdf('../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_ExperimentalDesigns_Abunds.pdf', height=3.5,width=5)
par(mar=c(3,3,0.5,0.5), 
    mgp=c(1.5,0.3,0), 
    tcl=-0.2, 
    las=1, 
    cex=0.8,
    cex.lab=1.2)
plot(1,1,
  xlim = c(1, max(maxAbunds$Prey)),
  ylim = c(1, max(maxAbunds$Pred)),
  xlab = 'Prey abundances',
  ylab = 'Predator abundances',
  log = 'xy',
  type = 'n',
  axes = FALSE
)
eaxis(1)
axis(2, las = 2)
box(lwd = 1)

lapply(focal.designs, function(x) {
  hpts <- chull(x$design$Nprey,
                x$design$Npredator)
  polygon(
    x$design$Nprey[hpts],
    x$design$Npredator[hpts],
    border = NA,
    col = alpha('black', 0.1)
  )
})
lapply(focal.designs, function(x) {
  points(
    x$design$Nprey,
    x$design$Npredator,
    pch = 21,
    col = 'black',
    bg = 'white',
    cex = 0.4
  )
})
dev.off()


# Numbers of prey and predator levels per design
AbundLevels <- data.frame(t(sapply(focal.designs, function(x) {
  c(PreyLevels = length(unique(x$design$Nprey)),
    PredLevels = length(unique(x$design$Npredator)))
})))


cairo_pdf('../../../results/R/OnePredOnePrey_figs/OnePredOnePrey_ExperimentalDesigns_AbundLevelCount.pdf', height=3.5,width=5)
    x = AbundLevels$PreyLevels
    y = AbundLevels$PredLevels
    xlab = "Prey abundance levels"
    ylab = "Predator abundance levels"
    zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
    xhist = hist(x, plot=FALSE, breaks = seq(0, max(x), 1))
    yhist = hist(y, plot=FALSE, breaks = seq(0, max(y), 1))
    top = max(c(xhist$counts, yhist$counts))
    par(mar = c(3,3,0,0),
        bty = 'o',
        tcl = -0.3,
        mgp = c(1.75, 0.3, 0),
        cex.lab = 1)
    plot(x,y,
         xlim = c(0, max(x)),
         ylim = c(0, max(y)),
         xlab = xlab,
         ylab = ylab,
         pch=21,
         bg='grey')
    axis(3, labels = FALSE, tcl = 0.3)
    axis(4, labels = FALSE, tcl = 0.3)
    box(lwd=1)
    par(mar=c(0,3,1,1),
        bty='l')
    barplot(xhist$counts, 
            ylim=c(0, top), 
            space=0,
            ylab = 'Freq.')
    box(lwd=1)
    par(mar=c(3,0,1,1),
        bty='l')
    barplot(yhist$counts, 
            xlim=c(0, top), 
            space=0, 
            horiz=TRUE,
            xlab = 'Freq.')
    box(lwd=1)
    par(oma=c(3,3,0,0))
dev.off()