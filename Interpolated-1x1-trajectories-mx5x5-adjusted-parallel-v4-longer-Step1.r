## interpolate mx 5x5 trajectories from bayesPop (after adjusting median using wpp2015 values)
## compute 1x1 life tables (including truncated at open age 100+)
## source code optimized for speed and memory using parallel computing and data.table
# HS: add profiling
## Rprof(memory.profiling = TRUE, gc.profiling = TRUE, line.profiling = TRUE)
#memory.limit(size=320000)
## require(ffbase)
## memory.size(max = TRUE)
## --max-mem-size=500M
## --max-vsize=8000M

# install.packages("bayesPop")
# install.packages("reshape2")
# install.packages("readr")
# install.packages("signal") ## for pchip interpolation at age 90+ to prevent negative numbers
# ## install.packages("optimx") ## optimiser for ICM
# install.packages("BB") ## optimiser for ICM
# install.packages("foreach") ## for parallel processing (used for ICM)
# install.packages("doParallel") ## for parallel processing (used for ICM)
# install.packages("varhandle")
# install.packages("berryFunctions")


## library(bayesPop)
library(reshape2)
library(readr)
library(signal) ## for pchip interpolation at age 90+ to prevent negative numbers
## library(optimx) ## optimiser for ICM
library(BB)
library(foreach) ## for parallel processing (used for ICM)
library(doParallel) ## for parallel processing (used for ICM)
library(varhandle)
# library(DataCombine)
library(berryFunctions) ## lsMem() ## Show memory size of the biggest objects in MB
# HS: need version 1.9.7, see https://github.com/Rdatatable/data.table/wiki/Installation
library(data.table)
library(abind)  ## for multdimensional arrays

# Install development version of data.table
## install.packages("data.table", type = "source",
##    repos = "http://Rdatatable.github.io/data.table")

# Revert back to CRAN version
## remove.packages("data.table")         # First remove the current version
## install.packages("data.table")        # Then reinstall the CRAN version


# HS: directories adapted
## popdir <- 'i:/wpp2012/pop2012-small'
## popdir <- 'g:/PPP2015/Pop/results20160130'
## programdir <- '/Users/hana/bayespop/R/PatricksInterpolation'
#programdir <- 'C:/Users/Patrick/Dropbox/bayesPop/bayesPop2015/Pop/'
#programdir <- 'C:/Users/PGIB-Home/Dropbox/bayesPop/bayesPop2015/Pop/'
#programdir <- 'F:/PPP2015'
programdir <- getwd() # script is called from its directory
## popdir <- 'E:/PPP2015/Pop/results20160130'
popdir <- file.path(programdir, 'results') # output directory
#popdir <- 'E:/PPP2015/Pop/results20160130'
#popdir <- 'C:/PPP2015/Pop/results20160130'
#popdir <- 'F:/PPP2015/Pop/results20160130'
mx5x5.prefix <- file.path(programdir, "mx5x5/UN_PPP2015_mx5x5_") # input files
pop1x1.prefix <- file.path(programdir, "pop1x1/UN_PPP2015_Pop1x1_130_by_Time_")
#mx5x5.prefix <- file.path(popdir, "mx5x5/UN_PPP2015_mx5x5_")
#pop1x1.prefix <- file.path(popdir, "pop1x1/UN_PPP2015_Pop1x1_130_by_Time_")

## outdir <- paste0(popdir, "/mx1x1")
#outdir <- "c:/mx1x1"
outdir <- file.path(popdir, "mx1x1")
#outdir <- "h:/mx1x1"

# HS: added condition and setting wrkdir
if(!file.exists(outdir)) 
	dir.create(outdir)
wrkdir <- getwd()
setwd(outdir)

nr.traj <- 1000
#nr.traj <- 10 # HS: for testing purposes set to small number 
YearMin <- 1950
YearWPP <- 2015
YearMax <- 2100
prefix <- 'UN_PPP2015_lt1x1_130'
prefix100 <- 'UN_PPP2015_lt1x1'

Sex <- c("M", "F", "B")
SexID <- c(1, 2, 3)

oage.published <- 100

source(file.path(programdir, 'BeersMSplit_function.R'))

LocID.selection <- read.delim2(file.path(programdir, 'WPP2015_TFR_nosmall_input.txt'), header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, check.names = TRUE, blank.lines.skip = TRUE, strip.white = TRUE)
LocID.selection <- subset(LocID.selection, country_code <= 84)
country.codes <- LocID.selection$country_code
country.names <- sub('/', '-', LocID.selection$country)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Mortpak ICM function based on first term of Heligman-Pollard function
## reference about the function: 1982 UN model life tables and Mortpak software
## reimplemented here using non-linear solver and lower/upper bounds on params
	## ICM <- function(t, data)
	ICM <- function(t, x, y)
	        {    				
				Age <- x
				q5 <- y
				## ICM INITIAL ESTIMATES OF PARAMETERS
	      		## t <- as.vector(c(0.2, 1.45, 0.5))
	      		## t <- as.vector(c(0.0024, 0.0512, 0.0953))
	      		## Y = ln(-ln lqx) = ln(-ln tl) + t3 * ln (x + t2)
	      		Yx <- log(-log(t[1])) + (t[3] * log((Age + t[2])))
	      		## 1qx=exp(-exp(Y))
	      		qx.fit <- exp(-exp(Yx))
	      		## input initial 1q0 (known)
	      		## data$qx.fit[data$Age==0] <- data$q5[data$Age==0]
	      		lx.fit <- c(1, cumprod(1 - qx.fit[1:10]))
				q5.fit <- c(0,0,0)
	      		q5.fit[1] <- 1 - lx.fit[Age==1] / lx.fit[Age==0]
	      		q5.fit[2] <- 1 - lx.fit[Age==5] / lx.fit[Age==1]
	      		q5.fit[3] <- 1 - lx.fit[Age==10] / lx.fit[Age==5]
	      		
	      		## GoF <- sum of the squares of the proportionate difference between the fitted and the observed nqx (like in Heligman-Pollard paper)
	      		sum(((q5 / q5.fit)-1)^2, na.rm=TRUE)
			}

			
	KeyfitzICM <- function(y, x.fitted)
	        {    				
				## x.fitted <- Age1
				## y <- input$lx 
				lx <- y
				Age <- x.fitted
				
				## For under age 5
				## B = (5 * (l1 - l5)) / (4 + l5 - (5 * l1)) 
				# B <- ((1*5) * (lx[2] - lx[3])) / ((5 - 1) + (1 * lx[3]) - (5 * lx[2]))
				B <- (5 * (lx[2] - lx[3])) / (4 + lx[3] - (5 * lx[2]))
				## A = l1 * (1+B) - B
				# A <- ((lx[2] * (1 + B)) - B) / 1
				A <- lx[2] * (1 + B) - B
				
				## for Age 5-10
				## B = ((5*10) * (l5 - l10)) / (1 + l10 - (2 * l5))
				B2 <- ((5*10) * (lx[3] - lx[4])) / ((10 - 5) + (5 * lx[4]) - (10 * lx[3]))
				## A = ((l5 * (5 + B)) - B) / 5
				A2 <- ((lx[3] * (5 + B2)) - B2) / 5
				
				## lx = ((A * x) + B) / (x + B)
				lx.fit05 <- ((A * Age) + B) / (Age + B)
				lx.fit510 <- ((A2 * Age) + B2) / (Age + B2)

	      		## observed
				q5 <- (1 - (shift(lx, n=1, fill=NA, type="lead")) / lx)[1:3]

				## fitted	
				## under-5
				qx.fit1 <- (1 - (shift(lx.fit05, n=1, fill=NA, type="lead")) / lx.fit05)
				## age 5-10
				qx.fit2 <- (1 - (shift(lx.fit510, n=1, fill=NA, type="lead")) / lx.fit510)

				qx.fit <- qx.fit1
				## age 2-3
				qx.fit[3] <- 0.75 * qx.fit1[3] + 0.25 * qx.fit2[3]
				## age 3-4
				qx.fit[4] <- 0.5 * qx.fit1[4] + 0.5 * qx.fit2[4]
				## age 4-5
				qx.fit[5] <- 0.25 * qx.fit1[5] + 0.75 * qx.fit2[5]
				## age 5-6
				qx.fit[6] <- 0.1 * qx.fit1[6] + 0.9 * qx.fit2[6]
				## age 6+
				qx.fit[7:10] <- qx.fit2[7:10]
				
	      		lx.fit <- c(1, cumprod(1 - qx.fit))
	      		## lx.fit1 <- c(1, cumprod(1 - qx.fit1[1:10]))
	      		## lx.fit2 <- c(1, cumprod(1 - qx.fit2[1:10]))
                
				q5.fit <- rep(0,3)
	      		q5.fit[1] <- 1 - lx.fit[Age==1] / lx.fit[Age==0]
	      		q5.fit[2] <- 1 - lx.fit[Age==5] / lx.fit[Age==1]
	      		q5.fit[3] <- 1 - lx.fit[Age==10] / lx.fit[Age==5]
                ## 
				## q5.fit1 <- rep(0,3)
	      		## q5.fit1[1] <- 1 - lx.fit1[Age==1] / lx.fit1[Age==0]
	      		## q5.fit1[2] <- 1 - lx.fit1[Age==5] / lx.fit1[Age==1]
	      		## q5.fit1[3] <- 1 - lx.fit1[Age==10] / lx.fit1[Age==5]
                ## 
				## q5.fit2 <- rep(0,3)
	      		## q5.fit2[1] <- 1 - lx.fit2[Age==1] / lx.fit2[Age==0]
	      		## q5.fit2[2] <- 1 - lx.fit2[Age==5] / lx.fit2[Age==1]
	      		## q5.fit2[3] <- 1 - lx.fit2[Age==10] / lx.fit2[Age==5]
				## 
	      		## ## GoF <- sum of the squares of the proportionate difference between the fitted and the observed nqx (like in Heligman-Pollard paper)
	      		GoF <- sum(((q5 / q5.fit)-1)^2, na.rm=TRUE)
	      		## GoF1 <- sum(((q5 / q5.fit1)-1)^2, na.rm=TRUE)
	      		## GoF2 <- sum(((q5 / q5.fit2)-1)^2, na.rm=TRUE)
				## print(GoF)
				## print(GoF1)
				## print(GoF2)
                ## 
				## 				plot(Age1, log(qx.fit1), type="b")
				## 				lines(Age1, log(qx.fit2), type="b", col="red")
				## 				lines(Age1[1:10], log(qx.fit), type="b", col="green")
				## 				lines(lt1$Age, log(lt1$qx), type="b", col="blue")
				
				return(list(A1=A, B1=B, A2=A2, B2=B2, y.fitted = lx.fit, GoF = GoF))  
			}
			

i <- 0
# HS: use one country as a test 
## country.codes2 <- 4
## for(LocID in country.codes2) {
for(LocID in country.codes) {
	i <- i + 1
	
	# HS: check if mx5x5 file exists
	mx5.file.name <- paste0(mx5x5.prefix, LocID, '.csv')
	if(!file.exists(mx5.file.name)) next
	
	## LocID <- 4
	print(paste0("LocID: ", LocID, " (", round(100*i/length(country.codes),0), "% of countries)"))

	## mx5 <- read_csv('C:/PPP2015/Pop/results20160130/mx/WPP.wide.csv', col_names = TRUE, skip = 0, n_max = -1, progress = interactive())
	## mx5 <- subset(mx5, SexID < 3)
	## mx5 <- melt(mx5, id=c("LocID", "SexID", "Age"), variable="MidPeriod", value.name="mx")
	## mx5$AgeGrpStart <- mx5$Age
	## mx5$Age <- NULL
	## mx5$Trajectory <- 1
	## Sex <- unique(mx5$SexID)

	strt<-Sys.time()
	
	## mx5 <- data.table(read_csv(mx5.file.name, col_names = TRUE, skip = 0, n_max = -1, progress = interactive()))
	mx5 <- fread(paste0(mx5x5.prefix, LocID, '.csv'), header = TRUE, skip = 0, nrows=-1, showProgress=TRUE, data.table=TRUE)
	## mx5 <- subset(mx5, Trajectory <= nr.traj & SexID <= 2)
	mx5 <- mx5[Trajectory <= nr.traj & SexID <= 2]

	Sex <- unique(mx5$SexID)
	Age5 <- unique(mx5$AgeGrpStart)
	Trajectory <- unique(mx5$Trajectory)
	MidPeriod <- unique(mx5$MidPeriod)

## Step 1: ICM - Mortpak graduation for under age 10 (needs to convert into qx)
	mx10.long <- mx5[AgeGrpStart < 15]
	setnames(mx10.long, "AgeGrpStart", "Age")
	mx10.long[, c("AgeGrp", "AgeGrpSpan", "Period") := NULL]
	setkeyv(mx10.long, c("LocID", "SexID", "Age", "MidPeriod", "Trajectory"))
	
	## Step 2a. convert abridged mx into qx

	## ax by sex based on CD-West
	## warning: mx0 must be known and used in all equations below for conditions which is not the case right now!!!
	mx0 <- mx10.long[Age==0]
	setnames(mx0, "mx", "mx0")
	mx0[, c("Age") := NULL]

	mx10.long <- merge(mx10.long, mx0, by.x = c('LocID', 'SexID', 'MidPeriod', 'Trajectory'), by.y = c('LocID', 'SexID', 'MidPeriod', 'Trajectory'), all.x=TRUE, all.y=FALSE)
	mx0 <- NULL

	
	## Male
	## 1a0
	mx10.long[SexID==1 & Age==0 & mx0 < 0.107, ax := 0.045 + 2.684 * mx]
	mx10.long[SexID==1 & Age==0 & mx0 >= 0.107, ax := 0.33]

	## 4a1
	mx10.long[SexID==1 & Age==1 & mx0 < 0.107, ax := 1.651 - 2.816 * mx]
	mx10.long[SexID==1 & Age==1 & mx0 >= 0.107, ax := 1.352]


	## Female
	## 1a0
	mx10.long[SexID==2 & Age==0 & mx0 < 0.107, ax := 0.053 + 2.8 * mx]
	mx10.long[SexID==2 & Age==0 & mx0 >= 0.107, ax := 0.35]
	mx10.long[Age==0, n := 1]

	## 4a1
	mx10.long[SexID==2 & Age==1 & mx0 < 0.107, ax := 1.522 - 1.518 * mx]
	mx10.long[SexID==2 & Age==1 & mx0 >= 0.107, ax := 1.361]
	mx10.long[Age==1, n := 4]


	## 5a5
	mx10.long[Age>=5, c("ax", "n") := list(2.5, 5)]
	
	# nqx = n * mx / (1 + (n - ax) * mx)
	mx10.long[, qx := n * mx / (1 + (n  - ax) * mx)]

	Age1 <- as.vector(c(0:10))
	
##  Alternative longer form: fit ICM Mortpak model but slow - use analytical form (Keyfitz rationale function) only in case of poor fit

	#setup parallel backend to use no_cores processors
	## memory.limit()/memory.size() = max cores
	# Calculate the number of cores
	no_cores <- min(detectCores(), nr.traj) ## - 1
	## no_cores <- 3
 
	# Initiate cluster
	cl <- makeCluster(no_cores, type = "PSOCK")
	registerDoParallel(cl)

	#start time
	## strt<-Sys.time()

	## make the data/variables available to each cluster
	clusterExport(cl, c("mx10.long", "Trajectory", "MidPeriod"))

	## loop thru the set of trajectories by splitting automatically the set by # of clusters
	## recombine the results using rbind
	## note: must share data/ariables used in the clusters before calling function, and define the list of packages used (otherwise only base R available)
	## run ICM for each trajectory: for each sex, and then by period using default stating values for first period, and last period params for folllowing periods


	
	## tryCatch({
	u10lt <- NULL
	## u10lt <- foreach(traj = Trajectory, .packages=c("optimx", "BB"), .combine=rbind) %dopar% {
	u10lt <- foreach(traj = Trajectory, .packages=c("data.table", "BB"), .combine=rbind) %dopar% {
		## traj <- 2

		u10lt <- NULL
		## for (traj in Trajectory){
		for (j in 1:2) { ## by sex for male and female
			## j <- 1

				for (l in MidPeriod) {
						 ## l <- 2013
						 ## print(paste0("Step1: ICM - LocID: ", LocID, " SexID: ", j, " Trajectory: ", traj, " MidPeriod: ", l, " (", round(100*i/length(country.codes),0), "% of countries)"))
						 input <- mx10.long[SexID==j & Trajectory==traj & MidPeriod==l,]
						 input$lx <- c(1, cumprod(1 - input$qx[1:3]))

							lt1 <- data.frame(matrix(vector(), 11, 5,
  						              dimnames=list(c(), c("Age", "lx", "q5", "Yx", "q5.fit"))),
  						              stringsAsFactors=F)
							lt1$Age <- Age1
							lt1$lx[lt1$Age==0] <- input$lx[input$Age==0]
							lt1$lx[lt1$Age==1] <- input$lx[input$Age==1]
							lt1$lx[lt1$Age==5] <- input$lx[input$Age==5]
							lt1$lx[lt1$Age==10] <- input$lx[input$Age==10]
  						
							lt1$q5[lt1$Age==0] <- input$qx[input$Age==0]
							lt1$q5[lt1$Age==1] <- input$qx[input$Age==1]
							lt1$q5[lt1$Age==5] <- input$qx[input$Age==5]

							q5 <- c(input$qx[input$Age==0], input$qx[input$Age==1], input$qx[input$Age==5])
						 
						 	## Option #1 using Mortpak ICM function -- slow using to finding parameters
							## call ICM function (use 4 digit precision on params to speed up nlm, and use initial params from previous time period, except for first period)
							## results <- optimx(t, ICM, data = lt1, lower=c(0.0000001, 0.0000001, 0.013), upper=c(0.184, 3, 0.5), method="nlm", control=list(step.min=1e-4, rel.tol=1e-4))
							## results <- optimx(t, ICM, x=Age1, y=q5, lower=c(0.0000001, 0.0000001, 0.013), upper=c(0.184, 3, 0.5), method="nlm", control=list(step.min=1e-4, rel.tol=1e-4))
							## results <- optim(par=t, fn=ICM, x=Age1, y=q5, lower=c(0.0000001, 0.0000001, 0.013), upper=c(0.184, 3, 0.5), method="L-BFGS-B", 
							##			control=list(maxit=10000, factr=1e-4))
							## based on best fitted trajectories
							## results <- optimx(t, ICM, data = lt1, lower=c(0.0000005, 0.0000001, 0.02), upper=c(0.025, 0.15, 0.25), method="nlm", control=list(rel.tol=1e-4))
													
							## ICM INITIAL ESTIMATES OF PARAMETERS (at the beginning of the time series)
							## t <- as.vector(c(0.2, 1.45, 0.3))
							t <- as.vector(c(0.0024, 0.0512, 0.0953))
							results <- BBoptim(par=t, fn=ICM, x=Age1, y=q5, lower=c(0.0000001, 0.0000001, 0.013), upper=c(0.184, 3, 0.5), method = c(3, 2, 1), 
										control = list(gtol = 1e-02, trace=FALSE), quiet=TRUE)
							
							## use for next time period as initial params the values from the last period (to prevent potential jumps in time series only due to initial params)
							t <- as.vector(c(results$par[1], results$par[2], results$par[3]))
							## compute back ICM results using params
							lt1$Yx <- log(-log(t[1])) + (t[3] * log((lt1$Age + t[2])))
							lt1$qx <- exp(-exp(lt1$Yx))
                        
							lt1$lx <- c(1, cumprod(1 - lt1$qx[1:10]))
							## lt1$q5.fit[lt1$Age==0] <- 1 - lt1$lx[lt1$Age==1] / lt1$lx[lt1$Age==0]
							## lt1$q5.fit[lt1$Age==1] <- 1 - lt1$lx[lt1$Age==5] / lt1$lx[lt1$Age==1]
							## lt1$q5.fit[lt1$Age==5] <- 1 - lt1$lx[lt1$Age==10] / lt1$lx[lt1$Age==5]
							
							## check for basic constraints on qx
							## all qx > 0 (i.e., (lt1$qx<0))
							## and all 1qx(1+) > 1q0 (i.e., lt1$qx[lt1$Age>0] > lt1$qx[lt1$Age==0] )
							## if at least 1 true, then error or if results$value is a bad fit
							results2 <- list(A1=NA, B1=NA, A2=NA, B2=NA, y.fitted = NA, GoF = results$value)
							if ((sum(lt1$qx<0) > 0) | (sum(lt1$qx[lt1$Age>0] > lt1$qx[lt1$Age==0]) > 0) | (results$value > 0.001)) {
							
							# Option #2: speed up using only deterministic function.
							# use instead Keyfitz' rational function for child survival
								results2 <- KeyfitzICM(y=input$lx, x.fitted=Age1)
								lt1$lx2 <- results2$y.fitted[1:11]
								lt1$qx2 <- (1 - (shift(lt1$lx, n=1, fill=NA, type="lead")) / lt1$lx)
								
								## lt1$q5.fit[lt1$Age==0] <- 1 - lt1$lx[lt1$Age==1] / lt1$lx[lt1$Age==0]
								## lt1$q5.fit[lt1$Age==1] <- 1 - lt1$lx[lt1$Age==5] / lt1$lx[lt1$Age==1]
								## lt1$q5.fit[lt1$Age==5] <- 1 - lt1$lx[lt1$Age==10] / lt1$lx[lt1$Age==5]
								
							}
							if ((sum(lt1$qx<0) > 0) | (sum(lt1$qx[lt1$Age>0] > lt1$qx[lt1$Age==0]) > 0) | (results$value > results2$GoF)) {
								lt1$lx <- lt1$lx2
								lt1$qx <- lt1$qx2
							}
										
							## substitute input 1q0
							lt1$qx[lt1$Age==0] <- lt1$q5[lt1$Age==0]      				
							## recompute lx
							lt1$lx <- c(1, cumprod(1 - lt1$qx[1:10]))    				
							## compute ax
							lt1$ax[lt1$Age==0] <- input$ax[input$Age==0]
							lt1$ax[lt1$Age>0] <- 0.5
							## compute mx from qx							
							lt1$mx <- -lt1$qx / ((1 * lt1$qx) - (lt1$ax * lt1$qx) - 1)
							
							## lt1 <- cbind(lt1, par1=t[1], par2=t[2], par3=t[3], GoF=results$value)
							## u10lt <- rbind(u10lt, cbind(LocID, SexID=j, MidPeriod=l, Trajectory=traj, subset(lt1, select=c(Age, mx, qx, q5, q5.fit, par1, par2, par3, GoF))))
							u10lt <- rbind(u10lt, cbind(LocID, SexID=j, MidPeriod=l, Trajectory=traj, subset(lt1, Age < 10, select=c(Age, mx))))
					} # next year
			} # next sex
 		return(u10lt)	
	} # next trajectory
	## }, error = function(e) return(paste0("The trajectory '", traj, "'", " SexID ", j, " MidPeriod ", l, " caused the error: '", e, "'")))
	## print(Sys.time()-strt)
	## head(fits)
	
	stopCluster(cl)

	u10lt <- data.table(u10lt)
	rm("mx10.long")
	
	


	##  ## Tried Keyfitz rationale function as default - but SexID=2 and MidPeriod=2098 & Trajectory = 427 and 623 gives NA
    ##  
	##  ## reshape with Trajectory in column
	##  mx10.wide <- dcast(mx10.long, LocID + SexID + MidPeriod + Trajectory + mx0 ~ Age, value.var=c("mx", "ax", "qx"))	
	##  rm("mx10.long")
	##  setorder(mx10.wide, LocID, SexID, Trajectory, MidPeriod)
    ##  
	##  mx10.wide[, lx_0 := 1]
	##  mx10.wide[, lx_1 := lx_0 * (1-qx_0)]
	##  mx10.wide[, lx_5 := lx_1 * (1-qx_1)]
	##  mx10.wide[, lx_10 := lx_5 * (1-qx_5)]
    ##  
	##  ## For under age 5
	##  ## B = (5 * (l1 - l5)) / (4 + l5 - (5 * l1)) 
	##  mx10.wide[, B1 := (5 * (lx_1 - lx_5)) / (4 + lx_5 - (5 * lx_1))]
	##  ## A = l1 * (1+B) - B
	##  mx10.wide[, A1 := lx_1 * (1+B1) - B1]
	##  
	##  ## for Age 5-10
	##  ## B = ((5*10) * (l5 - l10)) / (1 + l10 - (2 * l5))
	##  mx10.wide[, B2 := ((5*10) * (lx_5 - lx_10)) / ((10 - 5) + (5 * lx_10) - (10 * lx_5))]
	##  ## A = ((l5 * (5 + B)) - B) / 5
	##  mx10.wide[, A2 := ((lx_5 * (5 + B2)) - B2) / 5]
	##  
    ##  
	##  ## populate fitted life table
	##  ID <- mx10.wide[,1:4, with=FALSE]
	##  u10lt <- NULL
	##  for (i in Age1) {
	##  		u10lt <- rbind(u10lt, cbind(ID, Age=i))
	##  }
	##  ID <- NULL
    ##  
	##  u10lt <- merge(u10lt, mx10.wide, by.x = c('LocID', 'SexID', 'MidPeriod', 'Trajectory'), by.y = c('LocID', 'SexID', 'MidPeriod', 'Trajectory'), all.x=TRUE, all.y=FALSE)
	##  u10lt[, c("mx_0", "mx_1", "mx_5", "mx_10", "ax_1", "ax_5", "ax_10", "qx_1", "qx_5", "qx_10", "lx_0", "lx_1") := NULL]
	##  setorder(u10lt, LocID, SexID, Trajectory, MidPeriod, Age)
	##  
	##  
	##  ## lx = ((A * x) + B) / (x + B)
	##  u10lt[, lx.fit1 := ((A1 * Age) + B1) / (Age + B1)]
	##  u10lt[, lx.fit2 := ((A2 * Age) + B2) / (Age + B2)]
    ##  
	##  ## fitted	
	##  ## under-5
	##  u10lt[, qx.fit1 := (1 - (shift(lx.fit1, n=1, fill=NA, type="lead")) / lx.fit1), by=.(LocID, SexID, MidPeriod, Trajectory)]
	##  ## age 5-10
	##  u10lt[, qx.fit2 := (1 - (shift(lx.fit2, n=1, fill=NA, type="lead")) / lx.fit2), by=.(LocID, SexID, MidPeriod, Trajectory)]
    ##  
	##  ## combined
	##  u10lt[, qx := qx.fit1]
	##  ## Age 0 substitute with original value
	##  u10lt[Age == 0, qx := qx_0]
	##  ## age 2-3
	##  u10lt[Age==2, qx := 0.75 * qx.fit1 + 0.25 * qx.fit2]
	##  ## age 3-4
	##  u10lt[Age==3, qx := 0.5 * qx.fit1 + 0.5 * qx.fit2]
	##  ## age 4-5
	##  u10lt[Age==4, qx := 0.25 * qx.fit1 + 0.75 * qx.fit2]
	##  ## age 5-6
	##  u10lt[Age==5, qx := 0.1 * qx.fit1 + 0.9 * qx.fit2]
	##  ## age 6+
	##  u10lt[Age > 5, qx := qx.fit2]
	##  
	##  ## deal with exception if lx==1, i.e., qx=0 then B1-A2 = NaN and qx = NaN
	##  u10lt[Age < 5 & lx_5==1, qx := 0]
	##  u10lt[lx_10==1, qx := 0]
    ##  
    ##  
	##  ## recompute lx
	##  u10lt[, lx := c(1, cumprod(1 - qx)), by=.(LocID, SexID, MidPeriod, Trajectory)]
	##  ## compute ax
	##  u10lt[, ax := 0.5]
	##  u10lt[Age == 0, ax := ax_0]
	##  ## compute mx from qx							
	##  u10lt[, mx := -qx / ((1 * qx) - (ax * qx) - 1), by=.(LocID, SexID, MidPeriod, Trajectory)]
	##  ## Age 0 substitute with original value
	##  u10lt[Age == 0, mx := mx0]
	##  ## keep only under age 10
	##  u10lt <- u10lt[Age < 10, ]
	##  ## keep only mx
	##  u10lt[, c("mx0", "ax_0", "qx_0", "B1", "A1", "B2", "A2", "lx_5", "lx_10", "lx.fit1", "lx.fit2", "qx.fit1", "qx.fit2", "qx", "lx", "ax") := NULL]
	##  setorder(u10lt, LocID, SexID, Trajectory, MidPeriod, Age)

	print(paste("step: 1", Sys.time()-strt))



## Step 2: interpolate abridged mx by sex over age using pchip (fitted under age 10 using ICM fitted mx)
	strt<-Sys.time()
## 2a. combine with ICM for mx age 10 and under

	## use original mx for last open age group
	oage <- max(mx5$AgeGrpStart)
	mx5.long.oage <- mx5[AgeGrpStart==oage,]
	setnames(mx5.long.oage, "AgeGrpStart", "Age")
	mx5.long.oage[, c("AgeGrp", "AgeGrpSpan", "Period") := NULL]
	
	## use mx for 5-year age groups from 10 onward
	mx5icm <- mx5[AgeGrpStart >= 10,]
	mx5icm[, MidAge := AgeGrpStart + 2.5]
	mx5icm[, c("AgeGrp", "AgeGrpStart", "AgeGrpSpan", "Period") := NULL]
	rm("mx5")
	
	## use ICM interpolated values up to age 10
	setnames(u10lt, "Age", "MidAge")	
	u10lt[, MidAge := MidAge + 0.5]
	
	## append ICM data for age 10 and under
	mx5icm <- rbindlist(list(mx5icm, u10lt), use.names=TRUE)
	
	rm("u10lt")
	
 	## order rows
	setorder(mx5icm, LocID, SexID, MidPeriod, Trajectory,  MidAge)
	

	MidAge <- unique(mx5icm$MidAge)
	Age1 <- c(min(MidAge):max(MidAge))

##	## try to run interp1 directly on data.table without loop -- but not clear how to filter subset for other data.table within function
##	## populate fitted life table
##	ID <- mx5icm[MidAge==0.5, 1:4, with=FALSE]
##	mx5x1.long <- NULL
##	for (i in Age1) {
##			mx5x1.long <- rbind(mx5x1.long, cbind(ID, Age=i))
##	}
##	ID <- NULL
##	setorder(mx5x1.long, LocID, SexID, Trajectory, MidPeriod, Age)
##	
##	test <- mx5x1.long[Trajectory ==1 & SexID==1 & MidPeriod==2013,]
##	setorder(test, LocID, SexID, Trajectory, MidPeriod, Age)
##	input <- mx5icm[Trajectory ==1 & SexID==1 & MidPeriod==2013,]
##	setorder(input, LocID, SexID, Trajectory, MidPeriod, MidAge)
##  
##	test[, yinterpol := mapply(
##    function(x, y, x2) interp1(x, log(y), x2, method = 'pchip', extrap = TRUE),
##    input[, MidAge],
##    input[, mx],
##    test[, Age])
##    ,	by=.(LocID, SexID, MidPeriod, Trajectory)]

	#start time
	## strt<-Sys.time()

	# Initiate cluster
	no_cores <- min(detectCores(), nr.traj) ## - 1
	cl <- makeCluster(no_cores, type = "PSOCK")
	registerDoParallel(cl)
	## make the data/variables available to each cluster
	clusterExport(cl, c("LocID", "mx5icm", "Trajectory", "MidPeriod", "Age1"))

	mx5x1.long <- NULL
	mx5x1.long <- foreach(k = Trajectory, .packages=c("data.table", "signal"), .combine=rbind) %dopar% {
			mx5x1.long <- NULL
			for (j in 1:2) { ## by sex for male and female
				for (l in MidPeriod) {
					 ## input <- mx5icm[SexID==2 & MidPeriod==2098 & Trajectory==427,]
					 input <- mx5icm[SexID==j & Trajectory==k & MidPeriod==l,]
					 ## deal with exception if exactly 0 and log
					 input$mx[input$mx<=0] <- (1/100000000)
					 yinterpol <- interp1(input$MidAge, log(input$mx), Age1, 'pchip', extrap = TRUE)
					 
					 mx5x1.long <- rbind(mx5x1.long, cbind(LocID, SexID=j, MidPeriod=l, Age=Age1-0.5, Trajectory=k, mx=exp(yinterpol)))
				}
		}
		return(mx5x1.long)
	} # next trajectory
	
	stopCluster(cl)
	## print(Sys.time()-strt)

	
	mx5x1.long <- data.table(mx5x1.long)
	rm("mx5icm")
	## subset(mx5x1.long, SexID==1 & Trajectory==1 & MidPeriod==2013)

	
	## use original mx for last open age group
	mx5x1.long.oage <- rbindlist(list(mx5x1.long[Age < oage,], mx5.long.oage), use.names=TRUE)
	rm("mx5x1.long", "mx5.long.oage")
	gc()

	## sort records
	setorder(mx5x1.long.oage, LocID, SexID, MidPeriod, Trajectory, Age)

	## reshape with Trajectory in column
	mx5x1.wide <- dcast(mx5x1.long.oage, LocID + SexID + Age + MidPeriod ~ Trajectory, value.var=c("mx"))	
	setorder(mx5x1.wide, LocID, SexID, Age, MidPeriod)
	rm("mx5x1.long.oage")
	# mx5x1.wide$MidPeriod <- as.numeric(mx5x1.wide$MidPeriod)
	Age1 <- unique(mx5x1.wide$Age)

	# test
	# mx5x1.wide2 <- mx5x1.wide
    # mx5x1.wide <- mx5x1.wide2[, 1:54, with=FALSE]
	print(paste("step: 2", Sys.time()-strt))
	
## Step 3: interpolate mx5x1 by sex over time using Beers Modified six-term formula
	
	## process matrix by sex (male and female)

	strt<-Sys.time()
	no_cores <- min(detectCores(), nr.traj) ## - 1
	cl <- makeCluster(no_cores, type = "PSOCK")
	registerDoParallel(cl)
	clusterExport(cl, c("mx5x1.wide", "Age1"))
	mx1x1 <- NULL
	mx1x1 <- foreach(k = Age1, .packages=c("data.table"), .combine=rbind) %dopar% {
				mx1x1 <- NULL
				for (j in 1:2) { ## by sex for male and female
					input <- as.data.frame(mx5x1.wide[SexID==j & Age==k,])
				
					## interpolate data using Beers Modified
					data <- input[,5:ncol(input)]
					rownames(data) <- input$MidPeriod
					output.wide <- BeersMSplit(data)   ## interpolated over time and get back mid-year results
					output.wide <- as.data.frame(output.wide) * 5  ## to rescale back (per default function is for counts, not rates and split input by 1/5)
							
					## combine back IDs and column headers
					Year <- c((min(input$MidPeriod-3)):(max(input$MidPeriod+1)))
					output.wide <- cbind(LocID=LocID, SexID=j, Year, Age=k, output.wide)
					mx1x1 <- rbind(mx1x1, output.wide)
				}
				return(mx1x1)
			}
	stopCluster(cl)

	
	mx1x1 <- data.table(mx1x1)
	rm("mx5x1.wide")

	mx1x1.long <- melt(mx1x1, id=c("LocID", "SexID", "Year", "Age"), variable="Trajectory", value.name="mx")
	rm("mx1x1")
	# HS: use data.table for ordering
	#mx1x1.long <- mx1x1.long[order(mx1x1.long$LocID, mx1x1.long$SexID, mx1x1.long$Year, mx1x1.long$Trajectory, mx1x1.long$Age),]
	#stop('')
	setorder(mx1x1.long, LocID, SexID, Year, Trajectory, Age)
	
	## censor year before start of pop. projections (i.e., pop open age group = 100)
	mx1x1.long <- mx1x1.long[Year >= YearWPP, ]

	## reshape with Year in column
	mx1x1.wide <- dcast(mx1x1.long, LocID + SexID + Trajectory + Age ~ Year, value.var=c("mx"))
	rm("mx1x1.long")
	# HS: use data.table for ordering
	#mx1x1.wide <- mx1x1.wide[order(mx1x1.wide$LocID, mx1x1.wide$SexID, mx1x1.wide$Trajectory, mx1x1.wide$Age),]
	setorder(mx1x1.wide, LocID, SexID, Trajectory, Age)
	
	## mx1x1.long <- melt(mx1x1.wide, id=c("LocID", "SexID", "Age", "Trajectory"), variable="Year", value.name="mx")


	## clean memory ls() and rm() and gc() 
	## list=ls()
	## rm("mx10.long", "mx1x1", "mx1x1.long", "mx5", "mx5.long.oage", "mx5icm", "mx5x1.long", "mx5x1.long.oage", "mx5x1.wide", "u10lt", "u10lt.wide")
	gc()
	print(paste("step: 3", Sys.time()-strt))

	## mx1x1.wide <- cbind(mx1x1.wide[,1:4], mx1x1.wide[,10:ncol(mx1x1.wide)])
	
	
	
## Step 4: compute mx 1x1 for both sexes

	strt<-Sys.time()

## step 1: compute population exposure
## read adjusted pop 1x1
	## strt<-Sys.time()
	pop1x1 <- data.table(read_csv(paste0(pop1x1.prefix, LocID, '.csv'), col_names = TRUE, skip = 0, n_max = -1, progress = interactive()))
	## read_csv faster than fread for this bigger file
	## pop1x1 <- fread(paste0(pop1x1.prefix, LocID, '.csv'), header = TRUE, skip = 0, nrows=-1, showProgress=TRUE, data.table=TRUE)
	## print(Sys.time()-strt)

	pop1x1 <- pop1x1[SexID<3 & Trajectory <= nr.traj, ]
	setorder(pop1x1, LocID, SexID, Trajectory, Age)
	## pop1x1

	Exposure1x1 <- cbind(pop1x1[,(1:4), with=FALSE], (pop1x1[,5:(ncol(pop1x1)-1), with=FALSE] + pop1x1[,6:ncol(pop1x1), with=FALSE]) / 2)
	## head(Exposure1x1)
	rm("pop1x1")

## Step 2: compute mx for both sexes
	## if 0 then impute a minimum lower value to avoid divide by zero
	nc <- ncol(Exposure1x1)
	data <- Exposure1x1[, 5:nc, with=FALSE]
	data[data==0] <- 0.001
	Exposure1x1 <- cbind(Exposure1x1[, 1:4, with=FALSE], data)
	rm("data")
	gc()

	## strt<-Sys.time()
	## print(Sys.time()-strt)
	## lsMem()
	
	## initialize BS recordset using male recordset
	mx1x1.wide.BS <- mx1x1.wide[SexID==1,]
	## impute SexID
	mx1x1.wide.BS[, SexID := 3]
	## compute mx for both sexes using population exposure by sex
	## with mx and exposure already pre-sorted
	# HS: extracting columns from data.table needs with=FALSE
	#mx1x1.wide.BS[,5:nc] <- ((mx1x1.wide[mx1x1.wide$SexID==1,5:nc] * Exposure1x1[Exposure1x1$SexID==1,5:nc]) + (mx1x1.wide[mx1x1.wide$SexID==2,5:nc] * Exposure1x1[Exposure1x1$SexID==2,5:nc])) / (Exposure1x1[Exposure1x1$SexID==1,5:nc] + Exposure1x1[Exposure1x1$SexID==2,5:nc])
	
	mx1x1.wide.BS <- cbind(mx1x1.wide.BS[, 1:4, with=FALSE], ((as.matrix(mx1x1.wide[SexID==1, 5:nc, with=FALSE]) * Exposure1x1[SexID==1, 5:nc, with=FALSE]) + (as.matrix(mx1x1.wide[SexID==2, 5:nc, with=FALSE]) * Exposure1x1[SexID==2, 5:nc, with=FALSE])) / (Exposure1x1[SexID==1, 5:nc, with=FALSE] + Exposure1x1[SexID==2, 5:nc, with=FALSE]))
	
	mx1x1.wide <- rbind(mx1x1.wide, mx1x1.wide.BS)
	rm("mx1x1.wide.BS", "Exposure1x1")
	gc()
	print(paste("step: 4", Sys.time()-strt))

	
## Step	5: compute complete life table lt1x1
## note: experimented with different inner loops sequences and parallel processing to speed up and deal with memory limitation issues
## right now by Sex > Years and LT computed for all trajectories in columns at once and for each LT indicator use different data frames
## parallel processing for each of these sex x years allows to speed up and clean-up memory more efficiently
## Afterward reshaping of dataset into long format to get back indicators in columns and each trajectory in row
## export results by sex to keep file size manageable
## create 2 versions: untruncated and full precision, and public version rounded to significant digits and truncated at open age group (oage.published=100)
	
	strt<-Sys.time()

	## reshape life table into long format
	lt1.long <- melt(mx1x1.wide, id=c("LocID", "SexID", "Age", "Trajectory"), variable="Year", value.name="mx")
	rm("mx1x1.wide")
	# HS
	#lt1.long <- lt1.long[order(lt1.long$LocID, lt1.long$SexID, lt1.long$Year, lt1.long$Trajectory, lt1.long$Age),]
	setorder(lt1.long, LocID, SexID, Year, Trajectory, Age)

	## reshape with Trajectory in column
	mx1x1.wide <- dcast(lt1.long, LocID + SexID + Year + Age ~ Trajectory, value.var=c("mx"))
	rm("lt1.long")
	# HS
	#mx1x1.wide <- mx1x1.wide[order(mx1x1.wide$LocID, mx1x1.wide$SexID, mx1x1.wide$Year, mx1x1.wide$Age),]
	setorder(mx1x1.wide, LocID, SexID, Year, Age)
	gc()

	Year <- unique(mx1x1.wide$Year)
	print(paste("step: 5a", Sys.time()-strt))
	
	## mx1x1.wide <- as.data.frame(mx1x1.wide)
	for (j in SexID) {
		## loop in parallel thru the years
		strt<-Sys.time()
		no_cores <- min(detectCores(), nr.traj) ## - 1
		cl <- makeCluster(no_cores, type = "PSOCK")
		registerDoParallel(cl)
		clusterExport(cl, c("mx1x1.wide", "j", "Year"))
		lt1 <- NULL
		lt1 <- foreach(k = Year, .packages=c("data.table"), .combine=rbind) %dopar% {
			lt1 <- NULL
			lt <- mx1x1.wide[SexID==j & Year==k,]
			## for memory/speed efficiency, works with data frame for each LT indicator to be able to apply more efficient R on the whole frames and avoid loops
			lt.mx <- as.data.frame(lt[, 5:ncol(lt), with=FALSE])
	        
			## mx for open age /* reciprocal of 1/ex for open age group*/ 
	        
			## max age
			nmax <- max(lt$Age)
			
			## compute ax
			lt.ax <- lt.mx
			lt.ax[,] <- 0.5 # (default like HMD)
			## for under age 1: CD-West
			## Male
			if (j==1){
				lt.ax[1,lt.mx[1,] < 0.107] <- 0.045 + 2.684 * lt.mx[1, lt.mx[1,] < 0.107]
				lt.ax[1,lt.mx[1,] >= 0.107] <- 0.33
			}
			## Female & Both sexes
			if (j>=2){
				lt.ax[1,lt.mx[1,] < 0.107] <- 0.053 + 2.8 * lt.mx[1, lt.mx[1,] < 0.107]
				lt.ax[1,lt.mx[1,] >= 0.107] <- 0.35
			}
			## for open age group ax = ex
			
			# n width of the age intervals
			lt.n <- lt.mx
			lt.n[,] <- c(diff(lt$Age),-1)
			lt.n[lt$Age == nmax,] <- -1
			
			# nqx = n * mx / (1 + (n - ax) * mx)
			lt.qx <- lt.mx
			lt.qx <- (lt.n * lt.mx) / (1 + (lt.n  - lt.ax) * lt.mx)
			lt.qx[lt.qx>1] <- 1 # necessary for high nMx
			## open age /* cannot be computed due to trunaction (would be 1 is open age group is higher)*
			lt.qx[lt$Age == nmax,] <- NA 
			
			## Px <- (1 - Qx)
			
			# lx
			lt.lx <- lt.mx
			lt.lx[,] <- rbind(1, cumprod(1 - lt.qx[1:nrow(lt.lx),]))
			lt.lx <- lt.lx * 100000
			
			# dx
			lt.dx <- lt.lx * lt.qx
			## open age group dx = lx
			lt.dx[lt$Age == nmax,] <- lt.lx[lt$Age == nmax,]
			
			# Lx
			lt.Lx <- (lt.n * lt.lx) - (lt.ax * lt.dx)       # equivalent to n*l(x+n) + (n-nax)*ndx
			## open age group Lx = Tx
			lt.Lx[lt$Age==nmax,] <- lt.lx[lt$Age==nmax,] * lt.ax[lt$Age==nmax,]
			
			# Sx
			# for open age group, e.g., 85+ 
			## penultimate age group -> Last entry of S(x,n) is S( 80+,5) = T( 85) / T( 80)
			## for open age group iteself: Sx cannot be computed due to trunaction Sx <- NA
			
			# Tx
			lt.Tx <- lt.Lx
			rownames(lt.Tx) <- lt$Age
			lt.Tx <- lt.Tx[rev(rownames(lt.Tx)),]					
			lt.Tx <- cumsum(lt.Tx)
			lt.Tx <- lt.Tx[rev(rownames(lt.Tx)),]
			rownames(lt.Tx) <- lt$Age
			
			# ex
			lt.ex <- lt.Tx / lt.lx
			## for open age group mx=1/ex and ax = ex
			## lt.mx[lt$Age==nmax,] <- 1 / lt.ex[lt$Age==nmax,]
			## lt.ex[lt$Age==nmax,] <- 1 / lt.mx[lt$Age==nmax,]
			lt.ax[lt$Age==nmax,] <- lt.ex[lt$Age==nmax,]
				
			lt.all <- rbind(
			cbind(lt[,3:4, with=FALSE], Indicator="mx", lt.mx),
			cbind(lt[,3:4, with=FALSE], Indicator="qx", lt.qx),
			cbind(lt[,3:4, with=FALSE], Indicator="lx", lt.lx),
			cbind(lt[,3:4, with=FALSE], Indicator="dx", lt.dx),
			cbind(lt[,3:4, with=FALSE], Indicator="Lx", lt.Lx),
			cbind(lt[,3:4, with=FALSE], Indicator="Tx", lt.Tx),
			cbind(lt[,3:4, with=FALSE], Indicator="ex", lt.ex),
			cbind(lt[,3:4, with=FALSE], Indicator="ax", lt.ax))
			
			lt1 <- rbind(lt1, lt.all)
			return(lt1)
		} # next year
		stopCluster(cl)	
			
		## rm("mx1x1.wide")
		gc()
		

		## reshape life table into long format (note: this gets pretty big - need > 16GB virtual RAM)\
		# HS: use data.table for melting and ordering
		lt1 <- data.table(lt1)
		lt1.long <- melt(lt1, id=c("Year", "Age", "Indicator"), variable="Trajectory", value.name="Value")
		#lt1.long <- lt1.long[order(lt1.long$Year, lt1.long$Trajectory, lt1.long$Indicator, lt1.long$Age),]
		setorder(lt1.long, Year, Trajectory, Indicator, Age)
		
		## lsMem()
		## object.size(lt1.long)
		rm("lt1")
		gc()
		
		## convert factors into numeric variables in data.table
		## lt1.long$Year <- as.numeric(levels(lt1.long$Year)[lt1.long$Year])
		## lt1.long$Trajectory <- as.numeric(levels(lt1.long$Trajectory)[lt1.long$Trajectory])
		## changeCols <- colnames(lt1.long)[which(as.vector(lt1.long[,lapply(.SD, class)]) == "numeric")]
		## changeCols <- c("Age", "Year", "Trajectory")
		## for (col in changeCols) lt1.long[, (col) := as.integer(lt1.long[[col]])]
		
		
		
		## ## reshape with Indicators in column (problem with memory limit)
		## ## loop in parallel thru the trajectories
		## strt<-Sys.time()
		## no_cores <- detectCores() - 1
		## cl <- makeCluster(no_cores, type = "PSOCK")   ## , outfile=""
		## registerDoParallel(cl)
		## clusterExport(cl, c("lt1.long", "Trajectory"))
		## lt1x1.wide <- foreach(k = Trajectory, .packages=c("reshape2"), .combine=rbind) %dopar% {
		## 			lt <- subset(lt1.long, Trajectory==k)					
		## 			lt <- dcast(lt, Year + Trajectory + Age ~ Indicator, value.var=c("Value"))
		## 			return(lt)
		## 			} # next trajectory
		## stopCluster(cl)	
		## print(Sys.time()-strt)
		
		lt1x1.wide <- dcast(lt1.long, Year + Trajectory + Age ~ Indicator, value.var=c("Value"))
	
		rm("lt1.long")
		gc()
		
		# HS: data.table ordering of rows
		#lt1x1.wide <- lt1x1.wide[order(lt1x1.wide$Year, lt1x1.wide$Trajectory, lt1x1.wide$Age),]
		setorder(lt1x1.wide, Year, Trajectory, Age)
		
		## record columns and add LocID back
		
		## attach(lt1x1.wide, warn.conflicts = FALSE)
		# HS: added with=FALSE (data.table syntax)
		lt1x1.wide <- lt1x1.wide[ , c("Year", "Trajectory", "Age", "mx", "qx", "lx", "dx", "Lx", "Tx", "ex", "ax"), with=FALSE]
		lt1x1.wide <- cbind(LocID=LocID, SexID=j, lt1x1.wide)
	
		## round to significant number of digits
		lt1x1.wide[, mx := round(mx, 5)]
		lt1x1.wide[, qx := round(qx, 5)]
		lt1x1.wide[, lx := round(lx, 0)]
		lt1x1.wide[, dx := round(dx, 0)]
		lt1x1.wide[, Lx := round(Lx, 0)]
		lt1x1.wide[, Tx := round(Tx, 0)]
		lt1x1.wide[, ex := round(ex, 3)]
		lt1x1.wide[, ax := round(ax, 3)]
	
		Sex <- c("M", "F", "B")

		## save untruncated life table
		#write_csv(lt1x1.wide, paste(prefix, "_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col_names = TRUE)
		#HS
		fwrite(lt1x1.wide, paste0(prefix, "_", Sex[j], "_", LocID, ".csv"), append = FALSE, col.names = TRUE, turbo = TRUE)
		
		## truncate life table with open age group at 100+ (oage.published)
		lt1x1.wide <- lt1x1.wide[Age <= oage.published, ]
	
		## mx for open age /* reciprocal of 1/ex for open age group*/ 
		lt1x1.wide[Age==oage.published, mx := (1 / lt1x1.wide[Age==oage.published, ex])]
		
		## qx
		lt1x1.wide[Age==oage.published, qx := NA] 
		
		# lx
		
		## open age group dx = lx
		lt1x1.wide[Age==oage.published, dx := lt1x1.wide[Age==oage.published, lx]]
		
		## open age group Lx = Tx
		lt1x1.wide[Age==oage.published, Lx := lt1x1.wide[Age==oage.published, Tx]]
		
		# Sx
		# for open age group, e.g., 85+ 
		## penultimate age group -> Last entry of S(x,n) is S( 80+,5) = T( 85) / T( 80)
		## for open age group iteself: Sx cannot be computed due to trunaction Sx <- NA
		
		# Tx
		
		# ex
	
		## for open age group ax = ex
		lt1x1.wide[Age==oage.published, ax := lt1x1.wide[Age==oage.published, ex]]
	
		## round to significant number of digits (based on radix for lx = 100000 - default used in Mortpak & HMD life tables)	
		lt1x1.wide[, mx := round(mx, 5)]
		lt1x1.wide[, qx := round(qx, 5)]
		lt1x1.wide[, lx := round(lx, 0)]
		lt1x1.wide[, dx := round(dx, 0)]
		lt1x1.wide[, Lx := round(Lx, 0)]
		lt1x1.wide[, Tx := round(Tx, 0)]
		lt1x1.wide[, ex := round(ex, 3)]
		lt1x1.wide[, ax := round(ax, 3)]

	
		## last step export mx
		#write_csv(lt1x1.wide, paste(prefix100, "_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col_names = TRUE)
		#HS
		fwrite(lt1x1.wide, paste(prefix100, "_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		## changeCols <- c("Year")
		## for (col in changeCols) lt1x1.wide[, (col) := as.factor(lt1x1.wide[[col]])]
		
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("mx"))
		setorder(data, LocID, SexID, Trajectory, Age)
		
		# HS: all writes below can be written like this (instead of using the if statement):
		# fwrite(data, paste0("UN_PPP2015_mx1x1_by_Time_", LocID,".csv"), append = j>1, col.names = j==1)
		fwrite(data[SexID==j], paste("UN_PPP2015_mx1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)

	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("qx"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_qx1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("lx"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_lx1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("Tx"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_Tx1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("Lx"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_Yx1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("ex"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_ex1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		data <- dcast(lt1x1.wide, LocID + SexID + Trajectory + Age ~ Year, value.var=c("ax"))
		setorder(data, LocID, SexID, Trajectory, Age)
		fwrite(data[SexID==j], paste("UN_PPP2015_ax1x1_by_Time_", Sex[j], "_", LocID, ".csv",sep=""), append = FALSE, col.names = TRUE, turbo = TRUE)
	
		rm("lt1x1.wide", "data")
		gc()
		print(paste("step: 5 - Sex:", Sex[j], "   ", Sys.time()-strt))
	
	} # next sex
	rm("mx1x1.wide")
	gc()
} # next country

# HS: 
setwd(wrkdir)
## Rprof(NULL)
## print(summaryRprof(lines="show", memory="both"))