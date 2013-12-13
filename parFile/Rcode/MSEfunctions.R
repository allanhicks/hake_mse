simulateSurveyIndex.fn <- function(dat,yr,batage,selex,qSurv,M,SElog,SEtotal,verbose=FALSE) {
#SElog is the individual year standard error (in log space)
#SEtotal is the total standard error that generated the index (SElog + additional)
    survB <- sum(batage*selex)
    survI <- rlnorm(1,log(survB*qSurv*exp(-0.5*M)),SEtotal)   #midyr (but I didn't take out any catch)
    if(verbose) cat("yr                      = ",yr,"\n",
                    "survB                   = ",survB,"\n",
                    "qSurv                   = ",qSurv,"\n",
                    "M                       = ",M,"\n",
                    "survB*qSurv*exp(-0.5*M) = ",survB*qSurv*exp(-0.5*M),"\n",sep="")
    dat$N_cpue <- dat$N_cpue+1
    dat$CPUE <- rbind(dat$CPUE,c(yr,1,2,survI,SElog))
    return(dat)
}

simulateAgeComp.fn <- function(dat,yr,natage,selex,M,nAge,nAgeAdj,fleet,ageingError) {
    #simulate age comp data from a multinomial
    #nAge is the assumed sample size for the data
    #nAgeAdj is the tuning factor that SS3 uses in the Variance adjustment factors
    #Therefore, the effective sample size is nAgeAdj*nAge
    #Ageing error is a matrix of dimension 21 X 15, 
    # where 15 is the plus group for data, and 20 is plus group for population (includes age zero)
    #The numbers at age to sample from are nAge*nAgeAdj
    #The probabilities at age are natage*selex%*%AgeingError
    natage <- as.numeric(natage)  #I think that it is read in as a factor
    selex <- as.numeric(selex)  #I think that it is read in as a factor
    probs <- matrix(natage*selex,nrow=1)%*%ageingError    #Ageing error should be done after the multinomial to mimic the true process!!!!! (Future fix)
    AgeComp <- rmultinom(1,nAge*nAgeAdj,probs)[,1]   #midyr doesn't matter since scaled to sum to 1, unless take out some catch)
    #AgeComp[16] <- sum(AgeComp[16:21])
    #AgeComp <- AgeComp[-c(1,17:21)] #remove age zeros and ages beyond plus group
    #AgeComp <- AgeComp/sum(AgeComp)  #normalize, although it doesn't need to be
    dat$N_agecomp <- dat$N_agecomp + 1
    dat$agecomp <- rbind(dat$agecomp,c(yr,1,fleet,0,0,40,-1,-1,nAge,AgeComp))
    return(dat)
}

assess.fn <- function(i,yr,cmd,convCrit,verbose=T,maxIter=3,jit=1e-5,jitMult=10) {
    #expects you to be in the proper working directory
    tmp <- shell(cmd,intern=T)
    #Do convergence checks
    x <- readLines("ss3.par",n=1)
    x <- strsplit(x,"=")
    finalGrad <- as.numeric(x[[1]][length(x[[1]])])
    iter<-0
    while(finalGrad>convCrit & iter<maxIter) {
        if(verbose) cat("redoing assessment due to non-convergence. Sim",i,"Yr",yr,"\n")
        iter <- iter+1
        starter <- SS_readstarter(verbose=F)
        #put in some jitter
        starter$jitter_fraction <- jit
        SS_writestarter(starter,overwrite=T,verbose=F)
        tmp <- shell(cmd,intern=T)
        x <- readLines("ss3.par",n=1)
        x <- strsplit(x,"=")
        finalGrad <- as.numeric(x[[1]][length(x[[1]])])
        jit <- jit*jitMult
    }
    flush.console()
    return(finalGrad)
}

removeAssessJunk.fn <- function() {
    file.rename("ss3.par","tmp.par")
    file.remove(dir()[grep("ss3.",dir())])
    file.rename("tmp.par","ss3.par")
    file.remove(dir()[grep("posterior",dir())])
    file.remove(dir()[grep("ss_new",dir())])
    file.remove(dir()[grep("admodel",dir())])
    file.remove("SS3.exe","covar.sso","echoinput.sso","SIS_table.sso","ParmTrace.sso","runnumber.ss","checkup.sso","CumReport.sso","eigv.rpt","fmin.log","rebuild.sso","variance")
}

controlFile.fn <- function() {
    #adds a year to the main recruit dev section, and bias correction upper ramp
    ctl <- readLines("2012_hake_control.ss")
    ind <- grep("End year standard recruitment devs",ctl)
    x <- as.numeric(strsplit(ctl[ind],"#")[[1]][1])
    ctl[ind] <- paste(x+1,"# End year standard recruitment devs")
    ind <- grep("Last year for full bias correction in_MPD",ctl)
    x <- as.numeric(strsplit(ctl[ind],"#")[[1]][1])
    ctl[ind] <- paste(x+1,"# Last year for full bias correction in_MPD",sep="\t")
    ind <- grep("First_recent_yr_nobias_adj_in_MPD",ctl)
    x <- as.numeric(strsplit(ctl[ind],"#")[[1]][1])
    ctl[ind] <- paste(x+1,"# First_recent_yr_nobias_adj_in_MPD",sep="\t")
    
    writeLines(ctl,"2012_hake_control.ss")
}

parFile.fn <- function(yr) {
    #modify par file from previous year to work for current year
    #add first fore dev into main recr dev
    #shift fore devs and add a zero at end
    #Make sure to call this function in the folder that you are going to do an assessment
    #it will look a level up to find previous year's assessment
    pars <- readLines(paste("../assess",yr-1,"/ss3.par",sep=""))
    ind1 <- grep("recdev1:",pars)
    ind2 <- grep("Fcast_recruitments:",pars)
    mainDevs <- as.numeric(strsplit(pars[ind1+1]," ")[[1]])[-1]
    foreDevs <- as.numeric(strsplit(pars[ind2+1]," ")[[1]])[-1]
    mainDevs <- c(mainDevs,foreDevs[1])
    foreDevs <- c(foreDevs[-1],0)
    pars[ind1+1] <- paste(mainDevs,collapse=" ")
    pars[ind2+1] <- paste(foreDevs,collapse=" ")
    
    writeLines(pars,"ss3.par")

    #modify starter file to read par file
    starter <- SS_readstarter(verbose=F)
    starter$init_values_src <- 1
    SS_writestarter(starter,overwrite=T,verbose=F)
}

doHakeMSE.fn <- function(theDir,maxYr=2030,
                         surveyPer,surveySElog=0.085,surveySEtotal=0.42,surveyNage=65,surveyNageAdj=1.0,
                         commNage=800,commNageAdj=0.12,minCatch=100,
                         wtAtAge=c(0.03,0.0885,0.2562,0.3799,0.4913,0.5434,0.5906,0.662,0.7215,0.791,0.8629,0.9315,0.9681,1.0751,1.0016,1.0202,1.0202,1.0202,1.0202,1.0202,1.0202),
                         nSims=999,convCrit=0.1,assessCmd="ss3 -nohess",theSeed,ageingError="AgeingError.csv")
{
    #the final gradient matrix needs to be saved. Nothing is currently done with it
    origWD <- getwd()
    set.seed(theSeed)
    setwd(theDir)

    selex <- read.table("posterior_mse.sso",header=T)
    selex <- split(selex,selex$Fleet)

    if(is.null(ageingError)) {
        ageingError <- matrix(0,ncol=15,nrow=21,dimnames=list(paste("a",0:20,sep=""),paste("a",1:15,sep="")))
        diag(ageingError[2:15,]) <- 1
        ageingError[1,1] <- 1   #age zeros assigned to age 1 (make sure that selectivity is working appropriately
        ageingError[16:21,15] <- 1  #plus group ages assigned to data plus group of 15
        print("Using no ageing error")
    } 
    if(file.exists(ageingError)) {
        ageingError <- t(as.matrix(read.csv("AgeingError.csv",row.names=1)))
    } else {
        stop("No ageingError matrix defined. Use NULL if you want no ageing error")
    }
    
    assessYrs <- 2012:maxYr
    
    if(length(nSims)>1) { #a vector means the user supplied the specific sim numbers (rows of MCMC)
        sims <- nSims
        cat("If you are continuing, you should copy the old catch.csv to combine with these results.\nCatch.csv will be overwritten\n")
    }
    if(length(nSims)==1) { #a single number means randomly sample that many sims
        sims <- sort(sample(2:1000,nSims))
    }
    cat("Sim Numbers:",sims,"\n")
    finalGradient <- matrix(NA,nrow=length(sims),ncol=length(assessYrs),dimnames=list(as.character(sims),as.character(assessYrs)))    #if nSims>1, probably want to read in last

    params <- read.table("posteriors.sso",header=T)  #posterior parameter should not change with each mceval, so read it in now
    catch <- matrix(NA,nrow=length(sims),ncol=maxYr-2012+2,dimnames=list(as.character(sims),paste(2012:(maxYr+1))))
    catch[,as.character(2012)] <- rep(251809,length(sims))   #The set catch for 2012

    surveyYrs=c(2012,seq(2013,maxYr,surveyPer))

    for(i in sims) {    #simulation of MCMC sample (first one is burn-in)
        cat("Starting sim",i,"\n")
        flush.console()
        #create directory and copy files to it
        dir.create(paste("sim",i,sep=""),showWarnings=F)
        file.copy(c("2012_hake_data.ss","2012_hake_control.ss","forecast.ss","forecastAssess.ss","starter.ss","SS3.exe","ss3.psv","wtatage.ss"),paste("sim",i,sep=""))
        setwd(paste("sim",i,sep=""))
        file.copy("2012_hake_data.ss","2012_hake_data.2011.ss")
        file.copy("2012_hake_control.ss","2012_hake_control.2011.ss")
        starter <- SS_readstarter(verbose=F)
        #create starter file with the thinning starting at the current sim (to save a little time)
        starter$MCMCthin <- i
        SS_writestarter(starter,overwrite=T,verbose=F)
        #setup forecast.ss to have no catch to get beginning of year 2012
        forecastSS <- readLines("forecast.ss")
        ind <- grep("Number of forecast catch levels to input",forecastSS)
        tmp <- forecastSS[ind]
        tmp <- strsplit(tmp,"#")[[1]]
        tmp[1] <- 0
        forecastSS[ind] <- paste(tmp,collapse=" #")
        ind <- grep("basis for input Fcast catch",forecastSS)
        forecastSS <- forecastSS[1:ind]
        forecastSS <- c(forecastSS,999)
        writeLines(forecastSS,"forecast.ss")
    
        #Run mceval just to make sure that it gets proper values for beginning of year 2012
        tmp <- shell("ss3.exe -mceval",intern=T)
        M <- params$NatM_p_1_Fem_GP_1[params$Iter==i]
        qSurv <- params$Q_2[params$Iter==i]
        #get selectivity for each fleet for this particular MCMC sample
        s1 <- selex[[1]][selex[[1]]$mceval==i,6:26]  #fishery selex
        s2 <- selex[[2]][selex[[2]]$mceval==i,6:26]  #survey selex
        
        #Loop over years (start in 2012 since 2012 catch is already set, but need 2012 data)
        for(yr in assessYrs) {
            dat <- SS_readdat("2012_hake_data.ss",verbose=F)
            natage <- read.table("posterior_natage.sso",header=T,nrows=19)  #only need the first 19 rows because we are working with only a particular sim (saves memory and time)
            natage <- natage[natage$mceval==i & natage$Yr==yr,-(1:6)]
            if(nrow(natage)!=1) stop("natage not read in correctly in sim",i,"and year",yr,"\n")
            batage <- natage*wtAtAge
            if(yr %in% surveyYrs) { #Simulate new survey data for year and add to data file
                dat <- simulateSurveyIndex.fn(dat,yr,batage,s2,qSurv,M,surveySElog,surveySEtotal)        
                #simulate age comp data
                dat <- simulateAgeComp.fn(dat,yr,natage,s2,M,surveyNage,surveyNageAdj,fleet=2,ageingError)
            }
            #simulate fishery age comp data for year
            if(catch[as.character(i),as.character(yr)] > 0) {
                dat <- simulateAgeComp.fn(dat,yr,natage,s1,M,commNage,commNageAdj,fleet=1,ageingError)
            }
            #add year catch to data file
            dat$endyr <- yr   #dat$endyr+1
            if(dat$endyr!=yr) stop("End year in dat file does not match the year being simulated\n")
            dat$N_catch <- dat$N_catch+1
            dat$catch <- rbind(dat$catch,c(catch[as.character(i),as.character(yr)],yr,1))
            SS_writedat(dat, outfile="2012_hake_data.ss", overwrite = TRUE, verbose = F)
            file.copy("2012_hake_data.ss",paste("2012_hake_data",yr,"ss",sep="."))
            #add year to recruitment dev vector in control file
            controlFile.fn()
            file.copy("2012_hake_control.ss",paste("2012_hake_control",yr,"ss",sep="."))
            #assess (in a separate folder) to determine catch for next year
            dir.create(paste("assess",yr,sep=""))
            file.copy(c("2012_hake_data.ss","2012_hake_control.ss","forecastAssess.ss","starter.ss","SS3.exe","wtatage.ss"),paste("assess",yr,sep=""))
              #I may want to decrese the nForecastYrs so that it only forecasts the one year that I need. I'm not sure how much time this will save
            setwd(paste("assess",yr,sep=""))
            file.rename("forecastAssess.ss","forecast.ss")
            #create the par file to start from (but not in start year
            if(yr != assessYrs[1]) {
                parFile.fn(yr)
            }
            #file.rename(paste("2012_hake_data",yr,"ss",sep="."),"2012_hake_data.ss")
            finalGradient[as.character(i),as.character(yr)] <- assess.fn(i,yr,assessCmd,convCrit)
            #done assessment, now get predicted catch (set catch less than 100 to 0
            forecastOut <- readLines("Forecast-report.sso")
            ind <- grep("FORECAST:_With_F_to_match_adjusted_catch",forecastOut)
            forecastOut <- read.table("Forecast-report.sso",skip=(ind[1]),nrows=1,header=T)
            catch[as.character(i),as.character(yr+1)] <- floor(forecastOut$Total_Catch)
            #set min catch. Do not put at 0 because unexpected things will happen
            if(catch[as.character(i),as.character(yr+1)] < minCatch) catch[as.character(i),as.character(yr+1)] <- 10
            #clean up to save space
            removeAssessJunk.fn()
            setwd("..")
            #add catch for this year to forecast file to generate true population at start of next year
            forecastSS <- readLines("forecast.ss")
            ind <- grep("Number of forecast catch levels to input",forecastSS)
            tmp <- forecastSS[ind]
            tmp <- strsplit(tmp,"#")[[1]]
            tmp[1] <- yr+1-2012
            forecastSS[ind] <- paste(tmp,collapse=" #")
            forecastSS <- forecastSS[1:(length(forecastSS)-1)]
            forecastSS <- c(forecastSS,paste(yr,"1 1",catch[as.character(i),as.character(yr)]),999)
            writeLines(forecastSS,"forecast.ss")
            #runmceval to get next years beginning biomass
            #make sure to use propoer data and control files, so that it reads in forecast catches
            file.copy("2012_hake_data.2011.ss","2012_hake_data.ss",overwrite=T)
            file.copy("2012_hake_control.2011.ss","2012_hake_control.ss",overwrite=T)
            tmp <- shell("ss3.exe -mceval",intern=T)
            #copy files back so that it can read in next year
            file.copy(paste("2012_hake_data",yr,"ss",sep="."),"2012_hake_data.ss",overwrite=T)
            file.copy(paste("2012_hake_control",yr,"ss",sep="."),"2012_hake_control.ss",overwrite=T)
        }
        file.remove(dir()[grep("ss_new",dir())])
        file.remove("Report.sso","CompReport.sso","Forecast-report.sso","SS3.exe","covar.sso","admodel.dep","admodel.cov","admodel.hes","echoinput.sso","SIS_table.sso","ParmTrace.sso","runnumber.ss","checkup.sso","CumReport.sso","eigv.rpt","fmin.log","rebuild.sso","variance")
        setwd("..")
        write.csv(catch,file="catch.csv")
    }
    setwd(origWD)
}


getSimResults.fn <- function(dirName,simnums=NULL){
    #gets the results from the operating model (truth) for the cases that did not do assessments (i.e., no fishing and omniscient)
    #simnums can either be a vector of numbers indicating the tow numbers of the MCMC to use. This is in case only a subset of the MCMC samples were used in other simulations.
    #it can take on three different meanings
    ###  NULL: It uses all rows from the MCMC
    ###  numeric vector: it uses only those row numbers. i.e., c(4,8,27)
    ###  character: it defines a directory to look in to get the sim numbers (the directory contains folders that begin with "sim")

    if(!(substring(dirName,nchar(dirName))=="/" | substring(dirName,nchar(dirName))=="\\")) {
        dirName <- paste(dirName,"/",sep="")
    }
    
    j <- tryCatch({
        suppressWarnings(
            x <- read.table(paste(dirName,"derived_posteriors.sso",sep=""),header=T)
        )
    }, error = function(err){
        cat(paste("Error - cannot find derived_posteriors.sso in",dirName,"\n\n"))
        stop()
    })
    rownames(x) <- as.character(x$Iter)

    if(is.character(simnums)) {
        if(!(substring(simnums,nchar(simnums))=="/" | substring(simnums,nchar(simnums))=="\\")) {
            simnums <- paste(simnums,"/",sep="")
        }

        dirL <- dir(simnums)

        suppressWarnings(
            simnums <- sort(as.numeric(substring(dirL[grep("sim",dirL)],4)))
        )

    }
    if(is.numeric(simnums)) {
        simnums <- as.character(simnums)
    }    
    if(is.null(simnums)) {
        simnums <- as.character(x$Iter)
    }

    SPBmat <- x[simnums,grep("SPB",names(x))]
    SPRmat <- x[simnums,grep("SPRratio",names(x))]
    Fmat   <- x[simnums,grep("F_",names(x))]
    catchesReal <- x[simnums,grep("ForeCatch",names(x))]
    F_Btgt <- x[simnums,"Fstd_Btgt"]
    F_SPRtgt <- x[simnums,"Fstd_SPRtgt"]
    names(F_Btgt) <- names(F_SPRtgt) <- simnums
    deplMat <- t(apply(SPBmat,1,function(x){x/x[1]}))[,-c(1,2)]
    Rmat   <- x[simnums,grep("Recr_",names(x))]

    return(list(SSB=SPBmat,SPRratio=SPRmat,Fvalue=Fmat,depl=deplMat,catchesReal=catchesReal,F_Btgt=F_Btgt,F_SPRtgt=F_SPRtgt,Rmat=Rmat))
}



getSimResultsAssess.fn <- function(dirName,assessment=F,assessYrs=2012:2030,recrYrs=1975:2030){
    #gets the results from the operating model (truth) for the cases that actually do assessments (i.e., annual and biennial surveys)
    
    if(!(substring(dirName,nchar(dirName))=="/" | substring(dirName,nchar(dirName))=="\\")) {
        dirName <- paste(dirName,"/",sep="")
    }

    dirList <- dir(dirName)
    suppressWarnings(
        simnums <- sort(as.numeric(substring(dirList[grep("sim",dirList)],4)))
    )

    j <- tryCatch({
        suppressWarnings(
            x <- read.table(paste(dirName,"sim",simnums[1],"/derived_posteriors.sso",sep=""),header=T,nrows=1)
        )
    }, error = function(err){
        cat(paste("Error - cannot find any simulations in",dirName,"\n\n"))
        stop()
    })

    SPBmat <- matrix(NA,ncol=length(grep("SPB",names(x))),nrow=length(simnums),dimnames=list(simnums,names(x)[grep("SPB",names(x))]))
    SPRmat <- matrix(NA,ncol=length(grep("SPRratio",names(x))),nrow=length(simnums),dimnames=list(simnums,names(x)[grep("SPRratio",names(x))]))
    Fmat   <- matrix(NA,ncol=length(grep("F_",names(x))),nrow=length(simnums),dimnames=list(simnums,names(x)[grep("F_",names(x))]))
    Rmat <- matrix(NA,ncol=length(grep("Recr",names(x))),nrow=length(simnums),dimnames=list(simnums,names(x)[grep("Recr_",names(x))]))
    catchesReal <- matrix(NA,ncol=length(grep("ForeCatch",names(x))),nrow=length(simnums),dimnames=list(simnums,names(x)[grep("ForeCatch",names(x))]))
    assSSB <- matrix(NA,ncol=length(assessYrs),nrow=length(simnums),dimnames=list(simnums,as.character(assessYrs+1)))
    assDepl<- matrix(NA,ncol=length(assessYrs),nrow=length(simnums),dimnames=list(simnums,as.character(assessYrs+1)))
    if(!is.null(recrYrs)) {
        assRecr<- array(NA,dim=c(length(recrYrs),length(assessYrs),length(simnums)),dimnames=list(paste("recr",recrYrs,sep=""),paste("ass",assessYrs+1,sep=""),simnums))
    }
    F_Btgt <- rep(NA,length(simnums))
    F_SPRtgt <- rep(NA,length(simnums))
    names(F_Btgt) <- names(F_SPRtgt) <- simnums

    for(i in simnums) {
        derPost <- read.table(paste(dirName,"sim",i,"/derived_posteriors.sso",sep=""),header=T,nrows=1)
        if(derPost[1,1] != i) {stop(paste("Iter number did not match sim number in sim",i))}

        SPBmat[as.character(i),] <- unlist(derPost[,grep("SPB",names(x))])
        SPRmat[as.character(i),] <- unlist(derPost[,grep("SPRratio",names(x))])
        catchesReal[as.character(i),] <- unlist(derPost[,grep("ForeCatch",names(x))])
        Fmat[as.character(i),] <- unlist(derPost[,grep("F_",names(x))])
        F_Btgt[as.character(i)] <- derPost[,"Fstd_Btgt"]
        F_SPRtgt[as.character(i)] <- derPost[,"Fstd_SPRtgt"]
        Rmat[as.character(i),] <- unlist(derPost[,grep("Recr_",names(x))])

        if(assessment) {
            xx <- getAssessmentResults.fn(paste(dirName,"sim",i,"/",sep=""),assessYrs,recrYrs=recrYrs)
            assSSB[as.character(i),] <- xx$SSB
            assDepl[as.character(i),] <- xx$depl
            if(!is.null(recrYrs)) {
                assRecr[,,as.character(i)] <- xx$recr
            }
        }
    }

    deplMat <- t(apply(SPBmat,1,function(x){x/x[1]}))[,-c(1,2)]

    catches <- read.csv(paste(dirName,"catch.csv",sep=""))
    rownames(catches) <- as.character(catches[,1])
    catches <- catches[,-1]

    if(assessment) {
        if(!is.null(recrYrs)) {
            return(list(SSB=SPBmat,SPRratio=SPRmat,Fvalue=Fmat,depl=deplMat,catches=catches,catchesReal=catchesReal,Rmat=Rmat,F_Btgt=F_Btgt,F_SPRtgt=F_SPRtgt,assessSSB=assSSB,assessDepl=assDepl,assessRecr=assRecr))
        }
        if(is.null(recrYrs)) {
            return(list(SSB=SPBmat,SPRratio=SPRmat,Fvalue=Fmat,depl=deplMat,catches=catches,catchesReal=catchesReal,Rmat=Rmat,F_Btgt=F_Btgt,F_SPRtgt=F_SPRtgt,assessSSB=assSSB,assessDepl=assDepl))
        }
    }
    
    return(list(SSB=SPBmat,SPRratio=SPRmat,Fvalue=Fmat,depl=deplMat,catches=catches,catchesReal=catchesReal,Rmat=Rmat,F_Btgt=F_Btgt,F_SPRtgt=F_SPRtgt))
    
    #catches contains the catches set by the TAC
    #catchesReal contains the catches realized (sometimes the model cannot take the set catch)
}


getAssessmentResults.fn <- function(direc,yrs,recrYrs=NULL,nLines=800) {
    #gets assessment result for one simulation (i.e, sim4)
    #predicted biomass for start of year after assessment is done (what catch is based on)
    brat <- ssb <- NULL
    if(!is.null(recrYrs)) {recr <- matrix(NA,nrow=length(recrYrs),ncol=length(yrs),dimnames=list(recrYrs,yrs))}
    if(is.null(recrYrs)) {recr <- NULL}
    for(i in yrs) {
        out <- readLines(paste(direc,"assess",i,"/Report.sso",sep=""),n=nLines)
        tmp <- out[grep(paste("Bratio",i+1,sep="_"),out)]   #get next year bratio since that is what is used to determine catch
        brat <- c(brat,as.numeric(strsplit(tmp," ")[[1]][3]))
        tmp <- out[grep(paste("SPB",i+1,sep="_"),out)]   #get next year bratio since that is what is used to determine catch
        ssb <- c(ssb,as.numeric(strsplit(tmp," ")[[1]][3]))
        if(!is.null(recrYrs)) {
            tmp <- out[grep(" Recr_",out)]   #get next year bratio since that is what is used to determine catch
            tmp <- as.data.frame(strsplit(tmp," "))
            tmpYrs <- suppressWarnings(as.numeric(substr(unlist(tmp[2,]),6,9)))
            indYrs <- tmpYrs %in% recrYrs
            recr[,as.character(i)] <- as.numeric(as.character(unlist(tmp[3,indYrs])))
        }
    }
    list(SSB=ssb,depl=brat,recr=recr)
}


doHakeMSE2.fn <- function(theDir,maxYr=2030,
                         surveyPer,surveySElog=0.085,surveySEtotal=0.42,surveyNage=65,surveyNageAdj=1.0,
                         commNage=800,commNageAdj=0.12,minCatch=100,
                         wtAtAge=c(0.03,0.0885,0.2562,0.3799,0.4913,0.5434,0.5906,0.662,0.7215,0.791,0.8629,0.9315,0.9681,1.0751,1.0016,1.0202,1.0202,1.0202,1.0202,1.0202,1.0202),
                         nSims=200,convCrit=0.1,assessCmd="ss3 -nohess",theSeed,origSeed=716)
{
 #samples a different set than the origianl 200 sims
    origWD <- getwd()
    set.seed(origSeed)
    x <- 2:1000
    origSims <- sort(sample(x,nSims))
    cat("Original Sim Numbers:",origSims,"\n")


    set.seed(theSeed)

    setwd(theDir)

    selex <- read.table("posterior_mse.sso",header=T)
    selex <- split(selex,selex$Fleet)

    assessYrs <- 2012:maxYr
    sims <- sort(sample(x[!x%in%origSims],nSims))
    cat("Sim Numbers:",sims,"\n")
    finalGradient <- matrix(NA,nrow=nSims,ncol=length(assessYrs),dimnames=list(as.character(sims),as.character(assessYrs)))

    params <- read.table("posteriors.sso",header=T)  #posterior parameter should not change with each mceval, so read it in now
    catch <- matrix(NA,nrow=nSims,ncol=maxYr-2012+2,dimnames=list(as.character(sims),paste(2012:(maxYr+1))))
    catch[,as.character(2012)] <- rep(251809,nSims)   #The set catch for 2012

    surveyYrs=c(2012,seq(2013,maxYr,surveyPer))

    for(i in sims) {    #simulation of MCMC sample (first one is burn-in)
        cat("Starting sim",i,"\n")
        flush.console()
        #create directory and copy files to it
        dir.create(paste("sim",i,sep=""),showWarnings=F)
        file.copy(c("2012_hake_data.ss","2012_hake_control.ss","forecast.ss","starter.ss","SS3.exe","ss3.psv","wtatage.ss"),paste("sim",i,sep=""))
        setwd(paste("sim",i,sep=""))
        file.copy("2012_hake_data.ss","2012_hake_data.2011.ss")
        file.copy("2012_hake_control.ss","2012_hake_control.2011.ss")
        starter <- SS_readstarter(verbose=F)
        #create starter file with the thinning starting at the current sim (to save a little time)
        starter$MCMCthin <- i
        SS_writestarter(starter,overwrite=T,verbose=F)
        #setup forecast.ss to have no catch to get beginning of year 2012
        forecastSS <- readLines("forecast.ss")
        ind <- grep("Number of forecast catch levels to input",forecastSS)
        tmp <- forecastSS[ind]
        tmp <- strsplit(tmp,"#")[[1]]
        tmp[1] <- 0
        forecastSS[ind] <- paste(tmp,collapse=" #")
        ind <- grep("basis for input Fcast catch",forecastSS)
        forecastSS <- forecastSS[1:ind]
        forecastSS <- c(forecastSS,999)
        writeLines(forecastSS,"forecast.ss")
    
        #Run mceval just to make sure that it gets proper values for beginning of year 2012
        tmp <- shell("ss3.exe -mceval",intern=T)
        M <- params$NatM_p_1_Fem_GP_1[params$Iter==i]
        qSurv <- params$Q_2[params$Iter==i]
        #get selectivity for each fleet for this particular MCMC sample
        s1 <- selex[[1]][selex[[1]]$mceval==i,6:26]  #fishery selex
        s2 <- selex[[2]][selex[[2]]$mceval==i,6:26]  #survey selex
        
        #Loop over years (start in 2012 since 2012 catch is already set, but need 2012 data)
        for(yr in assessYrs) {
            dat <- SS_readdat("2012_hake_data.ss",verbose=F)
            natage <- read.table("posterior_natage.sso",header=T,nrows=19)  #only need the first 19 rows because we are working with only a particular sim (saves memory and time)
            natage <- natage[natage$mceval==i & natage$Yr==yr,-(1:6)]
            if(nrow(natage)!=1) stop("natage not read in correctly in sim",i,"and year",yr,"\n")
            batage <- natage*wtAtAge
            if(yr %in% surveyYrs) { #Simulate new survey data for year and add to data file
                dat <- simulateSurveyIndex.fn(dat,yr,batage,s2,qSurv,M,surveySElog,surveySEtotal)        
                #simulate age comp data
                dat <- simulateAgeComp.fn(dat,yr,natage,s2,M,surveyNage,surveyNageAdj,fleet=2)
            }
            #simulate fishery age comp data for year
            if(catch[as.character(i),as.character(yr)] > 0) {
                dat <- simulateAgeComp.fn(dat,yr,natage,s1,M,commNage,commNageAdj,fleet=1)
            }
            #add year catch to data file
            dat$endyr <- yr   #dat$endyr+1
            if(dat$endyr!=yr) stop("End year in dat file does not match the year being simulated\n")
            dat$N_catch <- dat$N_catch+1
            dat$catch <- rbind(dat$catch,c(catch[as.character(i),as.character(yr)],yr,1))
            SS_writedat(dat, outfile="2012_hake_data.ss", overwrite = TRUE, verbose = F)
            file.copy("2012_hake_data.ss",paste("2012_hake_data",yr,"ss",sep="."))
            #add year to recruitment dev vector in control file
            controlFile.fn()
            file.copy("2012_hake_control.ss",paste("2012_hake_control",yr,"ss",sep="."))
            #assess (in a separate folder) to determine catch for next year
            dir.create(paste("assess",yr,sep=""))
            file.copy(c("2012_hake_data.ss","2012_hake_control.ss","forecast.ss","starter.ss","SS3.exe","wtatage.ss"),paste("assess",yr,sep=""))
              #I may want to decrese the nForecastYrs so that it only forecasts the one year that I need. I'm not sure how much time this will save
            setwd(paste("assess",yr,sep=""))
            #create the par file to start from (but not in start year
            if(yr != assessYrs[1]) {
                parFile.fn(yr)
            }
            #file.rename(paste("2012_hake_data",yr,"ss",sep="."),"2012_hake_data.ss")
            finalGradient[as.character(i),as.character(yr)] <- assess.fn(i,yr,assessCmd,convCrit)
            #done assessment, now get predicted catch (set catch less than 100 to 0
            forecastOut <- readLines("Forecast-report.sso")
            ind <- grep("FORECAST:_With_F_to_match_adjusted_catch",forecastOut)
            forecastOut <- read.table("Forecast-report.sso",skip=(ind[1]),nrows=1,header=T)
            catch[as.character(i),as.character(yr+1)] <- floor(forecastOut$Total_Catch)
            #set min catch. Do not put at 0 because unexpected things will happen
            if(catch[as.character(i),as.character(yr+1)] < minCatch) catch[as.character(i),as.character(yr+1)] <- 10
            #clean up to save space
            removeAssessJunk.fn()
            setwd("..")
            #add catch for this year to forecast file to generate true population at start of next year
            forecastSS <- readLines("forecast.ss")
            ind <- grep("Number of forecast catch levels to input",forecastSS)
            tmp <- forecastSS[ind]
            tmp <- strsplit(tmp,"#")[[1]]
            tmp[1] <- yr+1-2012
            forecastSS[ind] <- paste(tmp,collapse=" #")
            forecastSS <- forecastSS[1:(length(forecastSS)-1)]
            forecastSS <- c(forecastSS,paste(yr,"1 1",catch[as.character(i),as.character(yr)]),999)
            writeLines(forecastSS,"forecast.ss")
            #runmceval to get next years beginning biomass
            #make sure to use propoer data and control files, so that it reads in forecast catches
            file.copy("2012_hake_data.2011.ss","2012_hake_data.ss",overwrite=T)
            file.copy("2012_hake_control.2011.ss","2012_hake_control.ss",overwrite=T)
            tmp <- shell("ss3.exe -mceval",intern=T)
            #copy files back so that it can read in next year
            file.copy(paste("2012_hake_data",yr,"ss",sep="."),"2012_hake_data.ss",overwrite=T)
            file.copy(paste("2012_hake_control",yr,"ss",sep="."),"2012_hake_control.ss",overwrite=T)
        }
        file.remove(dir()[grep("ss_new",dir())])
        file.remove("Report.sso","CompReport.sso","Forecast-report.sso","SS3.exe","covar.sso","admodel.dep","admodel.cov","admodel.hes","echoinput.sso","SIS_table.sso","ParmTrace.sso","runnumber.ss","checkup.sso","CumReport.sso","eigv.rpt","fmin.log","rebuild.sso","variance")
        setwd("..")
        write.csv(catch,file="catch.csv")
    }
    setwd(origWD)
}

randWalkSelex.fn <- function(pars) {
    #calculates the selectivity from the random walk parameters in SS (option 17)
    #-1000 means to set equal to 0
    #assumes that this is all pars from age 0 to max age
    
    logS <- rep(NA,length(pars))
    logS[1] <- 0 #first value is never estimated (age 0)
    for(a in 2:length(pars)) {
        ifelse(pars[a] == -1000, logS[a] <- 0, logS[a] <- logS[a-1]+pars[a])
    }
    selex <- exp(logS-max(logS))
    selex[pars== -1000] <- 0
    return(selex)
}

randWalkSelex.fn(c(-1000,-1000,0,0.295698,0.170968,-0.441426,0.564627,0,0,0,0,0,0,0))
randWalkSelex.fn(c(-1000,-1000,0,-0.0605365,0.141217,0.474321,-0.00714344,0,0,0,0,0,0,0))


randWalkSelexPars.fn <- function(selex) {
    #determines the parameter values of the random walk selectivity in SS (option 17) given the selectivity vector
    #it assumes the selectivity vector begins at age 0
    logS <- log(selex)
    maxLogS <- max(logS)
    n <- length(selex)
    
    p <- logS[2:n] - logS[1:(n-1)] - maxLogS
    p[is.infinite(p)] <- -1000
    return(p)
}

0   0.554401    0.74515 0.884086    0.568572    1   1   1   1   1   1.00E+00    1   1.00E+00    1.00E+00    1

0   0.574071    0.54035 0.622307    1   0.992882    0.992882    0.992882    0.992882    0.992882    0.992882    0.992882    0.992882    0.992882    0.992882
