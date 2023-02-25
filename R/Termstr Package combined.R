## Bond pricing function

bond_prices <- function(method="ns", beta, m, cf, lambda,mythology = TRUE){
  # calculate spot rates
  spot_rates <- spotrates(method,beta,m,lambda)/100
  
  # replace NaNs by zeros
  spot_rates[is.nan(spot_rates)] <- 0        
  
  if (mythology){# calculate discount factors
    discount_factors <- exp(-m*spot_rates)
  } else {
    discount_factors <- 1/(1+spot_rates)^m
  }
  
  # calculate bond prices
  bond_prices <- apply(cf*discount_factors, 2, sum)  
  
  # return spot rates, discount factors and bond prices
  return(list(spot_rates=spot_rates,discount_factors=discount_factors,
              bond_prices=bond_prices))
}  


## Bond pricing function

bond_prices2 <- function(method="ns", beta, m, cf, lambda){
  # calculate spot rates
  spot_rates <- spotrates(method,beta,m,lambda)/100
  
  # replace NaNs by zeros
  spot_rates[is.nan(spot_rates)] <- 0        
  
  # calculate discount factors
  discount_factors <- exp(-m*spot_rates)
  
  # calculate bond prices
  bond_prices <- apply(cf*discount_factors, 2, sum)  
  
  # return spot rates, discount factors and bond prices
  return(list(spot_rates=spot_rates,discount_factors=discount_factors,
              bond_prices=bond_prices))
}  

## Calculation of bond yields

bond_yields <- function(cashflows, m, m2, freque = 2, searchint=c(-1,9), tol=1e-10, mythology = TRUE, Annulized=TRUE) {
  # convert input data to matrices if necessary
  if (!is.matrix(cashflows))
    cashflows <- as.matrix(cashflows)
  if (!is.matrix(m))
    m <- as.matrix(m)
  if (!is.matrix(m2))
    m2 <- as.matrix(m2)
  # create empty bond yields matrix in appropriate size
  bondyields <- matrix(0, nrow=ncol(cashflows), ncol=2)
  # put maturity of every bond into first column of bond yields matrix
  bondyields[,1] <- apply(m, 2, max)
  # traverse list of bonds
  
  if (mythology) {
    for (i in seq_len(ncol(cashflows))) {
      # present value of cash flows for root finding 
      pvcashflows <- function(y) {
        t(cashflows[,i])%*%exp(-m[,i]*y)
      }
      # calculate bond yields
      bondyields[i,2] <- uniroot(pvcashflows, searchint, tol = tol,maxiter=3000)$root 
    }
  } else {
    for (i in seq_len(ncol(cashflows))) {
      # present value of cash flows for root finding 
      pvcashflows2 <- function(y) {
        t(cashflows[,i])%*%(1/((1+y/freque)^m2[,i]))
      }
      # calculate bond yields
      bondyields[i,2] <- uniroot(pvcashflows2, searchint, tol = tol,maxiter=3000)$root 
      if (Annulized){
        bondyields[i,2] <- ((1+bondyields[i,2]/2)^2)-1
      }
    }
  }
  
  # return calculated bond yields matrix
  rownames(bondyields) <- colnames(cashflows)
  colnames(bondyields) <- c("Maturity","Yield")
  bondyields
}

## Calculation of the duration,modified duration and duration based weights 
duration <- function (cf_p,m_p,y,m_p_2,freque,mythology,wmethod){
  
  y <- matrix(rep(y,nrow(m_p)),ncol=ncol(m_p),byrow=TRUE)
  
  if (mythology){
    # cont compounded
    # mac cauly duration (Definition wie Filipovic)
    d <- apply(cf_p*m_p*exp(-y*m_p),2,sum)/-cf_p[1,]
    # See Literature ModD = MacD
    md <- d
  } else {
    # Macaulay duration nach TradeWeb/cBonds
    d <- apply(m_p_2*cf_p/((1+y/freque)^m_p_2),2,sum)/(-cf_p[1,]*freque)
    # Modified Duration
    md <- d/(1+y[1,])
  }
  
  # Weights
  if (wmethod =="rf") {
    omega <- (1/d)/sum(1/d)
    text1 <- paste("Weights (calculated after)",wmethod,sep=" ")
  } else if (wmethod =="Yallup"){
    omega <- (100/(-cf_p[1,]*md))^2
    text1 <- paste("Weights (calculated after)",wmethod,sep=" ")
  } else if (wmethod =="Anderson"){
    #Implies Prices near nominal
    omega <- 1/(d^2)
    text1 <- paste("Weights (calculated after)",wmethod,sep=" ")
  } else {
    omega <- 1/d
    text1 <- "Weights"
  }
  
  dur <- cbind(d,md,omega)
  colnames(dur) <- c("Macaulay Duration","Modified duration", text1)
  dur
  
}

## Loss function: mean squared (weighted error)
loss_function <- function(p,phat,omega) {
  sum(omega*((p-phat)^2))
}


## Limit a bond data set to a certain maturity range 
maturity_range <- function(bonddata,lower,upper,version=0,SettlementCheck=FALSE) {
  m <- lapply(bonddata,version=version,SettlementCheck=SettlementCheck,create_maturities_matrix)
  colmax <- function(m) apply(m,2,max)
  m_max <- lapply(m,colmax)
  
  bonds_in_range <- function(group) names(group[which(
    group>lower & group<upper)])
  
  # list with ISINs of bonds in range
  isins_range <- lapply(m_max,bonds_in_range)
  index_set <- which(unlist(lapply(isins_range,length))>0)
  
  # list with positions of bonds in isins_range
  isins_range_pos <- list()
  
  for (i in seq_along(bonddata)) {
    isins_range_pos[[i]] <- which(bonddata[[i]]$ISIN %in% 
                                    isins_range[[i]])
  }
  names(isins_range_pos) <- names(bonddata)
  
  # first part of bonddata for filtering
  
  N <- which(names(bonddata[[1]]) %in% c("ISIN","MATURITYDATE","ISSUEDATE","COUPONRATE","PRICE","ACCRUED"))
  
  first <- function(lst) lst[N]
  filtered <- lapply(bonddata,first)
  
  bonddata_range <- list()
  for (i in seq_along(bonddata)) {
    bonddata_range[[i]] <- as.list(as.data.frame(filtered[[i]])
                                   [isins_range_pos[[i]],])
    # convert to character                       
    bonddata_range[[i]][["ISIN"]] <- as.character(bonddata_range[[i]][["ISIN"]])
    bonddata_range[[i]][["MATURITYDATE"]] <- as.character(bonddata_range[[i]][["MATURITYDATE"]])
    bonddata_range[[i]][["ISSUEDATE"]] <- as.character(bonddata_range[[i]][["ISSUEDATE"]])
  }
  names(bonddata_range) <- names(bonddata)
  
  # list with positions of cashflows in isins_range
  
  isins_range_pos <- list()
  for (i in seq_along(bonddata)) {
    isins_range_pos[[i]] <- which(bonddata[[i]][["CASHFLOWS"]]
                                  [["ISIN"]]%in%isins_range[[i]])
  }
  names(isins_range_pos) <- names(bonddata)
  
  for (i in seq_along(bonddata)) {
    CASHFLOWS <- as.list(as.data.frame(bonddata[[i]]
                                       [["CASHFLOWS"]])[isins_range_pos[[i]],])
    CASHFLOWS$ISIN <- as.character(CASHFLOWS$ISIN)
    CASHFLOWS$DATE <- as.character(CASHFLOWS$DATE)             
    bonddata_range[[i]][["CASHFLOWS"]] <- list()
    bonddata_range[[i]][["CASHFLOWS"]] <- CASHFLOWS
  }
  names(bonddata_range) <- names(bonddata)
  bonddata_range
  
  # add TODAY from bonddata
  for (i in seq_along(bonddata)) {
    bonddata_range[[i]][["TODAY"]] <- bonddata[[i]][["TODAY"]]
    bonddata_range[[i]][["SETTLEMENT DAY"]] <- bonddata[[i]][["SETTLEMENT DAY"]] 
  }
  
  # delete countries where no bonds are available
  bonddata_range <- bonddata_range[index_set]
  bonddata_range
}

rm_bond <- function(bonddata, group, ISIN) UseMethod("rm_bond")

## remove bonds from a static bonddata set  
rm_bond.couponbonds <- function(bonddata, group, ISIN){
  cf_isin_index <- which(bonddata[[group]]$CASHFLOWS$ISIN %in% ISIN)
  isin_index <- which(bonddata[[group]]$ISIN %in% ISIN)	
  bonddata[[group]]$ISIN <-  bonddata[[group]]$ISIN[-isin_index]
  bonddata[[group]]$MATURITYDATE <- bonddata[[group]]$MATURITYDATE[-isin_index]
  bonddata[[group]]$ISSUEDATE <- bonddata[[group]]$ISSUEDATE[-isin_index] 
  #  bonddata[[group]]$STARTDATE <- bonddata[[group]]$STARTDATE[-isin_index] 
  bonddata[[group]]$COUPONRATE <- bonddata[[group]]$COUPONRATE[-isin_index]
  bonddata[[group]]$PRICE <- bonddata[[group]]$PRICE[-isin_index]
  bonddata[[group]]$ACCRUED <- bonddata[[group]]$ACCRUED[-isin_index]
  bonddata[[group]]$CASHFLOWS$ISIN <- bonddata[[group]]$CASHFLOWS$ISIN[-cf_isin_index]
  bonddata[[group]]$CASHFLOWS$CF <- bonddata[[group]]$CASHFLOWS$CF[-cf_isin_index]
  bonddata[[group]]$CASHFLOWS$DATE <- bonddata[[group]]$CASHFLOWS$DATE[-cf_isin_index]
  
  class(bonddata) <- "couponbonds"
  bonddata
}

## remove bonds from a dynamic bonddata set 
rm_bond.dyncouponbonds <- function(bonddata, group, ISIN) {
  for (i in seq(length(bonddata))) {
    cf_isin_index <- which(bonddata[[i]][[group]]$CASHFLOWS$ISIN %in% ISIN)
    isin_index <- which(bonddata[[i]][[group]]$ISIN %in% ISIN)
    bonddata[[i]][[group]]$ISIN <- bonddata[[i]][[group]]$ISIN[-isin_index]
    bonddata[[i]][[group]]$MATURITYDATE <- bonddata[[i]][[group]]$MATURITYDATE[-isin_index]
    bonddata[[i]][[group]]$ISSUEDATE <- bonddata[[i]][[group]]$ISSUEDATE[-isin_index] #Geändert von Startdate zu Issuedate
    #    bonddata[[i]][[group]]$STARTDATE <- bonddata[[i]][[group]]$STARTDATE[-isin_index]
    bonddata[[i]][[group]]$COUPONRATE <- bonddata[[i]][[group]]$COUPONRATE[-isin_index]
    bonddata[[i]][[group]]$PRICE <- bonddata[[i]][[group]]$PRICE[-isin_index]
    bonddata[[i]][[group]]$ACCRUED <- bonddata[[i]][[group]]$ACCRUED[-isin_index]
    bonddata[[i]][[group]]$CASHFLOWS$ISIN <- bonddata[[i]][[group]]$CASHFLOWS$ISIN[-cf_isin_index]
    bonddata[[i]][[group]]$CASHFLOWS$CF <- bonddata[[i]][[group]]$CASHFLOWS$CF[-cf_isin_index]
    bonddata[[i]][[group]]$CASHFLOWS$DATE <- bonddata[[i]][[group]]$CASHFLOWS$DATE[-cf_isin_index]
  }
  class(bonddata) <- "dyncouponbonds"
  bonddata
}



## preprocess a static bonddataset,i.e., sort data, calculate cashflows, maturity matrix, yields 
prepro_bond <- function(group,
                        bonddata,
                        matrange="all", version = 0, SettlementCheck=FALSE,freque=2, mythology=TRUE, Annulized=TRUE, Actual = 0,wmethod="rf"){
  
  # select given group from bonddata
  bonddata <- bonddata[group]
  
  # select data according to chosen maturity range
  if (length(matrange)==1) {bonddata <- bonddata
  }else{bonddata <- maturity_range(bonddata,matrange[1],matrange[2],version, SettlementCheck) } ###   Version + SettlementCheck Term`
  
  # number of groups 
  n_group <- length(bonddata) 
  
  # group sequence
  sgroup <- seq(n_group)
  
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)
  
  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE),
                 sgroup,SIMPLIFY=FALSE)
  
  # create maturities matrix
  m <- lapply(bonddata, version=version,SettlementCheck=SettlementCheck,create_maturities_matrix) ###   Version Term
  
  # create maturities matrix for calculation of quoted yield including zeros (2.5.4 / Tradeweb)  
  # m2 <- lapply(bonddata, actual=actual,SettlementCheck=SettlementCheck,maturities_matrix_1825)  
  m_p_2 <- mapply(function(k) maturities_matrix_1825(bonddata[[k]],include_price=TRUE,Actual = Actual,SettlementCheck=SettlementCheck),sgroup,SIMPLIFY=FALSE)
  
  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE,version = version,SettlementCheck=SettlementCheck),sgroup,SIMPLIFY=FALSE) 
  ###   Version + SettlementCheck Term
  
  # calculate dirty prices
  p <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)
  
  # extract accrued interest
  ac <- mapply(function(k) bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)
  
  # assign ISIN 
  for(k in sgroup) {
    names(p[[k]]) <- bonddata[[k]]$ISIN
    names(ac[[k]]) <- bonddata[[k]]$ISIN
    
  }
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  for (i in seq(length(cf))) {
    if (is.null(dim(cf[[i]]))){
      name1 <- names(cf[[i]])
      cf[[i]] <- matrix(cf[[i]], nrow = 1)
      colnames(cf[[i]]) <- name1
    }
  }
  
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  
  ##MATRIX immer hinschreiben 
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  for (i in seq(length(m))) {
    if (is.null(dim(m[[i]]))){
      name1 <- names(m[[i]])
      m[[i]] <- matrix(m[[i]], nrow = 1)
      colnames(m[[i]]) <- name1
    }
  }
  
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p_2 <- mapply(function(k) m_p_2[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  
  p <- mapply(function(k) p[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  ac <- mapply(function(k) ac[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]],m_p_2[[k]], freque=freque, mythology=mythology, Annulized=Annulized),
              sgroup,SIMPLIFY=FALSE)
  # calculate duration
  duration <- mapply(function(k) duration(cf_p[[k]],m_p[[k]],y[[k]][,2],m_p_2[[k]],freque,mythology,wmethod),
                     sgroup,SIMPLIFY=FALSE)
  
  res <- list(n_group=n_group,sgroup=sgroup,cf=cf,cf_p=cf_p,m=m,m_p=m_p,p=p,ac=ac,y=y,duration=duration,m_p_2=m_p_2,timestamp=bonddata[[1]]$TODAY, settlementstamp=bonddata[[1]]$`SETTLEMENT DAY`)
  res
}

## Bonddata postprocessing function                      

postpro_bond <- function(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,m_p_2,method,lambda,freque, mythology, Annulized){
  
  # theoretical bond prices with estimated parameters
  phat <- mapply(function(k) bond_prices(method,opt_result[[k]]$par,
                                         m[[k]],cf[[k]],lambda,mythology)$bond_prices,sgroup,SIMPLIFY=FALSE)
  
  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
  
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]],m_p_2[[k]],freque=freque, mythology=mythology, Annulized=Annulized),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"
  
  
  # maturity interval
  if(min(mapply(function(i) min(y[[i]][,1]), sgroup))>0.24){
    minvalue <- 0.20
  } else {
    minvalue<- round(min(mapply(function(i) min(y[[i]][,1]), sgroup)),2)
  }
  if(max(mapply(function(i) max(y[[i]][,1]), sgroup))<30){
    maxvalue <- 30
  } else {
    maxvalue<- ceiling(max(mapply(function(i) max(y[[i]][,1]), sgroup)))
  }
  
  t <- seq(minvalue,maxvalue,0.01)
  
  
  # calculate zero coupon yield curves  
  zcy_curves <- mapply(function(k) cbind(t,spotrates(method,opt_result[[k]]$par,t,lambda)/100),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
  
  # calculate spread curves              	 
  if(n_group != 1) {  
    s_curves <- mapply(function(k) cbind(t,(zcy_curves[[k]][,2] - zcy_curves[[1]][,2])),sgroup,
                       SIMPLIFY=FALSE)
  } else s_curves = "none"
  
  for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
  class(s_curves) <- "s_curves"
  
  # calculate extrapolation point                        
  expoints <- mapply(function(k) which(zcy_curves[[k]][,1] > 
                                         mapply(function(i) max(y[[i]][,1]), seq(n_group))[k])[1],sgroup, SIMPLIFY=FALSE )  
  
  # calculate forward rate curves 
  fwr_curves <-  mapply(function(k) cbind(t,forwardrates(method,opt_result[[k]]$par,t,lambda)/100),sgroup,SIMPLIFY=FALSE)                   
  
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"
  
  
  # calculate discount factor curves 
  if (mythology){
    df_curves <- mapply(function(k) cbind(zcy_curves[[k]][,1],exp(-zcy_curves[[k]][,1]*zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  } else {
    df_curves <- mapply(function(k) cbind(zcy_curves[[k]][,1],1/((1+zcy_curves[[k]][,2])^zcy_curves[[k]][,1])),sgroup,SIMPLIFY=FALSE)
  }
  
  for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
  class(df_curves) <- "df_curves"
  
  res <- list(phat=phat,perrors=perrors,yhat=yhat,yerrors=yerrors,t=t,zcy_curves=zcy_curves,
              s_curves=s_curves,expoints=expoints,fwr_curves=fwr_curves,df_curves=df_curves,opt_result=opt_result)
  res
  
}


## Cashflows matrix creation function 

create_cashflows_matrix <- function(group,include_price=FALSE) {
  
  
  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN),maxsum=100000)
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  
  # the number of rows of the matrix is the number of 
  # cashflows of the bond with the maximum of cashflows
  # all missing elements of the matrix are filled up with zeros 
  
  CASHFLOWMATRIX <-
    mapply(function(i) c(group$CASHFLOWS$CF[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  if (include_price == TRUE) {CASHFLOWMATRIX <- rbind(-(group[["PRICE"]] +
                                                          group[["ACCRUED"]] ),CASHFLOWMATRIX)}     
  
  #colnames(CASHFLOWMATRIX) <- group$ISIN
  
  
  if (is.null(dim(CASHFLOWMATRIX))){
    # CASHFLOWMATRIX <- rbind(CASHFLOWMATRIX[1:length(CASHFLOWMATRIX)],rep(0,n_of_bonds))
    CASHFLOWMATRIX <- matrix(CASHFLOWMATRIX, nrow = 1)
    #CASHFLOWMATRIX <- CASHFLOWMATRIX[-1,]
  } 
  
  colnames(CASHFLOWMATRIX) <- group$ISIN
  CASHFLOWMATRIX
}

## Maturities matrix creation function 
create_maturities_matrix <- function(group,include_price=FALSE, version=0, SettlementCheck=FALSE) {
  
  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN),maxsum=100000) # Anzahl der CF pro Bond
  n_of_bonds <- length(n_of_cf) #Anzahl der ISIN Nummern / Anzahl der Bonds
  max_cf <- max(n_of_cf) #Max Anzahl an kommenden CF 
  pos_cf <- c(0,cumsum(n_of_cf))
  
  # Usage of Actual/365 $$$$OLD CODE$$$$
  #year_diff <- as.numeric(difftime(as.Date(group$CASHFLOW$DATE),as.Date(group$TODAY),units="days"))/365 
  
  # Using FUN newdaycount to calculate Maturities $$$$NEW CODE$$$
  # Using Settlement Day instead of Today Date
  
  if (SettlementCheck) {
    year_diff <- newdaycount(as.Date(group$`SETTLEMENT DAY`),as.Date(group$CASHFLOW$DATE),version)
  } else {
    year_diff <- newdaycount(as.Date(group$TODAY),as.Date(group$CASHFLOW$DATE),version)
  }
  
  
  # the number of rows of the matrix is the number of 
  # maturity dates of the bond with the longest maturity
  # all missing elements of the matrix are filled up with zeros 
  
  MATURITYMATRIX <-
    mapply(function(i) c(year_diff[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  if (include_price == TRUE) {MATURITYMATRIX <- rbind(rep(0,n_of_bonds),
                                                      MATURITYMATRIX)}  
  
  #colnames(MATURITYMATRIX) <- group$ISIN
  
  if (is.null(dim(MATURITYMATRIX))){
    MATURITYMATRIX <- matrix(MATURITYMATRIX, nrow = 1)
    #MATURITYMATRIX <- rbind(MATURITYMATRIX[1:length(MATURITYMATRIX)],rep(0,n_of_bonds))
  } 
  colnames(MATURITYMATRIX) <- group$ISIN
  MATURITYMATRIX
}

gi <- function(t,T,i,s){
  g <- rep(0,length(t))
  for(j in seq_along(t)){
    if(i==1){
      if(T[i]<=t[j]&t[j]<T[i+1]){
        g[j] <- (T[i])^2/6 + ((T[i])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
      }
      if(t[j]>=T[i+1]){
        g[j] <- (T[i+1])*((2*T[i+1]-T[i])/6 + (t[j]-T[i+1])/2)
      }   
    }
    if(i>1&i<length(T)){
      if(t[j]<T[i-1]){
        g[j] <- 0
      }
      if(T[i-1]<=t[j]&t[j]<T[i]){
        g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
      }
      if(T[i]<=t[j]&t[j]<T[i+1]){
        g[j] <- (T[i]-T[i-1])^2/6 + ((T[i]-T[i-1])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
      }
      if(t[j]>=T[i+1]){
        g[j] <- (T[i+1]-T[i-1])*((2*T[i+1]-T[i]-T[i-1])/6 + (t[j]-T[i+1])/2)
      }
    }
    if(i==length(T)){
      if(t[j]<T[i-1]){
        g[j] <- 0
      }
      if(T[i-1]<=t[j]&t[j]<=T[i]){
        g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
      }
    } 
    if(i==s){
      g[j] <- t[j]
    }
  }
  g
}


#########################################################
### Cubic splines estimation method for 'couponbonds' after McCulloch ###
### Forumal and Function as in source 13              ###
#########################################################

estim_cs_mc <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",freque= 2) UseMethod("estim_cs_mc")
###   for Daycountversion

### Cubic spline term structure estimation 
estim_cs_mc.couponbonds <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",freque= 2) {
  
  ## data preprocessing 
  ###   for Daycountversion
  prepro <- prepro_bond(group=group,bonddata=bonddata,
                        matrange=matrange,version=version, SettlementCheck=SettlementCheck,freque=freque, mythology=mythology, Annulized=Annulized, Actual = Actual, wmethod=wmethod)
  
  
  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
  m_p_2=prepro$m_p_2
  
  
  # Choosing knot points (McCulloch)
  # number of bonds in each group 
  K <- mapply(function(k) ncol(m[[k]]),sgroup,SIMPLIFY=FALSE)  
  # number of basis functions
  s <-  mapply(function(k) round(sqrt(K[[k]])),sgroup,SIMPLIFY=FALSE)
  
  # only perform spline estimation if number of bonds per group >= 9
  if(sum(s>=3) != length(s))  stop(cat("Estimation aborted:
    For cubic splines estimation more than 9 observations per group are required","\n",
                                       "Check group(s):", group[which((s>=3)==FALSE)]),
                                   "\n" )
  
  # only used for knot point finding
  i <- mapply(function(k) 2:(max(2,(s[[k]]-2))),sgroup,SIMPLIFY=FALSE)  
  
  h <-  mapply(function(k) trunc(((i[[k]]-1)*K[[k]])/(s[[k]]-2)),sgroup,SIMPLIFY=FALSE)
  
  theta <- mapply(function(k)((i[[k]]-1)*K[[k]])/(s[[k]]-2)-h[[k]],sgroup,SIMPLIFY=FALSE)
  
  # calculate knot points
  T <- mapply(function(k) if(s[[k]]>3) c(floor(min(y[[k]][,1])),
                                         apply(as.matrix(m[[k]][,h[[k]]]),2,max)
                                         + theta[[k]]*(apply(as.matrix(m[[k]][,h[[k]]+1]),2,max)-apply(as.matrix(m[[k]][,h[[k]]]),2,max)),
                                         max(m[[k]][,ncol(m[[k]])])) else c(floor(min(y[[k]][,1])),max(m[[k]][,ncol(m[[k]])])),sgroup,SIMPLIFY=FALSE)
  
  # parameter estimation with OLS
  # dependent variable
  Y <- mapply(function(k) apply(cf_p[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
  
  
  # k ... group index
  # sidx ... index for spline function  
  # independetn variable	
  X <- mapply(function(k) mapply(function(sidx)  apply(cf[[k]]*mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]])),2,sum), 1:s[[k]] ),sgroup,SIMPLIFY=FALSE)
  
  # OLS parameter estimation
  regout <- mapply(function(k) lm(-Y[[k]]~X[[k]]-1,weights = unname(duration[[k]][,3])),sgroup,SIMPLIFY=FALSE) # parameter vector
  
  # estimated paramters  
  alpha <- lapply(regout, coef)
  
  # calculate discount factor matrix 
  dt <- list()
  for (k in sgroup){
    dt[[k]] <- matrix(1,nrow(m[[k]]),ncol(m[[k]]))
    for(sidx in 1:s[[k]]){
      dt[[k]] <- dt[[k]] + alpha[[k]][sidx]* mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]]))
    }
  }  
  
  # calculate estimated prices 
  phat <- mapply(function(k) apply(cf[[k]]*dt[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
  
  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]],m_p_2[[k]], freque=freque, mythology=mythology, Annulized=Annulized),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"
  
  # maturity interval
  t <- mapply(function(k) seq(max(round(min(T[[k]]),2),0.01),max(T[[k]]),0.01), sgroup,SIMPLIFY=FALSE)  
  
  # calculate mean and variance of the distribution of the discount function 
  mean_d <- mapply(function(k) apply(mapply(function(sidx) alpha[[k]][sidx]*
                                              gi(t[[k]],T[[k]],sidx,s[[k]]),1:s[[k]]),1,sum) +1, sgroup, SIMPLIFY=FALSE)
  
  #  variance covariance matrix for estimated parameters
  if(rse) Sigma <- lapply(regout,vcovHAC.default) else Sigma <- lapply(regout,vcov) 
  
  var_d <- mapply(function(k) apply(mapply(function(sidx) gi(t[[k]],T[[k]],
                                                             sidx,s[[k]]),1:s[[k]]),1,function(x) t(x)%*%Sigma[[k]]%*%x), sgroup, SIMPLIFY=FALSE) 
  
  # lower 95% confidence interval
  cl <- mapply(function(k) mean_d[[k]] + rep(qt(0.025,nrow(X[[k]])- ncol(X[[k]])),
                                             length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE)	
  # upper 95 % confidence interval	
  cu <- mapply(function(k) mean_d[[k]] + rep(qt(0.975,nrow(X[[k]])- ncol(X[[k]])),
                                             length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE) 
  
  # zero cupon yield curves for maturity interval t 
  zcy_curves <-  mapply(function(k)  cbind(t[[k]],-log(mean_d[[k]])/t[[k]],-log(cl[[k]])/t[[k]],
                                           -log(cu[[k]])/t[[k]]),sgroup, SIMPLIFY=FALSE )  
  
  
  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
  
  # calculate spread curves           	    
  if(n_group != 1) {     
    srange <- seq(max(unlist(lapply(t,min))),min(unlist(lapply(t,max))),0.01)
    s_curves <- mapply(function(k) cbind(srange,zcy_curves[[k]][c(which(zcy_curves[[k]][,1]== srange[1]): which(zcy_curves[[k]][,1] == srange[length(srange)])),2] 
                                         - zcy_curves[[1]][c(which(zcy_curves[[1]][,1]== srange[1]): which(zcy_curves[[1]][,1]== srange[length(srange)])),2]),sgroup, SIMPLIFY=FALSE) 
    
  } else s_curves = "none"
  for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
  class(s_curves) <- "s_curves"
  
  # create discount factor curves 
  df_curves <- mapply(function(k) cbind(t[[k]],mean_d[[k]]),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
  class(df_curves) <- "df_curves"	
  
  # calculate forward rate curves
  
  fwr_curves <- mapply(function(k) cbind(t[[k]],impl_fwr(zcy_curves[[k]][,1],zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"	
  
  # return list of results
  result <- list(  group=group,          # e.g. countries, rating classes
                   matrange=matrange,    # maturity range of bonds
                   n_group=n_group,      # number of groups
                   knotpoints=T,         # knot points
                   rse=rse,              # robust standard errors
                   spot=zcy_curves, 	# zero coupon yield curves
                   spread=s_curves,      # spread curves
                   discount=df_curves,	# forward rate curves
                   forward=fwr_curves,	# discount factor curves
                   cf=cf,                # cashflow matrix
                   m=m,                  # maturity matrix
                   p=p,                  # dirty prices
                   phat=phat,            # estimated prices
                   perrors=perrors,     	# price errors
                   y=y,                  # maturities and yields
                   yhat=yhat,            # estimated yields
                   yerrors=yerrors,	# yield errors
                   alpha=alpha,          # cubic splines parameters                             
                   regout=regout,        # OLS output
                   s=s,                   #Needed for LOOC
                   duration=duration
  )
  
  # assign names to results list 
  for ( i in 6:length(result)) names(result[[i]]) <- group 
  
  class(result) <- "termstrc_cs"
  result
  
}

estim_cs_mc.dyncouponbonds <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",freque= 2){
  solution_list <- vector(mode = "list", length = length(bonddata))
  for (i in 1:length(bonddata)){
    cs_res <- estim_cs_mc(bonddata[[i]], group=group, matrange=matrange,rse=rse, version = version, SettlementCheck=SettlementCheck, mythology = mythology, Annulized=Annulized, Actual = Actual,wmethod=wmethod,freque=freque)
    solution_list[[i]] <- cs_res
  }
  class(solution_list) <- "dyncouponbonds_cs"
  return(solution_list)
}


#### CHECK1 ####
##########################################################################
### Nelson/Siegel-type yield curve estimation method for 'couponbonds' ###
##########################################################################

estim_nss <- function(dataset, ...) UseMethod("estim_nss")

estim_nss.couponbonds <- function(dataset,                  # dataset (static)
                                  group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
                                  matrange="all",            # maturity range in years c(min, max) 
                                  method="ns",
                                  version = 0,              #   for Daycountversion
                                  SettlementCheck = FALSE,   #   for SettlementDay Usage
                                  freque=2,                  #   for Frequency of Coupon Payment
                                  mythology=TRUE,            #   for Different Yield Method 
                                  Annulized=TRUE,            #   for Yield Method 2 (Annulized/Semiannulized)
                                  Actual = 0,                #   for Daycount convention Yield Method 2
                                  wmethod="rf",              #   for Weight Method
                                  startparam=NULL,           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
                                  # otherwise globally optimal parameters are searched automatically
                                  lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                  tauconstr = NULL,          # constraints for tau parameter grid
                                  constrOptimOptions = list(control = list(maxit = 1000), outer.iterations = 200, outer.eps = 1e-03),...
) {
  
  ## data preprocessing
  ###   for Daycountversion
  prepro <- prepro_bond(group=group,bonddata=dataset,matrange=matrange, version=version, SettlementCheck=SettlementCheck,freque=freque, mythology=mythology, Annulized=Annulized, Actual = Actual, wmethod=wmethod)
  
  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
  m_p_2 = prepro$m_p_2
  
  ## automatically determine globally optimal start parameters
  spsearch <- list()
  length(spsearch) <- n_group
  
  ## default tau constraints (if not specified by user)
  if (is.null(tauconstr)) {
    tauconstr <- list()
    length(tauconstr) <- n_group
    names(tauconstr) <- group
    for (k in sgroup){
      tauconstr[[k]] <- c(min(m[[k]][1,]), max(m[[k]]), 0.25, 0.5)
      names(tauconstr[[k]]) <- c("tau_min", "tau_max", "gridstepsize", "tau_distance")
      if (method == "asv") {tauconstr[[k]][4] = 0.5}
      if (method == "ns") {tauconstr[[k]] = tauconstr[[k]][1:3]}
    }
    if (method!="dl"){
      print("The following constraints are used for the tau parameters:")
      print(tauconstr)}
  }
  
  if(is.null(startparam)){
    startparam <- matrix(ncol = 6, nrow = n_group)
    
    colnames(startparam) <- c("beta0","beta1","beta2","tau1","beta3","tau2")
    
    if (method == "dl") {startparam <- startparam[,1:3, drop=FALSE]}
    if (method == "ns") {startparam <- startparam[,1:4, drop=FALSE]}
    
    for (k in sgroup){
      ## check if grid size is too large
      if((tauconstr[[k]][1] + tauconstr[[k]][3]) > (tauconstr[[k]][2] - tauconstr[[k]][3])) warning("Grid step size is too large!")
      if((tauconstr[[k]][1] + 3*tauconstr[[k]][3]) > tauconstr[[k]][2]) warning("Grid step size is too large!")
      
      print(paste("Searching startparameters for ", group[k]))
      spsearch[[k]] <- findstartparambonds(p[[k]],m[[k]],cf[[k]], duration[[k]][,3],
                                           method, tauconstr[[k]])
      startparam[k,] <- spsearch[[k]]$startparam 
      print(startparam[k,])
    }
  }
  
  rownames(startparam) <- group
  
  ## objective function (weighted price error minimization) 
  objfct <- get_objfct_bonds(method)
  
  ## gradient objective function
  grad_objfct <- get_grad_objfct_bonds(method)
  
  ## calculate optimal parameter vectors
  constraints <- list()
  for (k in sgroup){
    constraints[[k]] <- get_constraints(method, tauconstr[[k]])
  }
  opt_result <- list()
  for (k in sgroup){
    opt_result[[k]] <- estimatezcyieldcurve(method, startparam[k,], lambda, objfct, grad_objfct, constraints[[k]],
                                            constrOptimOptions, m[[k]], cf[[k]], duration[[k]][,3], p[[k]]) 
  }
  
  ## data post processing 
  postpro <- postpro_bond(opt_result,m,cf,sgroup,n_group,y,p,ac,m_p,m_p_2,method,lambda,freque,mythology,Annulized)
  
  ## return list of results 
  result <- list(group=group,                   # e.g. countries, rating classes
                 matrange=matrange,             # maturity range of bonds
                 method=method,                 # estimation method
                 startparam=startparam,         # calculated startparameters
                 n_group=n_group,               # number of groups,
                 lambda=lambda,                 # lambda parameter for dl
                 spsearch = spsearch,           # detailed data from start param search
                 spot=postpro$zcy_curves,       # zero coupon yield curves
                 spread=postpro$s_curves,       # spread curves
                 forward=postpro$fwr_curves,    # forward rate curves
                 discount=postpro$df_curves,    # discount factor curves
                 expoints=postpro$expoints,     # extrapolation points
                 cf=cf,                         # cashflow matrix
                 m=m,                           # maturity matrix
                 duration=duration,             # duration, modified duration, weights
                 p=p,                           # dirty prices        
                 phat=postpro$phat,             # estimated dirty prices         
                 perrors=postpro$perrors,       # price errors
                 ac=ac,                         # accrued interest
                 y=y,                           # maturities and yields
                 yhat=postpro$yhat,             # estimated yields
                 yerrors=postpro$yerrors,       # yield errors
                 opt_result=opt_result          # optimisation results           
  )
  
  for ( i in 7:length(result)) names(result[[i]]) <- group
  class(result) <- "termstrc_nss"
  result
}

### Estimate zero-coupon yield curve

estimatezcyieldcurve <- function(method, startparam, lambda, objfct, grad_objfct, constraints, constrOptimOptions, m, cf, weights, p) {
  
  if(method=="dl") {
    opt_result <- constrOptim(theta = startparam,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              lambda, m, cf, weights, p)
  } else {
    opt_result <- constrOptim(theta = startparam,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              m, cf, weights, p)
  }
  opt_result
}

### Start parameter search routine for bond data

findstartparambonds <- function(p,m,cf, weights, method, tauconstr,
                                control = list(), outer.iterations = 30, outer.eps = 1e-03) {
  
  epsConst <- 0.0001 ## ensures that starting value for constrOptim can not be on the boundary
  
  if(method=="dl"){
    startparam = rep(1,3)
    tau = NULL
    fmin = NULL
    optind = NULL
  }
  
  if(method=="ns"){
    tau <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    fmin <- rep(NA, length(tau))
    lsbeta <- matrix(nrow = length(tau), ncol = 4)
    
    ui <- rbind(c(1,0,0),                 # beta0 > 0
                c(1,1,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)### HILFE###
    
    for (i in 1:length(tau)){
      
      lsparam <- constrOptim(theta = rep(1,3), # start parameters for D/L, objective function is convex
                             f = objfct_ns_bonds_grid,
                             grad = grad_ns_bonds_grid,
                             ui = ui,
                             ci = ci,
                             mu = 1e-04,
                             control = control,
                             method = "BFGS",
                             outer.iterations = outer.iterations,
                             outer.eps = outer.eps,
                             tau[i], m, cf, weights, p) ## additional inputs for f and grad
      beta <- c(lsparam$par,tau[i])
      fmin[i] <- lsparam$value
      lsbeta[i,] <- beta 
    }
    optind <- which(fmin == min(fmin, na.rm = TRUE))
    startparam <- lsbeta[optind,]
  }
  
  if(method=="sv") {
    
    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)### HILFE###
    
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    for (i in 1:length(tau1))
    {
      #print(i) # DEBUG
      for (j in 1:length(tau2))
      {
        
        if(tau1[i] + tauconstr[4] < tau2[j]) {
          #print(j) # DEBUG
          
          lsparam <- constrOptim(theta = rep(0.01,4),
                                 f = objfct_sv_bonds_grid,
                                 grad = grad_sv_bonds_grid,
                                 ui = ui,
                                 ci = ci,
                                 mu = 1e-04,
                                 control = control,
                                 method = "BFGS",
                                 outer.iterations = outer.iterations,
                                 outer.eps = outer.eps,
                                 c(tau1[i], tau2[j]), m, cf, weights, p) ## additional inputs for f and grad
          
          beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
          
          fmin[i,j] <- lsparam$value
          lsbeta[(i-1)*length(tau1)+j,] <- beta
        }
      }
    }
    
    optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
    startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]    
  }
  
  if(method=="asv") {
    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0)### HILFE###
    
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    for (i in 1:length(tau1))
    {
      for (j in 1:length(tau2))
      {
        
        if(tau1[i] < tau2[j]) {
          
          lsparam <- constrOptim(theta = rep(0.01,4),
                                 f = objfct_asv_bonds_grid,
                                 grad = grad_asv_bonds_grid,
                                 ui = ui,
                                 ci = ci,
                                 mu = 1e-04,
                                 control = control,
                                 method = "BFGS",
                                 outer.iterations = outer.iterations,
                                 outer.eps = outer.eps,
                                 c(tau1[i], tau2[j]), m, cf, weights, p) ## additional inputs for f and grad
          
          beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
          
          fmin[i,j] <- lsparam$value
          lsbeta[(i-1)*length(tau1)+j,] <- beta
        }
      }
    }
    
    optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
    startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]   
  }
  result <- list(startparam = startparam, tau = tau, fmin = fmin, optind = optind)
  class(result) <- "spsearch"
  result
}

### Startparameter grid search plots

plot.spsearch <- function(x, main = "Start parameter search", rgl = TRUE, ...) {
  
  if(is.matrix(x$tau)){
    image(x$tau[,1],x$tau[,2],log(x$fmin),xlab = expression(tau[1]), ylab = expression(tau[2]),main = "Log(Objective function)")
    contour(x$tau[,1],x$tau[,2],log(x$fmin),nlevels=10,xlab = expression(tau[1]), ylab = expression(tau[2]),main = "Log(Objective function)", add = TRUE)
    points(x$tau[x$optind[1],1],x$tau[x$optind[2],2],pch = 10, col = "steelblue")
    if (rgl) {
      open3d()
      persp3d(x$tau[,1], x$tau[,2], log(x$fmin), col = "green3", box = FALSE,xlab = "tau_1", ylab = "tau_2", zlab = "Log(Objective function)")
      points3d(x$tau[x$optind[1],1],x$tau[x$optind[2],2],min(log(x$fmin), na.rm = TRUE), col = "red")
    }
    else {
      par(ask = TRUE)
      persp(x$tau[,1], x$tau[,2], log(x$fmin), col = "green3", box = TRUE, xlab = "tau_1", ylab = "tau_2", zlab = "Log(Objective function)",
            shade = TRUE, ticktype = "detailed", border = NA, cex.lab = 1, cex.axis = 0.7,  theta = 0, phi = 25, r = sqrt(3),
            d = 1, scale = TRUE, expand = 1, ltheta = 135, lphi = 0)
    }
  } else {
    plot(x$tau,log(x$fmin),xlab = expression(tau[1]), ylab = "Log(Objective function)", type = "l", main = main)
    points(x$tau[x$optind],log(x$fmin[x$optind]),pch = 10, col = "red")
  }
}


#############################################################################
### Nelson/Siegel-type yield curve estimation method for 'dyncouponbonds' ###
#############################################################################
###   for Daycountversion
estim_nss.dyncouponbonds <- function(dataset, group, matrange="all",method="ns",version = 0,SettlementCheck=FALSE, freque=2, mythology=TRUE,Annulized=TRUE,Actual=0,wmethod="rf",
                                     lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                     tauconstr = NULL,              # tau constraints
                                     optimtype = "firstglobal", # 'firstglobal' of 'allglobal'
                                     constrOptimOptions = list(control = list(maxit = 1000), outer.iterations = 200, outer.eps = 1e-03), ...
) {
  
  res <- list()
  
  ## perform sequence of term structure estimations
  for (i in seq(length(dataset))) {
    if(i>1 && optimtype == "firstglobal"){
      ## use optimal parameters from previous period as start parameters
      b <- t(mapply(function(j) res[[i-1]]$opt_result[[j]]$par,  seq_along(group)))
      rownames(b) <- group                               
    } else b <- NULL
    
    ## static estimation
    bonddata <- dataset[[i]]
    ###   for Daycountversion
    res[[i]] <- estim_nss(bonddata, group, matrange, method=method,version=version,SettlementCheck=SettlementCheck,freque=freque, mythology=mythology,Annulized=Annulized,Actual=Actual,wmethod=wmethod, startparam=b, lambda=lambda,tauconstr,constrOptimOptions)
  }
  class(res) <- "dyntermstrc_nss"
  
  res
}



#########################################################################
### Nelson/Siegel-type yield curve estimation method for 'zeroyields' ###
#########################################################################
###   for Daycountversion
estim_nss.zeroyields <- function (dataset,
                                  method = "ns",
                                  version = 0,
                                  SettlementCheck = FALSE, 
                                  freque=2,
                                  mythology=TRUE,
                                  Annulized=TRUE,
                                  Actual = 0,
                                  lambda = 0.0609*12,
                                  tauconstr = NULL,
                                  optimtype = "firstglobal",
                                  constrOptimOptions = list(control = list(), outer.iterations = 200, outer.eps = 1e-03),...)
{
  obj <- dataset
  optparam <- matrix(nrow = nrow(obj$yields), ncol = length(get_paramnames(method))) 
  opt_result <- list()
  spsearch <- list()
  
  ## default tau constraints (if not specified by user)
  if (is.null(tauconstr) && method != "dl"){
    tauconstr <- c(min(obj$maturities), max(obj$maturities), 0.25, 0.5)
    if (method == "asv") {tauconstr[4] = 0.5}
    if (method == "ns") {tauconstr = tauconstr[1:3]}
    print("The following constraints are used for the tau parameters:")
    print(tauconstr)
  }
  
  objfct <- get_objfct(method)
  grad_objfct <- get_grad_objfct(method)
  
  sp_search <- findstartparamyields(obj$yields[1,],obj$maturities, method, tauconstr)
  startparam <- sp_search$startparam
  constraints <- get_constraints(method, tauconstr)
  
  ## Estimation loop
  
  for (i in 1:nrow(obj$yields)){  
    yields <- obj$yields[i,]
    
    if (i==1) {
      beta <- startparam
      spsearch[[i]] <- sp_search
    }
    
    if(i>1 && optimtype == "allglobal"){
      sp_search <- findstartparamyields(yields,obj$maturities, method, tauconstr)
      beta <- sp_search$startparam
      spsearch[[i]] <- sp_search
    }
    if(i>1 && optimtype == "firstglobal"){
      beta <- optparam[i-1,]
    }
    
    opt_result[[i]] <- estimateyieldcurve(method, yields, obj$maturities, beta, lambda, objfct,
                                          grad_objfct, constraints, constrOptimOptions)
    optparam[i,] <- opt_result[[i]]$par
  }
  
  colnames(optparam) <- get_paramnames(method)
  
  
  yhat <- t(apply(optparam,1, function(x) spotrates(method,x,obj$maturities,lambda)))
  
  result <- list(optparam = optparam, opt_result = opt_result, method = method,
                 maturities = obj$maturities, dates = obj$dates, spsearch = spsearch,
                 yields = obj$yields,yhat = yhat, lambda = lambda)
  class(result) <- "dyntermstrc_yields"
  result
}

estimateyieldcurve <- function(method, y, m, beta, lambda, objfct, grad_objfct, constraints, constrOptimOptions) {
  if(method=="dl") {
    opt_result <- constrOptim(theta = beta,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              lambda, m, y)
  } else {
    opt_result <- constrOptim(theta = beta,
                              f = objfct,
                              grad = grad_objfct,
                              ui = constraints$ui,
                              ci = constraints$ci,
                              mu = 1e-04,
                              control = constrOptimOptions$control,
                              method = "BFGS",
                              outer.iterations = constrOptimOptions$outer.iterations,
                              outer.eps = constrOptimOptions$outer.eps,
                              m,y)
  }     
}

findstartparamyields <- function(y,m, method, tauconstr, control = list(), outer.iterations = 30, outer.eps = 1e-03)
{
  epsConst <- 0.0001 ## ensures that starting value for constrOptim can not be on the boundary
  
  if(method=="dl"){
    startparam = rep(1,3)
    tau = NULL
    fmin = NULL
    optind = NULL
  }
  
  if(method=="ns"){
    tau <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    fmin <- rep(NA, length(tau))
    lsbeta <- matrix(nrow = length(tau), ncol = 4)
    
    ui <- rbind(c(1,0,0),                 # beta0 > 0
                c(1,1,0))                 # beta0 + beta1 > 0
    ci <- c(0,0) ### HILFE###
    
    for (i in 1:length(tau)){
      X <- cbind(rep(1,length(y)),
                 ((1 - exp(-m/tau[i]))/(m/tau[i])),
                 (((1 - exp(-m/tau[i]))/(m/tau[i])) - exp(-m/tau[i])))
      
      lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
      beta <- c(lsparam[1:3],tau[i])
      
      ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0)
      if(beta[1]>0 && ((beta[1]+beta[2])>0)){
        fmin[i] <- objfct_ns(beta, m, y)
      } else {
        ## switch to constrOptim if OLS violates constraints
        lsparam <- constrOptim(theta = rep(1,3), # start parameters for D/L, objective function is convex
                               f = objfct_ns_grid,
                               grad = grad_ns_grid,
                               ui = ui,
                               ci = ci,
                               mu = 1e-04,
                               control = control,
                               method = "BFGS",
                               outer.iterations = outer.iterations,
                               outer.eps = outer.eps,
                               tau[i], m, y) ## additional inputs for f and grad
        
        beta <- c(lsparam$par,tau[i])
        fmin[i] <- objfct_ns(beta, m, y)
      }
      lsbeta[i,] <- beta 
    }
    optind <- which(fmin == min(fmin))
    startparam <- lsbeta[optind,]
  }
  
  
  if(method=="sv"){
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    ## fminsolver <- matrix(nrow = length(tau1), ncol = length(tau2)) ## DEBUG
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    
    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0) ### HILFE###
    
    for (i in 1:length(tau1))
    {
      for (j in 1:length(tau2))
      {
        ## minimum tau distance constraint
        if(tau1[i] + tauconstr[4] < tau2[j]) {
          
          ## reparametrize to avoid nonsingular matrix
          if (i == j){
            X <- cbind(rep(1,length(y)),
                       ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                       - exp(-m/tau1[i]))
            
            lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
            beta <- c(lsparam[1],lsparam[2]-lsparam[3],lsparam[3]/2,tau1[i],lsparam[3]/2,tau2[j])
          } else
          {
            X <- cbind(rep(1,length(y)),
                       ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                       (((1 - exp(-m/tau1[i]))/(m/tau1[i])) - exp(-m/tau1[i])),
                       (((1 - exp(-m/tau2[j]))/(m/tau2[j])) - exp(-m/tau2[j])))
            
            lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
            beta <- c(lsparam[1:3],tau1[i],lsparam[4],tau2[j])
          }
          ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0, tau distance)
          ##if(beta[1]>0 && ((beta[1]+beta[2])>0  ) && (-tau1[i] + tau2[j]) > tauconstr[4]){
          if(beta[1]>0 && ((beta[1]+beta[2])>0)){
            fmin[i,j] <- objfct_sv(beta, m, y)
          } else {
            ## switch to constrOptim if OLS violates constraints
            lsparam <- constrOptim(theta = rep(0.01,4),
                                   f = objfct_sv_grid,
                                   grad = grad_sv_grid,
                                   ui = ui,
                                   ci = ci,
                                   mu = 1e-04,
                                   control = control,
                                   method = "BFGS",
                                   outer.iterations = outer.iterations,
                                   outer.eps = outer.eps,
                                   c(tau1[i], tau2[j]), m, y) ## additional inputs for f and grad
            
            beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
            fmin[i,j] <- lsparam$value
          }
          lsbeta[(i-1)*length(tau1)+j,] <- beta
        }
      } 
    }
    optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
    startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]
  }
  
  if(method=="asv"){
    tau1 <- seq(tauconstr[1] + epsConst, tauconstr[2] - epsConst, tauconstr[3])
    tau2 <- tau1
    tau <- cbind(tau1, tau2)
    fmin <- matrix(nrow = length(tau1), ncol = length(tau2))
    lsbeta <- matrix(nrow = length(tau1)*length(tau2), ncol = 6)
    
    ui <- rbind(c(1,0,0,0),                 # beta0 > 0
                c(1,1,0,0))                 # beta0 + beta1 > 0
    ci <- c(0,0) ### HILFE###
    
    for (i in 1:length(tau1))
    {
      for (j in 1:length(tau2))
      {
        if(tau1[i] < tau2[j]) {
          
          X <- cbind(rep(1,length(y)),
                     ((1 - exp(-m/tau1[i]))/(m/tau1[i])),
                     (((1 - exp(-m/tau1[i]))/(m/tau1[i])) - exp(-m/tau1[i])),
                     (((1 - exp(-m/tau2[j]))/(m/tau2[j])) - exp(-2*m/tau2[j])))
          
          lsparam <- solve(t(X)%*%X)%*%t(X)%*%y
          beta <- c(lsparam[1:3],tau1[i],lsparam[4],tau2[j])
          
          ## check parameter contraints (beta_0 > 0, beta_0 + beta_1 > 0, tau distance)
          ##if(beta[1]>0 && ((beta[1]+beta[2])>0) && (-tau1[i] + tau2[j]) > 0){
          if(beta[1]>0 && ((beta[1]+beta[2])>0)){
            fmin[i,j] <- objfct_sv(beta, m, y)
          } else {
            ##print(paste("OLS violated constraints, solving with constrOptim, tau_1 =",tau1[i],"and tau_2 = ", tau2[j])) ## DEBUG
            ## switch to constrOptim if OLS violates constraints
            lsparam <- constrOptim(theta = rep(0.01,4),
                                   f = objfct_asv_grid,
                                   grad = grad_asv_grid,
                                   ui = ui,
                                   ci = ci,
                                   mu = 1e-04,
                                   control = control,
                                   method = "BFGS",
                                   outer.iterations = outer.iterations,
                                   outer.eps = outer.eps,
                                   c(tau1[i], tau2[j]), m, y) ## additional inputs for f and grad
            
            beta <- c(lsparam$par[1:3],tau1[i],lsparam$par[4],tau2[j])
            fmin[i,j] <- lsparam$value
          }
          lsbeta[(i-1)*length(tau1)+j,] <- beta
        }
      } 
    }
    optind <- which(fmin == min(fmin, na.rm = TRUE),arr.ind=TRUE)
    startparam <- lsbeta[(optind[1]-1)*length(tau1) + optind[2],]
  }
  
  result <- list(startparam = startparam, tau = tau, fmin = fmin, optind = optind)
  class(result) <- "spsearch"
  result
}

## Root mean squared error
rmse <- function (actual,estimated) {
  e <- actual - estimated
  sqrt(mean(e^2))
}

rmsewe <-function (actual,estimated,weights){
  e <- actual - estimated	
  sqrt(mean(weights*(e^2)))
} 

##  Mean absolute error
aabse <-function (actual,estimated){
  e <- actual - estimated	
  mean(abs(e))
} 

## Symmetric mean absolute percentage error (Nach Wikipedia)
smape <- function(actual, estimated){
  if (length(which(actual==0))!=0){
    actual <- actual[-which(actual==0)]
  }
  if (length(which(estimated==0))!=0){
    estimated<- estimated[-which(estimated==0)]
  }
  
  e <- (abs(estimated - actual))/(abs(actual)+abs(estimated))
  mean(e)
}

smape1 <- function(actual, estimated){
  e <- (abs(estimated - actual))/(abs(actual)+abs(estimated))
  e
}

### Gradient of Svensson grid loss function for bonds
grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("grad_sv_bonds_gridCpp", beta, tau, m, cf, w, p)
}

### Gradient of Svensson loss function for bonds
grad_sv_bonds <- function(beta,m,cf,w,p){
  
  emt1 <- exp(-m/beta[4])
  emt2 <- exp(-m/beta[6])
  oemt1 <- (1-emt1)
  t1emt1 <- beta[4]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-emt2 + beta[6]*(1 - emt2)/m)
  emt1t <- emt1/beta[4]
  emt2t <- emt2/beta[6]
  
  a <- exp((-beta[1] - beta[3]*emt1tm - beta[5]*emt2tm - (beta[2]*t1emt1)/m)*m/100)
  
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  test3 <- dm*emt2tm
  test3[is.na(test3)] <- 0 #new
  test4 <- dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4]))
  test4[is.na(test4)] <- 0 #new
  test5 <-dm*beta[5]*( -emt2t + (1-emt2)/m - emt2t*m/beta[6])
  test5[is.na(test5)] <- 0 #ne
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2,parallel=TRUE, na.rm=FALSE)))
  gbeta5 <- sum(b*(-colsums(test3,parallel=TRUE, na.rm=FALSE)))     
  gbeta4 <- sum(b*(-colsums(test4,parallel=TRUE,na.rm=FALSE)))
  gbeta6 <- sum(b*(-colsums(test5,parallel=TRUE,na.rm=FALSE)))
  
  c(gbeta1,gbeta2,gbeta3,gbeta4,gbeta5,gbeta6)
  
}

### Gradient of Nelson/Siegel grid loss function for bonds
#check
grad_ns_bonds_grid <- function(beta, tau, m, cf, w, p){
  emt1 <- exp(-m/tau)
  oemt1 <- (1-emt1)
  t1emt1 <- tau*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt1t <- emt1/tau
  
  a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
                (beta[2]*t1emt1)/m)*m)/100)
  
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2, na.rm=FALSE)))    
  
  c(gbeta1,gbeta2,gbeta3)
}

### Gradient of Nelson/Siegel grid loss function for bonds
grad_ns_bonds <- function(beta, m, cf, w, p){
  emt1 <- exp(-m/beta[4])
  oemt1 <- (1-emt1)
  t1emt1 <- beta[4]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt1t <- emt1/beta[4]
  
  a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
                (beta[2]*t1emt1)/m)*m)/100)
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  #new
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  test4 <- dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4]))
  test4[is.na(test4)] <- 0 #new
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2,parallel=TRUE, na.rm=FALSE)))    
  gbeta4 <- sum(b*(-colsums(test4,parallel=TRUE,na.rm=FALSE)))
  
  c(gbeta1,gbeta2,gbeta3,gbeta4)
}

### Gradient of Adjusted Svensson grid loss function for bonds
grad_asv_bonds_grid <- function(beta, tau, m, cf, w, p){
  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  oemt1 <- (1-emt1)
  t1emt1 <- tau[1]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/tau[2]) + (tau[2]*(1 -emt2))/m)
  emt1t <- emt1/tau[1]
  
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
               beta[4]*emt2tm + 
               (beta[2]*t1emt1)/m)*m/100)
  
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  #new
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  test3 <- dm*emt2tm
  test3[is.na(test3)] <- 0 #new
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2,parallel=TRUE, na.rm=FALSE)))
  gbeta5 <- sum(b*(-colsums(test3,parallel=TRUE, na.rm=FALSE)))     
  
  c(gbeta1, gbeta2, gbeta3, gbeta5)
  
}

### Gradient of Adjusted Svensson loss function for bonds
grad_asv_bonds <-  function(beta, m, cf, w, p){
  emt1 <- exp(-m/beta[4])
  emt2 <- exp(-m/beta[6])
  oemt1 <- (1-emt1)
  t1emt1 <- beta[4]*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-exp(-2*m/beta[6]) + (beta[6]*(1 -emt2))/m)
  emt1t <- emt1/beta[4]
  
  a <- exp(-(beta[1] + beta[3]*emt1tm + 
               beta[5]*emt2tm + 
               (beta[2]*t1emt1)/m)*m/100)
  
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  #new
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  test3 <- dm*emt2tm
  test3[is.na(test3)] <- 0 #new
  test4 <- dm*(-beta[2]*emt1t + beta[2]*oemt1/m + beta[3]*(-emt1t+ oemt1/m - emt1t*m/beta[4]))
  test4[is.na(test4)] <- 0 #new
  test5 <-dm*beta[5]*(-emt2/beta[6] + (1-emt2)/m - 2*exp(-2*m/beta[6])*m/beta[6]^2 )
  test5[is.na(test5)] <- 0 #new
  
  # beta[5]*(-exp(-m/beta[6])/beta[6] + (1-exp(-m/beta[6]))/m - 2*exp(-2*m/beta[6])*m/beta[6]^2 
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2, parallel=TRUE,na.rm=FALSE)))
  gbeta4 <- sum(b*(-colsums(test4 ,parallel=TRUE,na.rm=FALSE)))
  gbeta5 <- sum(b*(-colsums(test3,parallel=TRUE, na.rm=FALSE)))     
  gbeta6 <- sum(b*(-colsums(test5,parallel=TRUE,na.rm=FALSE)))
  
  c(gbeta1, gbeta2, gbeta3, gbeta4, gbeta5, gbeta6)
  
}

### Gradient of Diebold/Li loss function for bonds
grad_dl_bonds <- function(beta, lambda, m, cf, w, p){
  tau <- 1/lambda 
  emt1 <- exp(-m/tau)
  oemt1 <- (1-emt1)
  t1emt1 <- tau*oemt1
  emt1tm <- (-emt1 + t1emt1/m)
  emt1t <- emt1/tau
  
  a <- exp(-((beta[1] + beta[3]*(-emt1 +t1emt1/m) + 
                (beta[2]*t1emt1)/m)*m)/100)
  
  acf <- a*cf
  acf1 <- a*cf
  acf1[is.na(acf1)] <- 0 #new
  b <- -2*w*(p-colsums(acf1,parallel=TRUE,na.rm=FALSE))
  d <- acf/100
  dm <- d*m
  
  dm1 <- d*m
  dm1[is.na(dm1)] <- 0 #new
  test1 <- d*t1emt1
  test1[is.na(test1)] <- 0 #new
  test2 <-dm*emt1tm
  test2[is.na(test2)] <- 0 #new
  
  gbeta1 <- sum(b*(-colsums(dm1,parallel=TRUE,na.rm=FALSE)))
  gbeta2 <- sum(b*(-colsums(test1,parallel=TRUE, na.rm=FALSE)))
  gbeta3 <- sum(b*(-colsums(test2,parallel=TRUE, na.rm=FALSE)))    
  
  c(gbeta1,gbeta2,gbeta3)
}

### Faster version of column sums
cSums <- function (x, na.rm = FALSE, dims = 1L) {
  dn <- dim(x)
  n <- prod(dn[1L:dims])
  dn <- dn[-(1L:dims)]
  z <-  .Internal(colSums(x, n, prod(dn), na.rm))
  z
}

### Gradient of Nelson/Siegel loss function for yields
grad_ns <- function(beta, m, y)
{
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
              (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
                                             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),
    
    sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
             beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))
        *(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
            (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
  )
}


grad_ns_grid <- function(beta, tau, m, y)
{
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
              (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
    
    sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
                                           (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),
    
    sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
             (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
  )
}



### Gradient of Svensson loss function for yields
grad_sv <- function(beta, m, y)
{
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
              beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
              (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
                                             beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
                                             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),
    
    sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
             beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(-2*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - m/(beta[6]^2*exp(m/beta[6])))*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp(-m/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
  )
}


grad_sv_grid <- function(beta, tau, m, y)
{
  .Call("grad_sv_gridCpp", beta, tau, m, y)
  
  ##         c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
  ##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
  ##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
  
  ##           sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
  ##         beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
  ##         (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),
  
  ##           sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
  ##     (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
  ##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
  ##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
  
  
  ##           sum(-2*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)*
  ##     (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
  ##       beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
  ##       (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
  ##           )
}


### Gradient of adjusted Svensson loss function for yields
grad_asv <- function(beta, m, y)
{
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
              beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
              (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum((-2*beta[4]*(1 - exp(-m/beta[4]))*(-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
                                             beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
                                             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))/m),
    
    sum(-2*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(2*(beta[2]/(beta[4]*exp(m/beta[4])) - (beta[2]*(1 - exp(-m/beta[4])))/m - 
             beta[3]*(-(1/(beta[4]*exp(m/beta[4]))) + (1 - exp(-m/beta[4]))/m - m/(beta[4]^2*exp(m/beta[4]))))*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(-2*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y)),
    
    sum(-2*beta[5]*(-(1/(beta[6]*exp(m/beta[6]))) + (1 - exp(-m/beta[6]))/m - 
                      (2*m)/(beta[6]^2*exp((2*m)/beta[6])))*
          (-beta[1] - beta[3]*(-exp(-m/beta[4]) + (beta[4]*(1 - exp(-m/beta[4])))/m) - 
             beta[5]*(-exp((-2*m)/beta[6]) + (beta[6]*(1 - exp(-m/beta[6])))/m) - 
             (beta[2]*beta[4]*(1 - exp(-m/beta[4])))/m + y))
  )
}



grad_asv_grid <- function(beta, tau, m, y)
{
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
              beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
              (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
    
    sum((-2*tau[1]*(1 - exp(-m/tau[1]))*(-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
                                           beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
                                           (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))/m),
    
    sum(-2*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
             beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
             (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y)),
    
    
    sum(-2*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m)*
          (-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - 
             beta[4]*(-exp((-2*m)/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - 
             (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m + y))
  )
}



grad_dl <- function(beta,lambda, m, y)
{
  tau <- 1/lambda
  c(sum(-2*(-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
              (beta[2]*tau*(1 - exp(-m/tau)))/m + y)),
    
    sum((-2*tau*(1 - exp(-m/tau))*(-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
                                     (beta[2]*tau*(1 - exp(-m/tau)))/m + y))/m),
    
    sum(-2*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m)*
          (-beta[1] - beta[3]*(-exp(-m/tau) + (tau*(1 - exp(-m/tau)))/m) - 
             (beta[2]*tau*(1 - exp(-m/tau)))/m + y))
  )
}


print.dyncouponbonds <- function(x, ...) {
  cat("This is a dynamic data set of coupon bonds.\n")
  cat(paste("There are",length(x), "observations between",x[[1]][[1]]$TODAY, "and",x[[length(x)]][[1]]$TODAY,".\n"))
}

print.couponbonds <- function(x, ...) {
  cat("This is a data set of coupon bonds for:\n")
  cat(names(x),",","\n")
  cat(paste("observed at ", x[[1]]$TODAY,".","\n",sep=""))
  
}



plot.ir_curve <- function(x,ylim=c(),xlim=c(),lwd=2, type="l",
                          xlab="Maturity (years)",ylab="Zero-coupon yields (in percent)", 
                          col="steelblue",lty=1, ...) 
{
  plot(x[,1] ,x[,2]*100, type=type, ylim=ylim, xlim=xlim, xlab=xlab,
       ylab=ylab,lwd=lwd,lty=lty,col=col, ... )
  
}

plot.spot_curves <- function(x,multiple= FALSE,
                             ylim= c(range(mapply(function(i) 
                               range(x[[i]][,2]),seq(x))))*100,xlim=c(),
                             type="l", lty=1, lwd=2, expoints=NULL, 
                             ylab= "Zero-coupon yields (percent)",
                             xlab= "Maturity (years)",main="Zero-coupon yield curve",
                             ...) {
  
  if(multiple) 
  { plot(x[[which.max(mapply(function(i) max(x[[i]][,1]),
                             seq(x)))]][,1], x[[which.max(mapply(function(i) 
                               max(x[[i]][,1]), seq(x)))]][,2]*100, 
         type=type,col=which.max(mapply(function(i) max(x[[i]][,1]),seq(x))),
         lty=lty,lwd=lwd,xlab=xlab,ylab=ylab,ylim=ylim, ... )
    
    for(k in c((seq(x))[-which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]))
    { lines(x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
                    else seq(nrow(x[[k]]))),1],x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
                                                       else seq(nrow(x[[k]]))),2]*100,col=k,lwd=lwd,lty=lty, ... )
      
      if(is.numeric(expoints))
      {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
             x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*100,col=k,lwd=lwd,lty=5, ... )
      }
      title(main)
      legend("bottom",legend=names(x),col=seq(x),lty=lty,lwd=lwd)
    }
  }
  else
  {   old.par <- par(no.readonly = TRUE)
  par(ask=TRUE)
  for(k in seq(x)) 
  { 
    plot.ir_curve(x[[k]],lty=lty,lwd=lwd,xlab=xlab,ylab=ylab,ylim=ylim)
    title(main)
    legend("bottom",legend=main,col=c("steelblue"), lty = 1 , pch=c(-1))
  }	
  on.exit(par(old.par))	
  }
  
  
}


plot.fwr_curves <- function(x,multiple= FALSE,
                            ylim= c(range(mapply(function(i) range(x[[i]][,2]),
                                                 seq(x))))*100,xlim=c(),type="l", lty=1, 
                            lwd=2, expoints=NULL, ylab= "Forward rate (percent)",
                            xlab= "Maturity (years)",main="Forward rate curve",...) 
  
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,
                   multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )
  
}


plot.s_curves <- function(x,xlim=c(range(mapply(function(i) 
  range(x[[i]][,1]),seq(x)))),
  ylim=c(range(mapply(function(i) range(x[[i]][,2]),
                      seq(x))))*10000,expoints=NULL, xlab="Maturity (years)", 
  ylab="Spread (basis points)", lwd=2,lty=1, main="Spread curve", ...)
  
{  if(!is.character(x))
{
  
  plot(0,0, type="n",xlab=xlab,ylab=ylab,xlim=xlim, ylim=ylim,...)
  
  for(k in c(2:length(x)))
  { lines(x[[k]][(if(is.numeric(expoints))
    seq(expoints[k]) else seq(nrow(x[[k]]))),1],
    x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
            else seq(nrow(x[[k]]))),2]*10000,col=k,lwd=lwd,lty=lty, ... )
    
    if(is.numeric(expoints))
    {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
           x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*10000,
           col=k,lwd=lwd,lty=5, ... )
      
    }
  } 
  title(main)
  legend("topleft",legend=names(x)[-1],col=seq(x)[-1],lty=1,lwd=lwd)
} else warning("No spread curves available")
} 


plot.error <- function(x,type="b",main="", mar= c(7,6,6,2) + 0.1, oma=c(4,2,2,2) +0.1,
                       ylab="Error", ...) {
  old.par <- par(no.readonly = TRUE)
  par(mar=mar, oma=oma, ... )
  
  plot(x[,1],x[,2],axes=FALSE,pch=19,lwd=c(1,2),xlab="", ylab=ylab,type=type, ...)
  axis(1,x[,1],rownames(x),las=3,...)
  axis(2,...)
  axis(3,x[,1],round(x[,1],2),...)
  mtext("Maturity (years)",3,line=2.5)
  lines(x[,1],rep(0,nrow(x)),lty=2,lwd=1,... )
  title(xlab="ISIN", outer=TRUE,main=main,...) 
  
  on.exit(par(old.par))
}   


plot.df_curves <- function(x,multiple= FALSE,
                           ylim= c(range(mapply(function(i) range(x[[i]][,2]),
                                                seq(x))))*100,xlim=c(),type="l", lty=1,
                           lwd=2, expoints=NULL, ylab="Discount factor (percent)",
                           xlab= "Maturity (years)",main="Discount factor curve",...) 
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main, multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )
  
}

summary.dyntermstrc_nss <- function(object, ...) {
  x <- object
  
  # extract convergence info
  sumry <- list()
  sumry$convergence <- t(mapply(function(i) summary(x[[i]])$convergencegroup, seq_along(x)))
  if (length(x[[1]]$group)>1){
    colnames(sumry$convergence) <- x[[1]]$group 
  } else {
    rownames(sumry$convergence) <- x[[1]]$group 
  }
  sumry$solvermsg <- t(mapply(function(i) summary(x[[i]])$convergence, seq_along(x)))
  if (length(x[[1]]$group)>1){
    colnames(sumry$solvermsg) <- x[[1]]$group 
  } else {
    rownames(sumry$solvermsg) <- x[[1]]$group 
  }
  
  # ## Split for different group lengths 
  # buckets <- mapply(function(l) mapply(function(j) dim(x[[j]][["y"]][[l]])[1],1:length(x)),1:length(x[[1]]$group),SIMPLIFY = FALSE)
  # 
  # uniquebuckets <- mapply(function(j) unique(buckets[[j]]),1:length(buckets),SIMPLIFY = FALSE)
  
  f_p_mrsme <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_maabse <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_mrsme <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_maabse <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_smape <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_smape <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_wrmse <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_wrmse <- vector(mode = "list", length = length(x[[1]]$group))
  
  for (z in 1:length(x)){
    subset <- x[[z]]
    RMSE_p <- mapply(function(i) rmse(subset$p[[i]],subset$phat[[i]]),seq(subset$n_group))
    AABSE_p <- mapply(function(i) aabse(subset$p[[i]],subset$phat[[i]]) ,seq(subset$n_group))
    RMSE_y <- mapply(function(i) rmse(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100) ,seq(subset$n_group))
    AABSE_y <- mapply(function(i) aabse(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100) ,seq(subset$n_group))
    
    sMAPE_p <- mapply(function(i) smape(subset$p[[i]],subset$phat[[i]]),seq(subset$n_group))
    sMAPE_y <- mapply(function(i) smape(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100),seq(subset$n_group))
    
    wRMSE_p <- mapply(function(i) rmsewe(subset$p[[i]],subset$phat[[i]],subset[["duration"]][[i]][,3]),seq(subset$n_group))
    wRMSE_y <- mapply(function(i) rmsewe(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100,subset[["duration"]][[i]][,3]),seq(subset$n_group))
    
    
    for (u in 1:length(x[[1]]$group)){
      f_p_mrsme[[u]] <- c(f_p_mrsme[[u]],RMSE_p[[u]])
      f_p_maabse[[u]] <- c(f_p_maabse[[u]],AABSE_p[[u]])
      f_y_mrsme[[u]] <- c(f_y_mrsme[[u]],RMSE_y[[u]])
      f_y_maabse[[u]] <- c(f_y_maabse[[u]],AABSE_y[[u]])
      f_p_smape[[u]] <- c(f_p_smape[[u]],sMAPE_p[[u]])
      f_y_smape[[u]]<- c(f_y_smape[[u]],sMAPE_y[[u]])
      f_p_wrmse[[u]] <- c(f_p_wrmse[[u]],wRMSE_p[[u]])
      f_y_wrmse[[u]] <- c(f_y_wrmse[[u]],wRMSE_y[[u]])
      
    }
  }
  
  p_mrsme <- mapply(function(k) mean(f_p_mrsme[[k]]), 1:length(x[[1]]$group))
  p_maabse <- mapply(function(k) mean(f_p_maabse[[k]]), 1:length(x[[1]]$group))
  p_smape <- mapply(function(k) mean(f_y_mrsme[[k]]), 1:length(x[[1]]$group))
  y_mrsme <- mapply(function(k) mean(f_y_maabse[[k]]), 1:length(x[[1]]$group))
  y_maabse <- mapply(function(k) mean(f_p_smape[[k]]), 1:length(x[[1]]$group))
  y_smape <- mapply(function(k) mean(f_y_smape[[k]]), 1:length(x[[1]]$group))
  p_wRMSE <- mapply(function(k) mean(f_p_wrmse[[k]]), 1:length(x[[1]]$group))
  y_wRMSE <- mapply(function(k) mean(f_y_wrmse[[k]]), 1:length(x[[1]]$group))
  
  sumry$gof <- rbind(p_mrsme,p_maabse,p_smape,p_wRMSE,y_mrsme,y_maabse,y_smape,y_wRMSE)
  colnames(sumry$gof) <- x[[1]]$group
  rownames(sumry$gof) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices","Weighted RMS-Prices", "RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMS-Yields (in %)")
  
  class(sumry) <- "summary.dyntermstrc_nss"
  sumry
}

print.summary.dyntermstrc_nss <- function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Goodness of fit:\n")
  cat("---------------------------------------------------\n")
  
  print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
  
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("Convergence information from optim ():\n")
  cat("---------------------------------------------------\n")
  
  print.default(x$convergence)
  
}

plot.dyntermstrc_nss <- function(x,range=c(0,20), ...) {
  
  # 3D plot of zero-coupon yield curves
  tsparam <- param.dyntermstrc_nss(x)
  X <- seq(if(range[1]==0) range[1]+0.1 else range[1],range[2],0.1)
  Y <- seq(nrow(tsparam[[1]]))
  
  for(j in seq(x[[1]]$n_group)){
    Z <-  mapply(function(i) spotrates(method=x[[1]]$method,tsparam[[j]][i,],X,x[[1]]$lambda), seq(nrow(tsparam[[j]])))
    open3d()
    persp3d(X,Y,Z,col = "green3",xlab="Maturity (years)", zlab="Zero-coupon yields (in %)",ylab="Time",box=FALSE)
  }
  
}

print.dyntermstrc_nss <- function(x,...){
  cat("---------------------------------------------------\n")
  cat("Estimated",get_realnames(x[[1]]$method), "parameters:")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  cat("Number of oberservations:",length(x),"\n")
  cat("\n")
  tsparam <- param.dyntermstrc_nss(x)
  print(lapply(tsparam,summary.default))
  cat("\n")
}

param <- function(x,...) UseMethod("param") 

param.dyntermstrc_nss <- function(x,...) {
  param <- list()
  for(i in seq(x[[1]]$n_group)) param[[i]] =  t(mapply(function(j) x[[j]]$opt_result[[i]]$par,seq_along(x)))
  names(param) <- x[[1]]$group                          
  class(param) <- "dyntermstrc_param"
  param
}

param.dyntermstrc_yields <- function(x,...){
  param <- list() 
  param[[1]] <- x$optparam
  class(param) <- "dyntermstrc_param"
  param
}

summary.dyntermstrc_param <- function(object,type="none",lags=1,selectlags="Fixed", ...) {
  x <- object
  sumry <- list()
  length(sumry) <- length(x) 
  for(i in seq_along(x)) {
    
    # Augmented Dickey Fuller Test for levels
    sumry[[i]]$adflevels <- apply(x[[i]],2,function(x) ur.df(x,type=type,lags=lags,selectlags=selectlags)) #alternatively use adf.testx  
    
    sumry[[i]]$adflevelsm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x[[i]]))
    
    for (j in 1:length(sumry[[i]]$adflevels)) {
      sumry[[i]]$adflevelsm[switch(type,"none"=1,"trend"=c(1:3)),j] <- sumry[[i]]$adflevels[[j]]@teststat # adf.test : $statisic
      sumry[[i]]$adflevelsm[switch(type,"none"=2,"trend"=4),j] <- sumry[[i]]$adflevels[[j]]@lags # adf.test: $parameter
      sumry[[i]]$adflevelsm[switch(type,"none"=3,"trend"=c(5:7)),j] <- sumry[[i]]$adflevels[[j]]@cval[,3] # adf.test: $p.value
    }
    rownames(sumry[[i]]$adflevelsm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
    colnames(sumry[[i]]$adflevelsm) <- colnames(x[[i]])
    
    
    # Augmented Dickey Fuller Test for first differences
    sumry[[i]]$adfdiff <- apply(x[[i]],2,function(x) ur.df(diff(x),type=type,lags=lags,selectlags=selectlags))
    sumry[[i]]$adfdiffm <- matrix(NA,nrow=switch(type,"none"=3,"trend"=7),ncol=ncol(x[[i]]))
    for (j in 1:length(sumry[[i]]$adflevels)) {
      sumry[[i]]$adfdiffm[switch(type,"none"=1,"trend"=c(1:3)),j] <- sumry[[i]]$adfdiff[[j]]@teststat # adf.test : $statisic
      sumry[[i]]$adfdiffm[switch(type,"none"=2,"trend"=4),j] <- sumry[[i]]$adfdiff[[j]]@lags # adf.test: $parameter
      sumry[[i]]$adfdiffm[switch(type,"none"=3,"trend"=c(5:7)),j] <- sumry[[i]]$adfdiff[[j]]@cval[,3] # adf.test: $p.value
    }
    rownames(sumry[[i]]$adfdiffm) <- switch(type,"trend"=c(rep("Test statistic",3), "Lag order", "p-value-5pct","p-value-5pct","p-value-5pct"),"none"=c("Test statistic", "Lag order", "p-value-5pct"))
    colnames(sumry[[i]]$adfdiffm) <- colnames(x[[i]])
    
    sumry[[i]]$paramcor <- cor(x[[i]])
    sumry[[i]]$diffparamcor <- cor(apply(x[[i]],2,diff))
    
  }
  names(sumry) <- names(x)   
  
  class(sumry) <- "summary.dyntermstrc_param"
  sumry
}


print.summary.dyntermstrc_param <- function(x, ...) {
  for(i in seq_along(x)) {
    cat("---------------------------------------------------\n")
    cat(paste("ADF for ",names(x)[[i]],": ",sep=""))
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    # Augmented Dickey Fuller Test for levels
    print.default(t(x[[i]]$adflevelsm))
    cat("\n")
    cat("---------------------------------------------------\n")
    cat(paste("ADF of differences for ",names(x)[[i]],": ",sep=""))
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    # Augmented Dickey Fuller Test for first differences
    print.default(t(x[[i]]$adfdiffm))
    cat("\n")
    # correlation matrix of parameters
    cat("---------------------------------------------------\n")
    cat(paste("Correlation of parameters for ",names(x)[[i]],": ",sep=""))
    cat("\n")
    cat("---------------------------------------------------\n")
    print.default(x[[i]]$paramcor)
    cat("\n")
    cat("---------------------------------------------------\n")
    cat(paste("Correlation of differences for ",names(x)[[i]],": ",sep=""))
    cat("\n")
    cat("---------------------------------------------------\n")
    print.default(x[[i]]$diffparamcor)
    cat("\n")
  }
}


plot.dyntermstrc_param <- function(x,type="param",...){
  old.par <- par(no.readonly = TRUE) 
  
  
  # 2D plot of parameters
  if(type=="param") {
    if(ncol(x[[1]])<=3) mfrow = c(1,3)
    if(ncol(x[[1]])==4) mfrow = c(2,2)
    if(ncol(x[[1]])>4 && ncol(x[[1]]) <= 6) mfrow = c(2,3)
    
    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)
    
    for(i in seq_along(x)){
      
      param <- x[[i]]
      
      
      plot(param[,1],type="l",xlab="Time",ylab=expression(hat(beta)[0]),
           col=1,lwd=2,... )
      grid()
      plot(param[,2],type="l",xlab="Time",ylab=expression(hat(beta)[1]),
           col=2,lwd=2,... )
      grid()
      plot(param[,3],type="l",xlab="Time",ylab=expression(hat(beta)[2]),
           col=3,lwd=2,... )
      grid()
      
      if(ncol(param)==4) {
        plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
             col=4,lwd=2,... )
        grid()
      }
      
      if(ncol(param)==6) {
        plot(param[,4],type="l",xlab="Time",ylab=expression(hat(tau)[1]),
             col=4,lwd=2,... )
        grid()
        plot(param[,5],type="l",xlab="Time",ylab=expression(hat(beta)[3]),
             col=5,lwd=2,... )
        grid()
        plot(param[,6],type="l",xlab="Time",ylab=expression(hat(tau)[2]),
             col=6,lwd=2,... )
        grid()
      }
    } 
  }
  
  # 2D plot of parameter differences
  if(type=="diffparam") {
    if(ncol(x[[1]])<=3) mfrow = c(1,3)
    if(ncol(x[[1]])==4) mfrow = c(2,2)
    if(ncol(x[[1]])>4 && ncol(x[[1]]) <= 6) mfrow = c(2,3)
    
    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)
    
    for(i in seq_along(x)){
      param <- x[[i]]
      diffparam <- apply(param,2,diff)
      
      for(i in seq(ncol(diffparam))) {
        
        plot(diffparam[,i],type="l",xlab="Time",
             ylab= paste("delta",colnames(diffparam)[i],sep=" "),col=i,lwd=2,... )
        grid()
      }
    }
  }
  
  # ACF/PCF
  if(type=="acf") {
    if(ncol(x[[1]])==3) mfrow = c(2,3)
    if(ncol(x[[1]])==4) mfrow = c(4,2)
    if(ncol(x[[1]])==6) mfrow = c(4,3)
    
    par(mfrow=mfrow,if(length(x)>1) ask=TRUE,...)
    for(i in seq_along(x)){
      param <- x[[i]]
      if(ncol(param) > 3 ){
        for(i in 1:(ncol(param)/2)) acf(param[,i],main=colnames(param)[i])
        for(i in 1:(ncol(param)/2)) pacf(param[,i],main=colnames(param)[i])
        
        for(i in (ncol(param)/2+ 1):ncol(param)) acf(param[,i],main=colnames(param)[i])
        for(i in (ncol(param)/2+ 1):ncol(param)) pacf(param[,i],main=colnames(param)[i])
      } else {
        
        for(i in 1:ncol(param)) acf(param[,i],main=colnames(param)[i])
        for(i in 1:ncol(param)) pacf(param[,i],main=colnames(param)[i])
        
      }
    } 
  }
  
  
  on.exit(par(old.par))
  
  
}

fcontrib <- function(x, method="ns",lambda=0.0609*12, index=1, m=1:10, ylim=NULL ,... ) UseMethod("fcontrib")

fcontrib.dyntermstrc_param <- function(x, method="ns",lambda=0.0609*12, index=1, m=1:10, ylim=NULL ,... ){
  par(if(length(x) > 1) mfrow=length(x),... )
  for(i in seq_along(x)){
    param <- x[[i]]
    if(ncol(param)==3) param <- cbind(param,1/lambda)
    fc1 <- param[,1]
    fc2 <- t(mapply(function(i) param[i,2]*((1-exp(-m/param[i,4]))/(m/param[i,4])), seq(nrow(param))))
    fc3 <- t(mapply(function(i) param[i,3]*(((1-exp(-m/param[i,4]))/(m/param[i,4]))-exp(-m/param[i,4])), seq(nrow(param))))
    if(ncol(param)==6) fc4 <- t(mapply(function(i) param[i,5]*(((1 - exp(-m/param[i,6]))/(m/param[i,6])) - exp(-m/param[i,6])), seq(nrow(param))))
    
    if(is.null(ylim)) ylim <- c(min(fc1,fc2,fc3),max(fc1,fc2,fc3))          
    
    plot(m,rep(fc1[index],length(m)), col=1,type="l",lty=1, ylim=ylim, xlab="Time to maturity", ylab="Factor contribution",lwd=2,main=get_realnames(method))
    
    
    ## beta_1*( )
    lines(m, fc2[index,],type="l",col=2,lty=3,lwd=2)
    ## beta_2*()
    lines(m, fc3[index,],lty=4,col=4,lwd=2)
    ## beta_3*()
    if(ncol(param)==6) lines(m, fc4[index,],lty=5,col=5,lwd=2)
    
    abline(h=0,lty=1, lwd = 1, col = "grey")
    
    legend("topright",bg='white', box.col="white", box.lwd = 0,
           legend=c(expression(beta[0]),
                    expression(beta[1]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1])))),
                    expression(beta[2]*(frac(1-exp(-frac(m,tau[1])),frac(m,tau[1]))-exp(-frac(m,tau[1])))),
                    if(method=="sv") expression(beta[3]*(frac(1-exp(-frac(m,tau[2])),frac(m,tau[2]))-exp(-frac(m,tau[2])))),
                    if(method=="asv")  expression(beta[3]*(frac(1-exp(-frac(m,tau[2])),frac(m,tau[2]))-exp(-frac(2*m,tau[2]))))                 
           ),
           lty=c(1,2,4,5), col=c(1,2,4,5),bty="o", cex = 0.9
    ) 
    box()  
  }
}

print.dyntermstrc_yields <- function(x, ...){
  cat("---------------------------------------------------\n")
  cat("Estimated",get_realnames(x$method), "parameters:")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  cat("Number of oberservations:",nrow(x$optparam),"\n")
  cat("\n")
  tsparam <- param.dyntermstrc_yields(x)
  print.default(lapply(tsparam,summary.default))
  cat("\n")
}

summary.dyntermstrc_yields <- function(object, ...){
  x <- object
  #
  #Change Mean and sqrt, da erst RMSE pro Zeile berechnet werden muss und dann der Mean 
  #y_mrsme <-  sqrt(mean(apply((x$yields-x$yhat)^2,1,mean))) # Original
  y_mrsme <-  mean(sqrt(apply((x$yields-x$yhat)^2,1,mean)))
  y_maabse <- mean(apply(abs(x$yields-x$yhat),1,mean))
  y_smape <- mean(apply(smape1(x$yields,x$yhat),1,mean))
  
  
  sumry <- list()
  sumry$gof <- rbind(y_mrsme,y_maabse,y_smape)
  rownames(sumry$gof) <- c("RMSE-Yields (in %)", "MAE-Yields (in %)",  "sMAPE-Yields (in %)")
  
  if (object$method != "dl") {
    ## extract convergence info
    for (i in (1:length(x$opt_result))) {
      sumry$convergence[i] <- x$opt_result[[i]]$convergence
      sumry$solvermsg <- x$opt_result[[i]]$message
    }
  }  
  class(sumry) <- "summary.dyntermstrc_yields"
  sumry
}

print.summary.dyntermstrc_yields <- function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Goodness of fit:\n")
  cat("---------------------------------------------------\n")
  
  print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
  
  if (length(x) > 1) {
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information from optim ():\n")
    cat("---------------------------------------------------\n")
    
    print.default(t(as.matrix(x$convergence)))
  }
}

plot.dyntermstrc_yields <- function(x,...)
{
  ## plot estimated yield curves in 3D
  
  Z <- matrix(nrow=nrow(x$optparam),ncol=length(x$maturities))# OK
  for (i in 1:nrow(x$optparam)){
    Z[i,] <- spotrates(x$method,x$optparam[i,],x$maturities, x$lambda)
  }
  
  X <- 1:nrow(Z)
  Y <- x$maturities
  
  open3d()
  persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Time", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
  
}

## print, summary and plot method for objects of the class "termstrc_cs"

print.termstrc_cs <- function(x,...) {
  cat("---------------------------------------------------\n")
  if(x$rse) cat("Estimated parameters and robust standard errors:\n") else
    cat("Estimated parameters and standard errors:\n") 
  cat("---------------------------------------------------\n")
  for(i in seq(x$n_group)) {
    print(paste(names(x$alpha)[[i]],":",sep=""))
    cs_coef <- if(x$rse) coeftest(x$regout[[i]],vcov=vcov) else summary(x$regout[[i]])
    if(x$rse) rownames(cs_coef) <- paste("alpha",c(seq_along(x$alpha[[i]]))) else
      rownames(cs_coef$coefficients) <- paste("alpha",c(seq_along(x$alpha[[i]])))
    print(cs_coef)
    cat("\n")
    x
  }
}

summary.termstrc_cs <-
  function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
    
    sMAPE_p <- mapply(function(i) smape(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    sMAPE_y <- mapply(function(i) smape(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
    
    wRMSE_p <- mapply(function(i) rmsewe(x$p[[i]],x$phat[[i]],x$duration[[i]][,3]),seq(x$n_group))
    wRMSE_y <- mapply(function(i) rmsewe(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100,x$duration[[i]][,3]),seq(x$n_group))
    
    gof <- rbind(RMSE_p,AABSE_p,sMAPE_p,wRMSE_p,RMSE_y,AABSE_y,sMAPE_y,wRMSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices","Weighted RMS-Prices", "RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMS-Yields (in %)")
    sumry <- list(gof)
    names(sumry) <- c("gof")
    class(sumry) <- "summary.termstrc_cs"
    sumry
  } 

print.summary.termstrc_cs <-
  function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
    cat("\n")
    x$gof
    
  }

plot.termstrc_cs <-
  function(x,matrange =c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                         max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))),
           multiple=FALSE, ctype="spot",
           lwd=2,lty=1,type="l",errors="none",inset=c(0.1,0.3),ask=TRUE, ...) {
    
    # min and max maturity of all bonds in the sample 
    samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                   max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))))     
    
    
    cdata <- switch(ctype, "spot" = x$spot,
                    "forward" = x$forward,
                    "discount" = x$discount
    )
    
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
                    "forward" = "Forward rate curve",
                    "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (multiple) {
      
      plot(x=cdata,multiple=multiple, expoints=NULL,lwd=lwd,type=type,...) }
    
    if (!multiple && ctype %in% c("spot", "forward", "discount")){
      old.par <- par(no.readonly = TRUE)
      if(x$n_group !=1) par(ask=ask)
      
      # plot each interest rate curve seperately
      for (k in seq(x$n_group)  ) 
      {
        
        plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
                      xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
                             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,ylab=cname, type=type,...)
        
        
        
        title(x$group[k])
        
        if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
          # lower ci         
          lines(cdata[[k]][,1],cdata[[k]][,3]*100, type="l", lty=3, col="steelblue" )   
          # upper ci 
          lines(cdata[[k]][,1],cdata[[k]][,4]*100, type="l", lty=3, col="steelblue")
          # knot points 
          abline(v=c(x$knotpoints[[k]]),lty=2, col="darkgrey")
          legend("bottom",legend=c("Zero-coupon yield curve",
                                   if(x$rse) "95 % Confidence interval (robust s.e.)" else "95 % Confidence interval" ,"Yield-to-maturity", "Knot points"),
                 col=c("steelblue","steelblue","red", "darkgrey"),
                 lty = c(1,3,-1,2), pch=c(-1,-1,21,-1))
          
        } else  legend("bottom",legend=cname,col=c("steelblue"), lty = lty , pch=(-1))
        
        
      }
      on.exit(par(old.par))
    }
    
    # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=NULL,
                                xlim= c(max(floor(samplemat[1]),matrange[1]),
                                        min(ceiling(samplemat[2]),matrange[2],max(mapply(function(i) 
                                          max(x$spread[[i]][,1]),seq(x$spread))))),lwd=lwd ,...) 
    }
    
    
    # plot errors 
    if(errors %in% c("price", "yield")){
      
      edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )
      
      if(x$n_group == 1) ask= FALSE
      
      for(k in seq(x$n_group)){
        plot.error(edata[[k]],ask=ask
                   ,main=x$group[k],ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
        
        legend("topright", legend=c(paste("RMSE",
                                          switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                                                 "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                                    paste("MAE",switch(errors,"price" = round(aabse(x$p[[k]],x$phat[[k]]),4),
                                                       "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),sep=": ")),bty="n", inset=c(-0.1,0),xpd=TRUE) 
        
      }
      
    }						
    
  } 
print.termstrc_nss <- 
  function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Estimated",get_realnames(x$method), "parameters:")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    parameters <- mapply(function(i) x$opt_result[[i]]$par,seq_along(x$opt_result))
    colnames(parameters) <- names(x$opt_result)
    rownames(parameters) <- get_paramnames(x$method)
    print.default(format(parameters,digits=6,scientific=FALSE),quote=FALSE)
    cat("\n")
    
  }

plot.termstrc_nss <-
  function(x,matrange=c(min(mapply(function(i) min(x$y[[i]][,1]),seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))))
           ,multiple=FALSE, expoints=unlist(x$expoints), ctype="spot",
           errors="none",
           lwd=2,lty=1,type="l",inset=c(0.8,0.1),ask=TRUE,
           ...) {
    
    # min and max maturity of all bonds in the sample 
    samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                   max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))) 
    
    
    cdata <- switch(ctype, "spot" = x$spot,
                    "forward" = x$forward,
                    "discount" = x$discount
    )
    
    
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
                    "forward" = "Forward rate curve",
                    "discount" = "Discount factor curve" )
    
    
    # plot all interest rate curves together
    if (multiple) {
      
      plot(x=cdata,multiple=multiple, expoints=expoints,lwd=lwd,type=type,...) }
    
    if (!multiple && ctype %in% c("spot", "forward", "discount")){
      old.par <- par(no.readonly = TRUE)
      
      if(x$n_group != 1) par(ask=ask)
      
      # plot each interest rate curve seperately
      for (k in seq(x$n_group)  ) 
      {
        
        plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
                      xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
                             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,ylab=cname, type=type, ...
        )
        
        title(names(x$opt_result)[k])
        
        if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
          legend("bottom",legend=c("Zero-coupon yield curve","Yield-to-maturity"),
                 col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))}
        else 	legend("bottom",legend=cname	,col=c("steelblue"), lty = lty , pch=(-1))
        
        
      }     
      on.exit(par(old.par))
    }
    
    # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=expoints,
                                xlim= c(max(floor(samplemat[1]),matrange[1]),
                                        min(ceiling(samplemat[2]),matrange[2])),lwd=lwd,
                                ...)
    }
    # plot errors 
    if(errors %in% c("price", "yield")){
      
      edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )
      if(x$n_group == 1) ask= FALSE
      for(k in seq(x$n_group)){
        
        plot.error(edata[[k]],ask=ask,main=x$group[k],
                   ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
        
        legend("topright", legend=c(paste("  RMSE",
                                          switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                                                 "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                                    paste("MAE",switch(errors,"price" = round(aabse(x$p[[k]],x$phat[[k]]),4),
                                                       "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),sep=": ")),bty="n", inset=c(-0.1,0),xpd=TRUE) 
        
      }
      
    }
    
    
    
  }  

summary.termstrc_nss <- function(object,...) {
  x <- object
  RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
  AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]) ,seq(x$n_group))
  RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100) ,seq(x$n_group))
  AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100) ,seq(x$n_group))
  
  sMAPE_p <- mapply(function(i) smape(x$p[[i]],x$phat[[i]]),seq(x$n_group))
  sMAPE_y <- mapply(function(i) smape(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100),seq(x$n_group))
  
  wRMSE_p <- mapply(function(i) rmsewe(x$p[[i]],x$phat[[i]],x$duration[[i]][,3]),seq(x$n_group))
  wRMSE_y <- mapply(function(i) rmsewe(x$y[[i]][,2]*100,x$yhat[[i]][,2]*100,x$duration[[i]][,3]),seq(x$n_group))
  
  gof <- rbind(RMSE_p,AABSE_p,sMAPE_p,wRMSE_p,RMSE_y,AABSE_y,sMAPE_y,wRMSE_y)
  colnames(gof) <- names(x$p)
  rownames(gof) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices", "Weighted RMSE-Prices" ,"RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMSE-Yields")
  convergencegroup <- as.matrix(mapply(function(i) x$opt_result[[i]]$convergence,
                                       seq_along(x$opt_result)))
  
  colnames(convergencegroup) <- "optim() convergence info"
  rownames(convergencegroup) <- x$group
  convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,seq_along(x$opt_result)))
  colnames(convergence) <- "optim() solver message"
  rownames(convergence) <- x$group
  sumry <- list(gof,convergencegroup,convergence,startparam=x$startparam)
  names(sumry) <- c("gof","convergencegroup","convergence","startparam")
  class(sumry) <- "summary.termstrc_nss"
  sumry
}

print.summary.termstrc_nss <-
  function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Startparameters:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(format(x$startparam,digits=6),quote=FALSE)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$convergencegroup)
    cat("\n")
    print.default(x$convergence)
    cat("\n")
    cat("\n")
    x
  }
### Nelson/Siegel and Svensson spot curve estimaton from zero yields 

zeroyields <- function(maturities, yields, dates)
{
  zy <- list(maturities = maturities, yields = yields, dates = dates)
  class(zy) <- "zeroyields"
  zy
}

print.zeroyields <- function(x, ...)
{
  cat("This is a data set of zero-coupon yields.\n")
  cat(paste("Maturities range from", min(x$maturities), "to", max(x$maturities),"years.\n"))
  cat(paste("There are",nrow(x$yields), "observations between",x$dates[1], "and",x$dates[length(x$dates)],".\n"))
}

summary.zeroyields <- function(object, ...)
{
  print(summary(object$yields))
}

plot.zeroyields <- function(x,...)
{
  Z <- as.matrix(x$yields)
  X <- 1:nrow(x$yields)
  Y <- x$maturities
  
  open3d()
  persp3d(X, Y, Z, col = "green3", box = FALSE,xlab = "Time", ylab = "Maturities (years)", zlab = "Zero-yields (%)")
}

### Spot rates functions

## Nelson/Siegel spot rate function
spr_ns <- function(beta, m){
  (beta[1] + beta[2]*((1-exp(-m/beta[4]))/(m/beta[4])) +
     beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4])))
}

## Svensson spot rate function 
spr_sv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
     beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
     beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))
}

## Adjusted Svensson spot rate function
spr_asv <- function(beta, m){
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
     beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
     beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-(2*m)/beta[6])))
}

## Diebold/Li spot rate function
spr_dl <- function(beta,m,lambda){
  (beta[1] + beta[2]*((1-exp(-m*lambda))/(m*lambda))+
     beta[3]*(((1-exp(-m*lambda))/(m*lambda))-exp(-m*lambda)))
}

## Spot rate wrapper function
spotrates <- function(method,beta,m,lambda = 0.0609*12){ 
  switch(method,
         "ns" = spr_ns(beta,m),
         "sv" = spr_sv(beta,m),
         "asv"= spr_asv(beta,m),
         "dl" = spr_dl(beta,m,lambda))
}

### Forward rate functions

## Nelson/Siegel forward rate function
fwr_ns <- function(beta,m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
     beta[3]*(m/beta[4]*exp(-m/beta[4])))
}

## Svensson forward rate function
fwr_sv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
     beta[3] *m/beta[4]*exp(-m/beta[4]) +
     beta[5] *m/beta[6]*exp(-m/beta[6]))
}

## Adjusted Svensson forward rate function
fwr_asv <- function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
     beta[3] *m/beta[4]*exp(-m/beta[4]) +
     beta[5] *(exp(-m/beta[6])+(2*m/beta[6] -1)*exp(-2*m/beta[6])))
}

## Diebold/Li forward rate function
fwr_dl <- function(beta, m,lambda) {
  (beta[1] + beta[2]*exp(-m*lambda)
   + beta[3]*(m*lambda*exp(-m*lambda)))
}

## Forward rates wrapper function
forwardrates <- function(method,beta,m,lambda){
  switch(method,
         "ns" = fwr_ns(beta,m),
         "sv" = fwr_sv(beta,m),
         "asv"= fwr_asv(beta,m),
         "dl"= fwr_dl(beta,m,lambda))
}

## Implied forward rates calculation
impl_fwr <- function(m,s) {
  impl_fwr <- c(s[1],(s[-1]*m[-1] - s[-length(s)]*m[-length(m)])/(diff(m)))
  impl_fwr[1] <- impl_fwr[2]
  impl_fwr	
}

get_paramnames <- function(method){
  names <- c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")
  switch(method,"ns"= names[1:4],"sv"=names,"asv"=names,"dl"=names[1:3])
}

get_realnames <- function(method){
  switch(method,"dl"="Diebold/Li","ns"="Nelson/Siegel","sv"="Svensson","asv"="Adj. Svensson")
}

### Loss function for parametric methods
get_objfct <- function(method) {
  objfct <- switch(method,
                   "dl" = objfct_dl,
                   "ns" = objfct_ns,
                   "sv" = objfct_sv,
                   "asv" = objfct_asv)
}

### Gradient of loss function for parametric methods
get_grad_objfct <- function(method) {
  grad_objfct <- switch(method,
                        "dl" = grad_dl,
                        "ns" = grad_ns,
                        "sv" = grad_sv,
                        "asv" = grad_asv)
}

### Diebold/Li loss function for yields
objfct_dl <- function(beta, lambda, m, y)
{
  sum((y - spr_dl(beta,m, lambda))^2)
}

### Nelson/Siegel loss function for yields
objfct_ns <- function(beta, m, y)
{
  sum((y - spr_ns(beta,m))^2)
}

### Nelson/Siegel grid loss function for yields
objfct_ns_grid <- function(beta, tau, m, y)
{
  betans <- c(beta, tau)
  sum((y - spr_ns(betans,m))^2)
}

### Svensson loss function for yields
objfct_sv <- function(beta, m, y)
{
  sum((y - spr_sv(beta,m))^2)
}


### Svensson grid loss function for yields
objfct_sv_grid <- function(beta, tau,  m, y)
{
  betasv <- c(beta[1:3], tau[1], beta[4], tau[2])
  sum((y - spr_sv(betasv,m))^2)
}

### Adjusted Svensson loss function for yields
objfct_asv <- function(beta, m, y)
{
  sum((y - spr_asv(beta,m))^2)
}


### Adjusted Svensson grid loss function for yields
objfct_asv_grid <- function(beta, tau, m, y)
{
  betasv <- c(beta[1:3], tau[1], beta[4], tau[2])
  sum((y - spr_asv(betasv,m))^2)
}

### Constraints for constrOptim()

get_constraints <- function(method, tauconstr) {
  
  ## tauconstr = c(upper, lower, gridsize, distance)
  
  ## Diebold/Li
  
  if (method == "dl") {
    ui <- rbind(c(1,0,0),               # beta0 > 0
                c(1,1,0))               # beta0 + beta1 > 0
    ci <- c(0,0)
  }
  
  ## Nelson/Siegel
  
  if (method == "ns") {
    ui <- rbind(c(1,0,0,0),             # beta0 > 0
                c(1,1,0,0),             # beta0 + beta1 > 0
                c(0,0,0,1),             # tau1 > tauconstr[1]
                c(0,0,0,-1))            # tau1 <= tauconstr[2]
    ci <- c(0,0,tauconstr[1],-tauconstr[2])### HILFE###
  }
  
  ## Svensson
  
  if (method == "sv") {
    ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                c(0,0,0,1,0,0),        # tau1 > tauconstr[1]
                c(0,0,0,-1,0,0),       # tau1 <= tauconstr[2]
                c(0,0,0,0,0,1),        # tau2 > tauconstr[1]
                c(0,0,0,0,0,-1),       # tau1 <= tauconstr[2]
                c(0,0,0,-1,0,1))       # tau2 - tau1 > tauconstr[4] # Minimal distance between both
    ci <- c(0,0,tauconstr[1],-tauconstr[2],tauconstr[1],-tauconstr[2],tauconstr[4]) ### HILFE###
  }
  
  ## Adjusted Svensson
  
  if (method == "asv") {
    ui <- rbind(c(1,0,0,0,0,0),        # beta0 > 0
                c(1,1,0,0,0,0),        # beta0 + beta1 > 0
                c(0,0,0,1,0,0),        # tau1 > tauconstr[1]
                c(0,0,0,-1,0,0),       # tau1 <= tauconstr[2]
                c(0,0,0,0,0,1),        # tau2 > tauconstr[1]
                c(0,0,0,0,0,-1),       # tau1 <= tauconstr[2]
                c(0,0,0,-1,0,1))       # tau2 - tau1 > 0 #tau2 on the right side of tau
    ci <- c(0,0,tauconstr[1],-tauconstr[2],tauconstr[1],-tauconstr[2],0) ### HILFE###
  }
  
  constraints <- list(ui = ui, ci = ci)
  constraints
}

### Loss function for parametric methods
get_objfct_bonds <- function(method) {
  objfct <- switch(method,
                   "dl" = objfct_dl_bonds,
                   "ns" = objfct_ns_bonds,
                   "sv" = objfct_sv_bonds,
                   "asv" = objfct_asv_bonds)
}

### Gradient of loss function for parametric methods
get_grad_objfct_bonds <- function(method) {
  grad_objfct <- switch(method,
                        "dl" = grad_dl_bonds,
                        "ns" = grad_ns_bonds,
                        "sv" = grad_sv_bonds,
                        "asv" = grad_asv_bonds)
}

### Diebold/Li loss function for bonds 
objfct_dl_bonds <- function(beta, lambda, m, cf, w, p) {
  phat <- bond_prices2("dl",beta,m,cf, lambda)$bond_prices
  sum(w*((p - phat)^2))
}

### Nelson/Siegel loss function for bonds 
objfct_ns_bonds <- function(beta, m, cf, w, p) {
  .Call("objfct_ns_bonds_Cpp", beta, m, cf, w, p)
}

### Nelson/Siegel grid loss function for bonds

objfct_ns_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("objfct_ns_bonds_gridCpp", beta, tau, m, cf, w, p)
}

### Svensson loss function for bonds 
objfct_sv_bonds <- function(beta, m, cf, w, p) {
  .Call("objfct_sv_bonds_Cpp", beta, m, cf, w, p)
}

### Svensson grid loss function for bonds
objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("objfct_sv_bonds_gridCpp", beta, tau, m, cf, w, p)
}

### Adjusted Svensson loss function for bonds 
objfct_asv_bonds <- function(beta, m, cf, w, p) {
  .Call("objfct_asv_bonds_Cpp", beta, m, cf, w, p)
}

### Adjusted Svensson grid loss function for bonds
objfct_asv_bonds_grid <- function(beta, tau, m, cf, w, p) {
  .Call("objfct_asv_bonds_gridCpp", beta, tau, m, cf, w, p)
}



### Differenz zwischen zwei Tagen
Datumsdiff <- function(wert, wert1) {
  if(wert < wert1){
    y = wert
    wert = wert1 
    wert1 = y}
  wert <- as.Date(wert, format="%Y-%m-%d")
  wert1 <- as.Date(wert1, format="%Y-%m-%d")
  save <- as.numeric(difftime(as.POSIXct(wert), as.POSIXct(wert1), units="days"))
  return(save)
}

##Erzeugung der Rohliste (Dynamisch)

Listeerstelltung <- function(wert,country){
  
  Datensatz <- list()
  
  #Erzeugung der generellen Listenstruktur
  init_liste <- list()
  init_liste[[""]] <- list()
  init_liste[[1]][[country]] <- list()
  init_liste[[1]][[country]]$"ISIN" <-"GB1111111111111"
  init_liste[[1]][[country]]$"MATURITYDATE" <- as.Date("2022-08-05")
  init_liste[[1]][[country]]$"ISSUEDATE" <- as.Date("2022-08-05")
  init_liste[[1]][[country]]$"COUPONRATE" <- 0.21
  init_liste[[1]][[country]]$"PRICE" <- 0.21
  init_liste[[1]][[country]]$"ACCRUED" <- 0.21
  init_liste[[1]][[country]]$"CASHFLOWS" <- list()
  init_liste[[1]][[country]][["CASHFLOWS"]]$"ISIN" <-"GB1111111111111"
  init_liste[[1]][[country]][["CASHFLOWS"]]$"CF" <- 100.00
  init_liste[[1]][[country]][["CASHFLOWS"]]$"DATE" <- as.Date("2022-08-05")
  init_liste[[1]][[country]]$"TODAY" <- as.Date("2022-08-05")
  init_liste[[1]][[country]]$"SETTLEMENT DAY" <- as.Date("2022-08-05")
  
  
  
  for (i in 1:wert){
    Datensatz[i] <- init_liste
    class(Datensatz[i]) <- "couponbonds"
  }
  
  class(Datensatz) <- "dyncouponbonds"
  return(Datensatz)
  
}

###Hinzufuegen einer Anleihe zum Datensatz

add_bond <- function(bonddata, group, ISIN, Maturitydate, issuedate, couponrate, price, accrued, valuationdate) UseMethod("add_bond")

## static data set
add_bond.couponbonds <- function(bonddata, group, ISIN, Maturitydate, issuedate, couponrate, price, accrued, valuationdate){
  
  high_index <- length(bonddata[[group]]$ISIN)
  bonddata[[group]]$ISIN[(high_index+1)] <- ISIN
  bonddata[[group]]$MATURITYDATE[(high_index+1)] <- Maturitydate
  bonddata[[group]]$ISSUEDATE[(high_index+1)] <- issuedate
  bonddata[[group]]$COUPONRATE[(high_index+1)] <- couponrate
  bonddata[[group]]$PRICE[(high_index+1)] <- price
  bonddata[[group]]$ACCRUED[(high_index+1)] <- accrued
  
  
  
  class(bonddata) <- "couponbonds"
  return(bonddata)
  #bonddata
}

## add bonds from a dynamic bonddata set 
add_bond.dyncouponbonds <- function(bonddata, group, ISIN, Maturitydate, issuedate, couponrate, price, accrued, valuationdate) {
  
  for (i in seq(length(bonddata))) {
    if (bonddata[[i]][[group]]$TODAY == valuationdate){
      search_index <- i
    }
  }
  high_index <- length(bonddata[[search_index]][[group]]$ISIN)
  bonddata[[search_index]][[group]]$ISIN[(high_index+1)] <- ISIN
  bonddata[[search_index]][[group]]$MATURITYDATE[(high_index+1)] <- Maturitydate
  bonddata[[search_index]][[group]]$ISSUEDATE[(high_index+1)] <- issuedate
  bonddata[[search_index]][[group]]$COUPONRATE[(high_index+1)] <- couponrate
  bonddata[[search_index]][[group]]$PRICE[(high_index+1)] <- price
  bonddata[[search_index]][[group]]$ACCRUED[(high_index+1)] <- accrued
  
  
  class(bonddata) <- "dyncouponbonds"
  return(bonddata)
  
}

### Erzeuge der Daten fuer zukünftige Cashflows
add_cashflows <- function(bonddata, group, ISIN, Maturitydate, couponrate, divday1, divmonth1, divday2, divmonth2, valuationdate, accrued, result) UseMethod("add_cashflows")

add_cashflows.couponbonds <- function(bonddata, group, ISIN, Maturitydate, couponrate, divday1, divmonth1, divday2, divmonth2, valuationdate, accrued, result){
  ### Bonddata only includes the data for one specific date, which makes it easier to go through the dataset 
  Y1 <- year(valuationdate)
  Y2 <- year(Maturitydate)
  diffyear <- Y2-Y1
  
  dates1 <- as.Date(sapply(c(Y1:Y2),getcoupondate, x=divday1,y=divmonth1))
  dates2 <- as.Date(sapply(c(Y1:Y2),getcoupondate, x=divday2,y=divmonth2))
  dates <- sort(c(dates1, dates2))
  
  if(Maturitydate < dates[length(dates)]){
    dates <- dates[-length(dates)]
    ISINvector <- rep(ISIN,2*(diffyear+1)-1)
    CFvector1 <- rep(100*couponrate/2,2*(diffyear+1)-2)
    CFvector2 <- c((couponrate*100)/2+100)
    CFvector <- c(CFvector1,CFvector2)
  } else {
    
    ISINvector <- rep(ISIN,2*(diffyear+1))
    CFvector1 <- rep(100*couponrate/2,2*(diffyear+1)-1)
    CFvector2 <- c((couponrate*100)/2+100)
    CFvector <- c(CFvector1,CFvector2)
  }
  
  
  if (valuationdate < getcoupondate(divday1,divmonth1,Y1)) {
    # In diesem Fall passiert nichts 
    a = 1
    
  } else if ((valuationdate >= getcoupondate(divday1,divmonth1,Y1))&(valuationdate < getcoupondate(divday2,divmonth2,Y1))) {
    # In diesem Fall muss das erste Datum geloescht werden
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
    
  } else if (valuationdate >= getcoupondate(divday2,divmonth2,Y1)) {
    # In diesem Fall muessen die ersten beiden Eintraege geloescht werden
    ISINvector <- ISINvector[-c(1:2)]
    dates <- dates[-c(1:2)]
    CFvector <- CFvector[-c(1:2)]
  }
  
  # In diesem Fall befinde ich mich zwischen Coupon Tag und Ex-Dividenden Tag. Coupon Tag ist eingeschlossen.
  if (accrued <=0){
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } else if ((valuationdate < getcoupondate(divday2,divmonth2,Y1))&(NextBuisnessday(valuationdate,delayday) > getcoupondate(divday2,divmonth2,Y1))){
    if (result$index=="Neither"){
      a = 1
    } else {
      a = 1
      if ((result$index == "Long") & (result$first=="Yes")){
        result$first="No"
      }
    } 
    
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } else if ((valuationdate < getcoupondate(divday1,divmonth1,Y1))&(NextBuisnessday(valuationdate,delayday) > getcoupondate(divday1,divmonth1,Y1))){
    if (result$index=="Neither"){
      a = 1
    } else {
      a = 1
      if ((result$index == "Long") & (result$first=="Yes")){
        result$first="No"
      }
    } 
    
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } 
  
  ### Short and long first Dividend Periods 
  if (result$index=="Neither"){
    a = 1
  } else {
    if (result$index == "Short"){
      CFvector[1] <- result$adjustedcashflow
    } else if ((result$index == "Long") & (result$first=="Yes")){
      ISINvector <- ISINvector[-1]
      dates <- dates[-1]
      CFvector <- CFvector[-1]
      CFvector[1] <- result$adjustedcashflow
    } else if ((result$index == "Long") & (result$first=="No")){
      CFvector[1] <- result$adjustedcashflow
    }
  }
  
  
  
  
  
  
  
  bonddata[[group]][["CASHFLOWS"]]$"ISIN" <- append(bonddata[[group]][["CASHFLOWS"]]$"ISIN", ISINvector)
  bonddata[[group]][["CASHFLOWS"]]$"CF" <- append(bonddata[[group]][["CASHFLOWS"]]$"CF", CFvector)
  bonddata[[group]][["CASHFLOWS"]]$"DATE" <- append(bonddata[[group]][["CASHFLOWS"]]$"DATE", dates)
  
  class(bonddata) <- "couponbonds"
  
  return(bonddata)
}

add_cashflows.dyncouponbonds <- function(bonddata, group, ISIN, Maturitydate, couponrate, divday1, divmonth1, divday2, divmonth2, valuationdate, accrued, result){
  ### Bonddata only includes the data for one specific date, which makes it easier to go through the dataset 
  Y1 <- year(valuationdate)
  Y2 <- year(Maturitydate)
  diffyear <- Y2-Y1
  
  dates1 <- as.Date(sapply(c(Y1:Y2),getcoupondate, x=divday1,y=divmonth1))
  dates2 <- as.Date(sapply(c(Y1:Y2),getcoupondate, x=divday2,y=divmonth2))
  dates <- sort(c(dates1, dates2))
  
  if(Maturitydate < dates[length(dates)]){
    dates <- dates[-length(dates)]
    ISINvector <- rep(ISIN,2*(diffyear+1)-1)
    CFvector1 <- rep(100*couponrate/2,2*(diffyear+1)-2)
    CFvector2 <- c((couponrate*100)/2+100)
    CFvector <- c(CFvector1,CFvector2)
  } else {
    
    ISINvector <- rep(ISIN,2*(diffyear+1))
    CFvector1 <- rep(100*couponrate/2,2*(diffyear+1)-1)
    CFvector2 <- c((couponrate*100)/2+100)
    CFvector <- c(CFvector1,CFvector2)
  }
  
  
  if (valuationdate < getcoupondate(divday1,divmonth1,Y1)) {
    # In diesem Fall passiert nichts 
    a = 1
    
  } else if ((valuationdate >= getcoupondate(divday1,divmonth1,Y1))&(valuationdate < getcoupondate(divday2,divmonth2,Y1))) {
    # In diesem Fall muss das erste Datum geloescht werden
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
    
  } else if (valuationdate >= getcoupondate(divday2,divmonth2,Y1)) {
    # In diesem Fall muessen die ersten beiden Eintraege geloescht werden
    ISINvector <- ISINvector[-c(1:2)]
    dates <- dates[-c(1:2)]
    CFvector <- CFvector[-c(1:2)]
  }
  
  if (accrued <=0){
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } else if ((valuationdate < getcoupondate(divday2,divmonth2,Y1))&(NextBuisnessday(valuationdate,delayday) > getcoupondate(divday2,divmonth2,Y1))){
    if (result$index=="Neither"){
      a = 1
    } else {
      a = 1
      if ((result$index == "Long") & (result$first=="Yes")){
        result$first="No"
      }
    } 
    
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } else if ((valuationdate < getcoupondate(divday1,divmonth1,Y1))&(NextBuisnessday(valuationdate,delayday) > getcoupondate(divday1,divmonth1,Y1))){
    if (result$index=="Neither"){
      a = 1
    } else {
      a = 1
      if ((result$index == "Long") & (result$first=="Yes")){
        result$first="No"
      }
    } 
    ISINvector <- ISINvector[-1]
    dates <- dates[-1]
    CFvector <- CFvector[-1]
  } 
  
  if (result$index=="Neither"){
    a = 1
  } else {
    if (result$index == "Short"){
      CFvector[1] <- result$adjustedcashflow
    } else if ((result$index == "Long") & (result$first=="Yes")){
      ISINvector <- ISINvector[-1]
      dates <- dates[-1]
      CFvector <- CFvector[-1]
      CFvector[1] <- result$adjustedcashflow
    } else if ((result$index == "Long") & (result$first=="No")){
      CFvector[1] <- result$adjustedcashflow
    }
  }
  
  for (i in seq(length(bonddata))) {
    if (bonddata[[i]][[group]]$TODAY == valuationdate){
      search_index <- i
    }
  }
  
  bonddata[[search_index]][[group]][["CASHFLOWS"]]$"ISIN" <- append(bonddata[[search_index]][[group]][["CASHFLOWS"]]$"ISIN", ISINvector)
  bonddata[[search_index]][[group]][["CASHFLOWS"]]$"CF" <- append(bonddata[[search_index]][[group]][["CASHFLOWS"]]$"CF", CFvector)
  bonddata[[search_index]][[group]][["CASHFLOWS"]]$"DATE" <- append(bonddata[[search_index]][[group]][["CASHFLOWS"]]$"DATE", dates)
  
  class(bonddata) <- "dyncouponbonds"
  return(bonddata)
}

adjustmentrate <- function(bonddata,ISIN,todaydate,issuedate,couponrate,divday1,divmonth1,divday2,divmonth2,settlementday){
  
  if (!((day(issuedate)==as.numeric(divday1)) & (month(issuedate)==as.numeric(divmonth1))) | ((day(issuedate)==as.numeric(divday2)) & (month(issuedate)==as.numeric(divmonth2)))){
    ##Check if we are still in the first ever coupon period: 
    Y1 <- year(todaydate)
    Y2 <- year(issuedate)
    diffyear <- Y2-Y1
    
    ## Abstand zwischen Issuedate und Todaydate muss kleiner als 2 Jahre sein
    if (abs(diffyear) < 2){
      
      dates1 <- as.Date(sapply(c(Y2:(Y2+1)),getcoupondate, x=divday1,y=divmonth1))
      dates2 <- as.Date(sapply(c(Y2:(Y2+1)),getcoupondate, x=divday2,y=divmonth2))
      dates <- sort(c(dates1, dates2))
      
      if (diffyear<0){
        if (issuedate < getcoupondate(divday1,divmonth1,Y2)) {
          # In diesem Fall passiert nichts
          a = 1
        } else if ((issuedate > getcoupondate(divday1,divmonth1,Y2))&(issuedate < getcoupondate(divday2,divmonth2,Y2))) {
          # In diesem Fall muss das erste Datum geloescht werden
          dates <- dates[-1]
        } else if (issuedate > getcoupondate(divday2,divmonth2,Y2)) {
          # In diesem Fall muessen die ersten beiden Eintraege geloescht werden
          dates <- dates[-c(1:2)]
        }
      } else {
        if (issuedate < getcoupondate(divday1,divmonth1,Y1)) {
          # In diesem Fall passiert nichts 
          a = 1
        } else if ((issuedate > getcoupondate(divday1,divmonth1,Y1))&(issuedate < getcoupondate(divday2,divmonth2,Y1))) {
          # In diesem Fall muss das erste Datum geloescht werden
          dates <- dates[-1]
        } else if (issuedate > getcoupondate(divday2,divmonth2,Y1)) {
          # In diesem Fall muessen die ersten beiden Eintraege geloescht werden
          dates <- dates[-c(1:2)]
        }
      }
      
      
      if (todaydate <= dates[1]){
        result <- list()
        #in First Coupon period
        search1 <- which(bonddata$`Close of Business Date`==offset(dates[1],-4,"QuantLib/UnitedKingdom/Exchange"))
        search2 <- which(bonddata$ISIN==ISIN)
        combine <- c(search1,search2)
        
        if(!duplicated(c(search1,search2))[which(duplicated(c(search1,search2))==TRUE)]){
          print("Fehler in der Routine")
        } else if (length(duplicated(c(search1,search2))[which(duplicated(c(search1,search2))==TRUE)])!=1){
          print("Mehr als ein Treffer")
        }
        
        if (bonddata$`Accrued Interest`[combine[which(duplicated(c(search1,search2))==TRUE)]]<0){
          #Hier liegt eine kurze Periode vor
          if (settlementday < dates[1]){
            
            search3 <- which(bonddata$`Close of Business Date`==todaydate)
            search4 <- which(bonddata$ISIN==ISIN)
            combine <- c(search3,search4)
            
            if(bonddata$`Accrued Interest`[combine[which(duplicated(c(search3,search4))==TRUE)]]>0){
              ####PSEUDO CODE
              # SUCHE dem zugehorigen heutigen Stichtag die accrued interest
              # Wenn dieser Wert negativ oder gar 0 ist, dann kannst du abbrechen
              # Ist der Wert groesser als 0, dann muss du die neue Couponrate ausrechnen
              
              result$index <- "Short"
              result$first <- "Yes"
              nextdays <- as.numeric(difftime(dates[1],issuedate,units="days"))
              
              if (month(dates[1])!=divmonth2){
                regularperiod <- getcoupondate(divday2,divmonth2,year(dates[1])-1)
              } else {
                regularperiod <- getcoupondate(divday1,divmonth1,year(dates[1]))
              }
              
              fulldays <- as.numeric(difftime(dates[1],regularperiod,units="days"))
              result$adjustedcashflow <- (100*couponrate/2)*(nextdays/fulldays)
            } else {
              result <- list()
              result$index <- "Neither"
              result$Type <- "Between settlement day and coupon day"
            }
          } else {
            result <- list()
            result$index <- "Neither"
            result$Type <- "Valuation Date can't be in the first 2 periods"
          }
          
        } else {
          #Long first period
          
          #hier kann kein negativer ACC auftreten, da ich in der langen Periode bin 
          ## Kann ueber print getestet werden ## Test entfällt
          
          if (todaydate == dates[1]){
            # Bin ich genau am eigentlichen Coupontag, taucht das Datum nicht mehr in den CF in der Matrix auf
            result$first <- "No"
          } else {
            result$first <- "Yes"
          }
          result$index <- "Long"
          
          if (month(dates[1])!=divmonth2){
            regularperiod <- getcoupondate(divday2,divmonth2,year(dates[1])-1)
          } else {
            regularperiod <- getcoupondate(divday1,divmonth1,year(dates[1]))
          }
          
          fulldays <- as.numeric(difftime(dates[1],regularperiod,units="days"))
          nextdays <- as.numeric(difftime(dates[1],issuedate,units="days"))
          
          result$adjustedcashflow <- (100*couponrate/2)+((100*couponrate/2)*(nextdays/fulldays))
        }
        
        
      } else if ((todaydate > dates[1])&(settlementday < dates[2])) {
        #Need to check if it is a regular period or are we in a long period
        
        search1 <- which(bonddata$`Close of Business Date`==offset(dates[1],-4,"QuantLib/UnitedKingdom/Exchange"))
        search2 <- which(bonddata$ISIN==ISIN)
        combine <- c(search1,search2)
        
        if(!duplicated(c(search1,search2))[which(duplicated(c(search1,search2))==TRUE)]){
          print("2: Fehler in der Routine")
        } else if (length(duplicated(c(search1,search2))[which(duplicated(c(search1,search2))==TRUE)])!=1){
          print("2: Mehr als ein Treffer")
        }
        
        #Check if we are in regular period or long first period 
        #if (...) true, then we are in an long period
        if (bonddata$`Accrued Interest`[combine[which(duplicated(c(search1,search2))==TRUE)]]>0){
          # In 2. period after issue date
          # Only need adjustment for the long period, since the short is already handled
          
          search3 <- which(bonddata$`Close of Business Date`==todaydate)
          search4 <- which(bonddata$ISIN==ISIN)
          combine <- c(search3,search4)
          
          if(bonddata$`Accrued Interest`[combine[which(duplicated(c(search3,search4))==TRUE)]]>0){
            ####PSEUDO CODE
            # SUCHE dem zugehorigen heutigen Stichtag die accrued interest
            # Wenn dieser Wert negativ oder gar 0 ist, dann kannst du abbrechen
            # Ist der Wert groesser als 0, dann muss du die neue Couponrate ausrechen
            result <- list()
            result$index <- "Long"
            result$first <- "No"
            
            if (month(dates[1])!=divmonth2){
              regularperiod <- getcoupondate(divday2,divmonth2,year(dates[1])-1)
            } else {
              regularperiod <- getcoupondate(divday1,divmonth1,year(dates[1]))
            }
            
            fulldays <- as.numeric(difftime(dates[1],regularperiod,units="days"))
            nextdays <- as.numeric(difftime(dates[1],issuedate,units="days"))
            
            result$adjustedcashflow <- (100*couponrate/2)+((100*couponrate/2)*(nextdays/fulldays))
            
          } else {
            result <- list()
            result$index <- "Neither"
            result$Type <- "Between settlement day and coupon day"
          }
          
        } else {
          result <- list()
          result$index <- "Neither"
          result$Type <- "Normal period"
        }
        
      } else {
        result <- list()
        result$index <- "Neither"
        result$Type <- "Normal period"
      }
      
    } else {
      result <- list()
      result$index <- "Neither"
      result$Type <- "Valuation Date can't be in the first 2 periods"
    }
    
  } else {
    result <- list()
    result$index <- "Neither"
    result$Type <- "Issue Date is Coupon Date"
  }
  return(result)
}

### Loeschen von Anleihen, die zwischen Ex-Dividenden Tag und Maturity Day liegen

remove_outliers <- function(Maturitydate,accruedinterest,valuationdate,issuedate){
  if((accruedinterest<=0)&(newdaycount(valuationdate,Maturitydate,version = 0)<0.25)){
    checkvalue<- TRUE
  } else if (valuationdate<issuedate) {
    checkvalue<- TRUE
  }else {
    checkvalue<- FALSE 
  }
  return(checkvalue)
}

### Settlement Day Calculation
NextBuisnessday <- function(day1, delayday) {
  Settlementdate <- day1 + delayday
  wert1 <- bizdays(day1, Settlementdate, "QuantLib/UnitedKingdom/Exchange")
  i <- 1
  while (wert1 < delayday + 1) {
    Settlementdate <- day1 + delayday + i
    wert1 <- bizdays(day1, Settlementdate, "QuantLib/UnitedKingdom/Exchange")
    i <- i + 1
  }
  return(Settlementdate - 1)
}

### Coupondate 
coupondate <- function(day,month,buisnessday){
  year_buisnessday <- as.character(year(buisnessday))
  day <- as.character(day)
  month <- as.character(month)
  wert <- as.Date(paste(year_buisnessday,month,day, sep="-"))
  return(wert)
}

### Daycount Conventions
newdaycount <- function(date1, date2, version = 0){
  #date1 : Startdatum Format - as.Date
  #date2 : Enddatum Format - as.Date
  #Version : Num als verschiedene Möglichkeiten:
  #Option 0 (default) : Actual/365
  #Option 1 : 30/360 ISDA Source: 2006 ISDA Definitions (Section 4.16(f))
  #Option 2 : Actual/360 Source: 2006 ISDA Definitions (Section 4.16(e))
  #Option 3 : Actual/Actual ISDA Source:  2006 ISDA Definitions (Section 4.16(b)) 
  #Option 4 : Actual/Actual (ISMA) Source: 2006 ISDA Definitions (Section 4.16(c)) # Nicht Implementiert
  
  #Ausgabe Numerischer Wert 
  Y2 <- year(date2)
  Y1 <- year(date1)
  M2 <- month(date2)
  M1 <- month(date1)
  D2 <- day(date2)
  D1 <- day(date1)
  
  if (version == 0){
    #Option 0 
    year_diff <- as.numeric(difftime(date2,date1,units="days"))/365
  } else if (version == 1) {
    #Option 1
    D1 <- mapply(FUN = rules, x = D1)
    D2 <- mapply(FUN = rules1, x = D2, y = D1)
    
    year_diff <- ((Y2 - Y1)*360 + (M2 - M1)*30 + (D2 - D1)) / 360
  } else if (version == 2) { 
    #Option 2
    year_diff <- as.numeric(difftime(date2,date1,units="days"))/360
  } else if (version == 3) { 
    #Option 3
    year_diff <- mapply(FUN = option3fun, x = Y1, y = Y2, wert1 = date1, wert2 = date2)
    
  } else if (version == 4) { 
    #Option 4
    ### Not implemented yet ###
  }
  
  return(year_diff)
  
}

### Function estimating the frac between two coupon dates using 182.5 Days between coupon days

newdaycount182 <- function(date1, date2, Actual = 0){
  #date1 : Startdatum Format - as.Date
  #date2 : Enddatum Format - as.Date
  #Version : Num als verschiedene Möglichkeiten:
  #Option 0 (default) : Actual/182.5
  #Option 1 : 30/182.5 (following 2006 ISDA Definitions (Section 4.16(f)))
  
  #Ausgabe Numerischer Wert 
  Y2 <- year(date2)
  Y1 <- year(date1)
  M2 <- month(date2)
  M1 <- month(date1)
  D2 <- day(date2)
  D1 <- day(date1)
  
  if (Actual == 0){
    #Option 0 
    year_diff <- as.numeric(difftime(date2,date1,units="days"))/182.5
  } else if (Actual == 1) {
    #Option 1
    D1 <- mapply(FUN = rules, x = D1)
    D2 <- mapply(FUN = rules1, x = D2, y = D1)
    
    year_diff <- ((Y2 - Y1)*360 + (M2 - M1)*30 + (D2 - D1)) / 182.5
  } 
  return(year_diff)
  
}

last_day <- function(x) {
  lastyear <- as.numeric(format(x, "%Y"))
  last_date <- as.Date(sprintf('%s-12-31',lastyear))
}
first_day <- function(x) {
  thisyear <- as.numeric(format(x, "%Y"))
  last_date <- as.Date(sprintf('%s-01-01',thisyear))
}
rules <- function(x){
  if (x == 31) {
    x <- 30
  } else {
    x <- x
  }
  return(x)
}
rules1 <- function(x,y){
  if ((x == 31)&(y==30 | y==31)){
    x <- 30
  } else {
    x <- x
  }
  return(x)
}
getcoupondate <- function(x,y,z){
  last_date <- as.Date(sprintf('%s-%s-%s',z,y,x))
}
option3fun <- function(x,y,wert1,wert2){
  if (x == y){
    if (leap_year(wert2)) {
      year_diff <- as.numeric(difftime(wert2,wert1,units="days"))/366
    } else {
      year_diff <- as.numeric(difftime(wert2,wert1,units="days"))/365
    }
  } else {
    isleap <- leap_year(c(x:y))
    #numberleap366 <- length(isleap[isleap==TRUE])
    #numberleap365 <- length(leap_year(c(Y1:Y2))) - numberleap366
    a <- last_day(wert1)
    b <- first_day(wert2)
    last_year_days <- as.numeric(difftime(a,wert1,units="days"))+1
    first_year_days <- as.numeric(difftime(wert2,b,units="days"))
    if (leap_year(wert1)) {
      sum1 <- last_year_days/366
    } else {
      sum1 <- last_year_days/365
    }
    if (leap_year(wert2)) {
      sum2 <- first_year_days/366
    } else {
      sum2 <- first_year_days/365
    }
    year_diff <- sum1 + sum2 + (length(isleap)-2)
  }
  return(year_diff)
}

### Export Kalender
nonworkingdays <- function(calender1="UnitedKingdom/Exchange", from1=as.Date("2002-11-01"),to1 = as.Date("2022-09-01")){
  savingkalender <- getHolidayList(calendar = calender1, from = from1, to = to1, includeWeekends = FALSE)
  savingkalender1 <- rep(0, times=length(savingkalender)) 
  savingkalender2 <- rep(0, times=length(savingkalender)) 
  for (k in 1:length(savingkalender)){
    savingkalender1[k] <- format(savingkalender[k],"%d-%m-%Y")
    savingkalender2[k] <- 0
  }
  savingkalender <- cbind.data.frame(savingkalender1,savingkalender2)
  colnames(savingkalender) <- c("Date","Number of Obs.")
  
  write.xlsx(savingkalender,file ="KalenderCountry.xlsx",overwrite = TRUE,asTable = TRUE)
}

### Yield to Maturity nach TradeWeb

yieldtrabeweb <- function(cashflows, m, freque = 2,  searchint=c(-1,4), tol=1e-15, Annulized=TRUE){
  
  if (!is.matrix(cashflows))
    cashflows <- as.matrix(cashflows)
  if (!is.matrix(m))
    m <- as.matrix(m)
  
  bondyields <- matrix(0, nrow=ncol(cashflows), ncol=2)
  # put maturity of every bond into first column of bond yields matrix
  bondyields[,1] <- apply(m, 2, max)
  
  for (i in seq_len(ncol(cashflows))) {
    # present value of cash flows for root finding 
    pvcashflows2 <- function(y) {
      t(cashflows[,i])%*%(1/((1+y/freque)^m[,i]))
    }
    # calculate bond yields
    bondyields[i,2] <- uniroot(pvcashflows2, searchint, tol = tol,maxiter=3000)$root 
    if (Annulized){
      bondyields[i,2] <- ((1+bondyields[i,2]/2)^2)-1
    }
    
  }
  
  rownames(bondyields) <- colnames(cashflows)
  colnames(bondyields) <- c("Maturity","Yield")
  bondyields
  
  
}

maturities_matrix_1825 <- function(group,include_price=FALSE, Actual=0, SettlementCheck=FALSE) {
  
  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN),maxsum=100000) # Anzahl der CF pro Bond
  n_of_bonds <- length(n_of_cf) #Anzahl der ISIN Nummern / Anzahl der Bonds
  max_cf <- max(n_of_cf) #Max Anzahl an kommenden CF 
  pos_cf <- c(0,cumsum(n_of_cf))
  
  if (SettlementCheck) {
    year_diff <- newdaycount182(as.Date(group$`SETTLEMENT DAY`),as.Date(group$CASHFLOW$DATE),Actual)
  } else {
    year_diff <- newdaycount182(as.Date(group$TODAY),as.Date(group$CASHFLOW$DATE),Actual)
  }
  
  MATURITYMATRIX <-
    mapply(function(i) c(year_diff[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  MATURITYMATRIX <- matrixchanger(MATURITYMATRIX)
  
  
  if (include_price == TRUE) {MATURITYMATRIX <- rbind(rep(0,n_of_bonds),
                                                      MATURITYMATRIX)}  
  
  
  if (is.null(dim(MATURITYMATRIX))){
    MATURITYMATRIX <- matrix(MATURITYMATRIX, nrow = 1)
  } 
  colnames(MATURITYMATRIX) <- group$ISIN
  MATURITYMATRIX
}

matrixchanger <- function(Matrix1){
  for (i in 2:dim(Matrix1)[1]){
    for (j in 1:dim(Matrix1)[2]){
      if (Matrix1[[i,j]]!=0){
        Matrix1[[i,j]] <- Matrix1[[1,j]] + 1*(i-1)
      }
    }
  }
  return(Matrix1)
}

## Convert Excel to readable file for program
convertdata <- function(bonddata,NA_Ex_Dividend_Date){
  # Umwandlung in Numerische Werte
  bonddata$Coupon <- as.numeric(gsub(",", ".", bonddata$Coupon))
  bonddata$`Clean Price` <- as.numeric(gsub(",", ".", bonddata$`Clean Price`))
  bonddata$`Dirty Price` <- as.numeric(gsub(",", ".", bonddata$`Dirty Price`))
  bonddata$Yield <- as.numeric(gsub(",", ".", bonddata$Yield))
  bonddata$`Accrued Interest` <- as.numeric(gsub(",", ".", bonddata$`Accrued Interest`))
  bonddata$`Mod Duration` <- as.numeric(gsub(",", ".", bonddata$`Mod Duration`))
  
  # Umwandlung in Datum
  bonddata$`Close of Business Date` <- as.Date(paste(gsub(" ","",paste("20",substr(bonddata$`Close of Business Date`,7,8))), substr(bonddata$`Close of Business Date`,4,5), substr(bonddata$`Close of Business Date`,1,2), sep="-"))
  bonddata$Maturity <- as.Date(paste(gsub(" ","",paste("20",substr(bonddata$Maturity,7,8))), substr(bonddata$Maturity,4,5), substr(bonddata$Maturity,1,2), sep="-"))
  
  list <- rep(Sys.Date(), dim(bonddata)[1])
  under40 <- which(substr(bonddata$`Issue Date`,7,8)<40)
  
  list[under40] <- as.Date(paste(gsub(" ","",paste("20",substr(bonddata$`Issue Date`[under40],7,8))), substr(bonddata$`Issue Date`[under40],4,5), substr(bonddata$`Issue Date`[under40],1,2), sep="-"))
  list[-under40] <- as.Date(paste(gsub(" ","",paste("19",substr(bonddata$`Issue Date`[-under40],7,8))), substr(bonddata$`Issue Date`[-under40],4,5), substr(bonddata$`Issue Date`[-under40],1,2), sep="-"))
  bonddata$`Issue Date` <- list
  
  
  ### Behandeln der N/A in Accrued Interest
  # Einlesen der Ex-Dividend Days
  NA_Ex_Dividend_Date$`Close of Business Date` <- as.Date(paste(gsub(" ","",paste("20",substr(NA_Ex_Dividend_Date$`Close of Business Date`,7,8))), substr(NA_Ex_Dividend_Date$`Close of Business Date`,4,5), substr(NA_Ex_Dividend_Date$`Close of Business Date`,1,2), sep="-"))
  NA_Ex_Dividend_Date$`Closest Ex-Dividend Date` <- as.Date(paste(gsub(" ","",paste("20",substr(NA_Ex_Dividend_Date$`Closest Ex-Dividend Date`,7,8))), substr(NA_Ex_Dividend_Date$`Closest Ex-Dividend Date`,4,5), substr(NA_Ex_Dividend_Date$`Closest Ex-Dividend Date`,1,2), sep="-"))
  na_list <- which(is.na(bonddata$`Accrued Interest`))
  
  for (i in 1:length(which(is.na(bonddata$`Accrued Interest`)))){
    ISIN_search <- bonddata$ISIN[[na_list[i]]]
    CloseDate <- bonddata$`Close of Business Date`[[na_list[i]]]
    closest_Ex_Div_Date <- filter(NA_Ex_Dividend_Date, NA_Ex_Dividend_Date$ISIN==ISIN_search & NA_Ex_Dividend_Date$`Close of Business Date`==CloseDate)$`Closest Ex-Dividend Date`
    date1 <- coupondate(bonddata$`Dividend Date 1 Tag`[[na_list[i]]],bonddata$`Dividend Date 1 Monat`[[na_list[i]]],CloseDate)
    date2 <- coupondate(bonddata$`Dividend Date 2 Tag`[[na_list[i]]],bonddata$`Dividend Date 2 Monat`[[na_list[i]]],CloseDate)
    if (abs(bizdays(closest_Ex_Div_Date, date1, "QuantLib/UnitedKingdom/Exchange")) > abs(bizdays(closest_Ex_Div_Date, date2, "QuantLib/UnitedKingdom/Exchange"))){
      rev_coupondate <- date2
    } else {
      rev_coupondate <- date1
    }
    if (CloseDate > closest_Ex_Div_Date){
      if (rev_coupondate > CloseDate) {
        settlementday <- NextBuisnessday(CloseDate,delayday)
        if (settlementday >= rev_coupondate){
          bonddata$`Accrued Interest`[[na_list[i]]] <- 0
        } else {
          cat("Kein Tausch der N/A Werte. \n")
          print(i)
        }
      } else {
        cat("Rev_coupondate ist nicht größer als nächstes CloseDate.\n")
      }
      
    } else {
      cat("CloseDate ist nicht größer als nächstes ExDivDate.\n")
    }
  }
  
  
  ### Ordnen der Daten erst nach Close of Business Date und dann nach Maturity
  # Dies wird benötigt um die Daten in der CF Matrix sinnvoll eintragen zu koennen. Es muss dann naemlich keine nachfolgende Ordnung der Daten geschehen
  # So kann die Funktion addcashflows ohne Ordnung aufgerufen werden
  bonddata <- bonddata[with(bonddata, order(bonddata$`Close of Business Date`,bonddata$Maturity)),]
  
  return(bonddata)
}

## Create Empty List to fill
subset <- function(bonddata, searchdates=NULL){
  if(!is.null(searchdates)) {
    indexofdates <- which(bonddata$`Close of Business Date` %in% searchdates)
    con <- bonddata[indexofdates,]
    bonddata <- con
    return(bonddata)
  } else{
    return(bonddata)
  }
}

emptylist <- function(bonddata, country='Great Britian'){
  
  ### Erstellen der leeren Liste
  liste_leng <- length(bizseq(min(bonddata$`Close of Business Date`), max(bonddata$`Close of Business Date`), "QuantLib/UnitedKingdom/Exchange"))
  leere_Liste <- Listeerstelltung(liste_leng, country)
  for (i in 1:liste_leng){
    class(leere_Liste[[i]]) <- "couponbonds"}
  
  #Werktage Eintragen
  werktage <- bizseq(min(bonddata$`Close of Business Date`), max(bonddata$`Close of Business Date`), "QuantLib/UnitedKingdom/Exchange")
  
  for (i in 1:liste_leng){
    leere_Liste[[i]][[country]][["TODAY"]] <- werktage[i]
    ### Settlement Day einfuegen
    leere_Liste[[i]][[country]][["SETTLEMENT DAY"]] <- NextBuisnessday(leere_Liste[[i]][[country]][["TODAY"]],delayday)
  }
  return(leere_Liste)
}

## Add the data to the list 
adddata <- function(bonddata,country,leere_Liste,unfiltereddata){
  ### Hinzufuegen der Eintraege
  for (i in 1:dim(bonddata)[1]){
    if ((isBusinessDay("UnitedKingdom/Exchange", bonddata$`Close of Business Date`[[i]]))){
      if (!remove_outliers(bonddata$Maturity[[i]],bonddata$`Accrued Interest`[[i]],bonddata$`Close of Business Date`[[i]],bonddata$`Issue Date`[[i]])){
        leere_Liste <- add_bond(leere_Liste, country, bonddata$ISIN[[i]],bonddata$Maturity[[i]],bonddata$`Issue Date`[[i]],bonddata$Coupon[[i]],bonddata$`Clean Price`[[i]],bonddata$`Accrued Interest`[[i]],bonddata$`Close of Business Date`[[i]])
        #Adjustment of Cashflows between Ex-Dividend Day and Coupon Day is used
        settlementday <- NextBuisnessday(bonddata$`Close of Business Date`[[i]],delayday)
        result <- adjustmentrate(unfiltereddata,bonddata$ISIN[[i]],bonddata$`Close of Business Date`[[i]],bonddata$`Issue Date`[[i]],bonddata$Coupon[[i]],bonddata$`Dividend Date 1 Tag`[[i]],bonddata$`Dividend Date 1 Monat`[[i]],bonddata$`Dividend Date 2 Tag`[[i]],bonddata$`Dividend Date 2 Monat`[[i]],settlementday)
        leere_Liste <- add_cashflows(leere_Liste, country, bonddata$ISIN[[i]],bonddata$Maturity[[i]],bonddata$Coupon[[i]],bonddata$`Dividend Date 1 Tag`[[i]],bonddata$`Dividend Date 1 Monat`[[i]],bonddata$`Dividend Date 2 Tag`[[i]],bonddata$`Dividend Date 2 Monat`[[i]],bonddata$`Close of Business Date`[[i]],bonddata$`Accrued Interest`[[i]],result)
      } 
    }
  }
  DynData <- rm_bond(leere_Liste,country,"GB1111111111111")
  
  return(DynData)
}

## Check if the Input is correct
numobs <- function(bonddata, group){
  obs1 <- rep(as.Date("2022-08-05"), times=length(bonddata)) 
  obs2 <- rep(0, times=length(bonddata)) 
  for (k in 1:length(bonddata)) {
    obs1[k] <- as.Date(bonddata[[k]][[group]][["TODAY"]])
    obs2[k] <- length(bonddata[[k]][[group]][["ISIN"]])
  }
  obs <- cbind.data.frame(obs1,obs2)
  colnames(obs) <- c("Date","Number of Obs.")
  class(obs) <- "dyndatanumobs"
  return(obs)
}

readcheck.dyndatanumobs <- function(object1, object2){
  if (length(object1$Date) != length(object2$Date)){
    cat("Die Anzahl der eingelesenen Tage stimmt nicht mit der tatsaechlichen Anzahl der Tage ueberein.\n")
  } else {
    z <- (object1$Date == object2$Date)
    if (sum(z, na.rm = TRUE) != length(object1$Date)){
      cat("Es liegt ein Fehler in den Werktagen vor. \n")
      z1 <- which(z %in% FALSE)
      z2 <- object2$Date[z1]
      print(z2)
    } else {
      z <- (object1$`Number of Obs.` == object2$`Number of Obs.`)
      if (sum(z, na.rm = TRUE) != length(object1$Date)){
        cat("Es liegt ein Fehler in der Anzahl der Beobachtungen an folgenden Tagen vor: \n")
        z1 <- which(z %in% FALSE)
        z2 <- object2$Date[z1]
        print(z2)
      } else {
        cat("Die Daten wurden korrekt eingelesen.\n")
      }
    }
  }
}

changenumobs <- function(data1,searchdates){
  data1$Date <- as.Date(paste(gsub(" ","",paste("20",substr(data1$Date,7,8))), substr(data1$Date,4,5), substr(data1$Date,1,2), sep="-"))
  
  if (!is.null(searchdates)){
    indexofdates <- which(data1$Date %in% searchdates)
    con <- data1[indexofdates,]
    data1 <- con
  }
  class(data1) <- "dyndatanumobs"
  return(data1)
}

changenumobs2 <- function(bonddata, Number_of_Obs){
  ### Loeschen der Eintraege
  for (i in 1:dim(bonddata)[1]){
    if (remove_outliers(bonddata$Maturity[[i]],bonddata$`Accrued Interest`[[i]],bonddata$`Close of Business Date`[[i]],bonddata$`Issue Date`[[i]])){
      Number_of_Obs[["Number of Obs."]][which(Number_of_Obs[["Date"]] == bonddata$`Close of Business Date`[[i]])] <- Number_of_Obs[["Number of Obs."]][which(Number_of_Obs[["Date"]] == bonddata$`Close of Business Date`[[i]])] -1
    }
  }
  return(Number_of_Obs)
}


## Check, ob jeder Arbeitstag Eintraege hat:
readcheck.dyncouponbonds <- function(bondlist,country){
  valuestodelete <- c()
  date2 <- rep(as.Date("2022-08-05"), times=1) 
  i <- 0
  for (k in 1:length(bondlist)){
    if (length(bondlist[[k]][[country]]$ISIN) == 0){
      i <- i + 1
      valuestodelete[i] <- k
      date2[i] <- as.Date(bondlist[[k]][[country]]$TODAY)
    }
  }
  if (!is.null(valuestodelete)){
    bondlist <- bondlist[-valuestodelete]
  } else {
    date2 <- c()
  }
  class(bondlist) <- "dyncouponbonds"
  if (!is.null(date2)){
    cat("Folgende Daten enthielten keine Eintraege. \n")
    print(date2)
  } else {
    cat("Jeder Eintrag enthielt Daten. \n")
  }
  return(bondlist)
}

### Export the Single Values
exportquality <- function(object,bondlist,country){
  x <- object
  values <- matrix(0, nrow=length(x), ncol=8)
  
  werktage <- rep(as.Date("2022-08-05"), times=length(x))
  
  for (i in 1:length(x)){
    werktage[i] <- bondlist[[i]][[country]][["TODAY"]]
  }
  
  # werktage <- bizseq(bondlist[[1]][[country]][["TODAY"]], bondlist[[length(bondlist)]][[country]][["TODAY"]], "QuantLib/UnitedKingdom/Exchange")
  
  if(class(x) == "dyncouponbonds_cs"){
    for (i in 1:length(x)){
      for (j in 1:8){
        values[i,j] <- summary.termstrc_cs(x[[i]])$gof[j]
      }
    }
  } else {
    for (i in 1:length(x)){
      for (j in 1:8){
        values[i,j] <- summary.termstrc_nss(x[[i]])$gof[j]
      }
    }
  }
  colnames(values) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices","Weighted RMSE-Prices", "RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMSE-Yields (in %)")
  rownames(values) <- as.character(werktage)
  values <- as.data.frame(values)
  write.xlsx(values,file ="QualMeasures.xlsx",overwrite = TRUE,asTable = TRUE,rowNames = TRUE)
  class(values) <- "exportqual"
  return(values)
}

plot.exportqual <- function(x,ylim=c(),xlim=c(),lwd=2, type="l", ctype="RMSE",
                            xlab="Evaluation Dates",ylab="RMSE-Error", 
                            col="steelblue",lty=1, ...) 
{
  x1 <- seq(from = 1, to = length(x$`RMSE-Prices`))
  cdata <- switch(ctype, "RMSE" = x$`RMSE-Prices`,
                  "MAE" = x$`MAE-Prices`,
                  "SMAPE" = x$`SMAPE-Prices`,
                  "wRMSE" = x$`Weighted RMSE-Prices`)
  ylab <- switch(ctype, "RMSE" = "RMSE-Error",
                 "MAE" = "MAE-Error",
                 "SMAPE" = "SMAPE-Error",
                 "wRMSE" = "wRMSE-Error")
  
  plot(x1 ,cdata, type=type, ylim=ylim, xlim=xlim, xlab=xlab,
       ylab=ylab,lwd=lwd,lty=lty,col=col, ... )
  
}

# yieldbinder_ytm <- function(dataset,                  # dataset (static)
#                         group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
#                         matrange="all",            # maturity range in years c(min, max) 
#                         version = 0,              #  for Daycountversion
#                         SettlementCheck = FALSE,   #  for SettlementDay Usage
#                         freque=2,                  #  for Frequency of Coupon Payment
#                         mythology=TRUE,            #   for Different Yield Method 
#                         Annulized=TRUE,            #   for Yield Method 2 (Annulized/Semiannulized)
#                         Actual = 0,...){
#   
#   ### 25 ist erstmal nur für die ersten Datensaetze
#   yields <- matrix(0, nrow = length(dataset), ncol = 41)
#   
#   for (i in 1:length(dataset)){
#     prepro <- prepro_bond(group=group,bonddata=dataset[[i]],matrange=matrange, version=version, SettlementCheck=SettlementCheck,freque=freque, mythology=mythology, Annulized=Annulized, Actual = Actual)
#     yields[i,] <- t(as.matrix(unname(prepro$y[[1]][,2])))
#   }
#   return(yields)
# }

## Creating and n x m Cashflow Matrix used for couponbonds
## n... Number of Bonds 
## m .. Number of unique Casflows
create_cashflow_matrix_2_0 <- function(group){
  cashflows <- matrix(0, nrow = length(group[["ISIN"]]), ncol = length(unique(group$CASHFLOWS$DATE)))
  for (i in 1: length(unique(group$CASHFLOWS$DATE))){
    index <- which(sort(unique(group$CASHFLOWS$DATE))[i] == group$CASHFLOWS$DATE)
    for (j in 1:length(index)){
      zeilenindex <- which(group$CASHFLOWS$ISIN[index[j]] == group$ISIN)
      cashflows[zeilenindex,i] <- group$CASHFLOWS$CF[index[j]]
    }
  }
  return(cashflows)
}


#########################################################
### Cubic splines estimation method for 'couponbonds' ###
### General Cubic B-splines for Discont function      ###
### Following the book of Filipovic                   ###
#########################################################

estim_cs <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",setdf=11,freque=2) UseMethod("estim_cs")
###   for Daycountversion

### Cubic spline term structure estimation 
estim_cs.couponbonds <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",setdf=11,freque=2) {
  
  ## data preprocessing 
  ###   for Daycountversion
  prepro <- prepro_bond(group=group,bonddata=bonddata,
                        matrange=matrange,version=version, SettlementCheck=SettlementCheck,freque=freque, mythology=mythology, Annulized=Annulized, Actual = Actual, wmethod=wmethod)
  
  
  n_group=prepro$n_group
  sgroup=prepro$sgroup
  cf=prepro$cf
  cf_p=prepro$cf_p
  m=prepro$m
  m_p=prepro$m_p
  p=prepro$p
  ac=prepro$ac
  y=prepro$y
  duration=prepro$duration
  m_p_2=prepro$m_p_2
  
  if (setdf==0){
    setdf = min(mapply(function(k) length(bonddata[[k]][["ISIN"]]),sgroup))
  }
  
  bonddata1 <- bonddata[group]
  cashflows <- lapply(bonddata1, create_cashflow_matrix_2_0)
  x <- mapply(function(k) newdaycount(bonddata[[k]]$TODAY,sort(unique(bonddata[[k]]$CASHFLOWS$DATE)),version = version),sgroup,SIMPLIFY=FALSE)
  # x <- lapply(newdaycount(bonddata$group$TODAY,sort(unique(group$CASHFLOWS$DATE)),version = version))
  
  price <- mapply(function(k) unname(p[[k]]),sgroup,SIMPLIFY=FALSE)
  regout <- mapply(function(k) lm(price[[k]]~cashflows[[k]]%*%bs(x[[k]],df =setdf, degree = 3,intercept = TRUE)-1,weights = unname(duration[[k]][,3])),sgroup,SIMPLIFY=FALSE)
  
  # estimated paramters  
  alpha <- lapply(regout, coef)
  
  # calculate estimated prices 
  phat <- mapply(function(k) cashflows[[k]]%*%(predict(bs(x[[k]],df =setdf, degree = 3,intercept = TRUE), x[[k]])%*%(unname(alpha[[k]]))),sgroup,SIMPLIFY=FALSE)
  # phat <- lapply(regout, fitted.values)
  
  # price errors
  # perrors <- lapply(regout, residuals)    
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(t(-phat[[k]]),cf[[k]]),m_p[[k]],m_p_2[[k]], freque=freque, mythology=mythology, Annulized=Annulized),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"
  
  # maturity interval
  
  #alt 
  t <- mapply(function(k) x[[k]], sgroup,SIMPLIFY=FALSE)
  
  #t <- mapply(function(k) seq(round(attr(bs(x[[k]],df =setdf, degree = 3,intercept = TRUE),"Boundary.knots")[1]+0.02,2),round(attr(bs(x[[k]],df =setdf, degree = 3,intercept = TRUE),"Boundary.knots")[2]-0.05,2),0.01), sgroup,SIMPLIFY=FALSE) 
  
  # calculate mean and variance of the distribution of the discount function 
  # test <- predict(bs(t[[k]],df =setdf, degree = 3,intercept = TRUE), newX)%*%(unname(alpha[[k]]))
  #mean_d <- mapply(function(k) (predict(bs(x[[k]],df =setdf, degree = 3,intercept = TRUE), t[[k]])%*%(unname(alpha[[k]]))), sgroup,SIMPLIFY=FALSE)
  
  #Alt
  mean_d <- mapply(function(k) bs(t[[k]],df =setdf, degree = 3,intercept = TRUE)%*%(unname(alpha[[k]])), sgroup,SIMPLIFY=FALSE)
  
  #  variance covariance matrix for estimated parameters
  if(rse) Sigma <- lapply(regout,vcovHAC.default) else Sigma <- lapply(regout,vcov) 
  
  var_d <- mapply(function(k) apply(bs(t[[k]],df =setdf, degree = 3,intercept = TRUE),1,function(x) t(x)%*%Sigma[[k]]%*%x), sgroup, SIMPLIFY=FALSE) 
  
  # lower 95% confidence interval
  ## 11 df 
  cl <- mapply(function(k) mean_d[[k]] + rep(qt(0.025,setdf),
                                             length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE)	
  for (j in 1:n_group){
    laenge <- length(which((cl[[j]]<=0)==TRUE))
    lliste <- which((cl[[j]]<=0)==TRUE)
    if (laenge != 0){
      for (i in 1:laenge){
        cl[[j]][lliste[i]] <- 0.0001
      }
    }
  }
  
  # upper 95 % confidence interval	
  cu <- mapply(function(k) mean_d[[k]] + rep(qt(0.975,setdf),
                                             length(mean_d[[k]]))*sqrt(var_d[[k]]),sgroup,SIMPLIFY=FALSE) 
  for (j in 1:n_group){
    laenge <- length(which((cu[[j]]<=0)==TRUE))
    lliste <- which((cu[[j]]<=0)==TRUE)
    if (laenge != 0){
      for (i in 1:laenge){
        cu[[j]][lliste[i]] <- 0.0001
      }
    }
  }
  
  # zero cupon yield curves for maturity interval t 
  zcy_curves <-  mapply(function(k)  cbind(t[[k]],-log(mean_d[[k]])/t[[k]],-log(cl[[k]])/t[[k]],-log(cu[[k]])/t[[k]]),sgroup, SIMPLIFY=FALSE)  
  
  
  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
  
  
  #calculate spread curves
  if(n_group != 1) {
    srange <- seq(max(unlist(lapply(t,min))),min(unlist(lapply(t,max))),0.01)
    s_curves <- mapply(function(k) cbind(srange,zcy_curves[[k]][c(which(zcy_curves[[k]][,1]== srange[1]): which(zcy_curves[[k]][,1] == srange[length(srange)])),2]
                                         - zcy_curves[[1]][c(which(zcy_curves[[1]][,1]== srange[1]): which(zcy_curves[[1]][,1]== srange[length(srange)])),2]),sgroup, SIMPLIFY=FALSE)
    
  } else s_curves = "none"
  for (k in sgroup) class(s_curves[[k]]) <- "ir_curve"
  class(s_curves) <- "s_curves"
  
  # create discount factor curves 
  df_curves <- mapply(function(k) cbind(t[[k]],mean_d[[k]]),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
  class(df_curves) <- "df_curves"	
  
  # calculate forward rate curves
  
  fwr_curves <- mapply(function(k) cbind(t[[k]],impl_fwr(zcy_curves[[k]][,1],zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"	
  
  knotpoints <- mapply(function(k) c(attr(bs(t[[k]],df =setdf, degree = 3,intercept = TRUE),"Boundary.knots"),attr(bs(t[[k]],df =setdf, degree = 3,intercept = TRUE),"knots")),sgroup,SIMPLIFY=FALSE)
  
  
  # return list of results
  result <- list(  group=group,          # e.g. countries, rating classes
                   matrange=matrange,    # maturity range of bonds
                   n_group=n_group,      # number of groups
                   knotpoints=knotpoints, # knot points
                   rse=rse,              # robust standard errors
                   spot=zcy_curves, 	# zero coupon yield curves
                   spread=s_curves,      # spread curves
                   discount=df_curves,	# forward rate curves
                   forward=fwr_curves,	# discount factor curves
                   cf=cf,                # cashflow matrix
                   m=m,                  # maturity matrix
                   p=p,                  # dirty prices
                   phat=phat,            # estimated prices
                   perrors=perrors,     	# price errors
                   y=y,                  # maturities and yields
                   yhat=yhat,            # estimated yields
                   yerrors=yerrors,	# yield errors
                   alpha=alpha,          # cubic splines parameters                             
                   regout=regout,         # OLS output
                   duration = duration,
                   x = x                # LOOC
                   
  )
  # assign names to results list 
  for ( i in 6:length(result)) names(result[[i]]) <- group 
  
  class(result) <- "termstrc_cs"
  result
  
}

############################################################
### Cubic splines estimation method for 'dyncouponbonds' ###
### Summary and Print summary                            ###
############################################################

estim_cs.dyncouponbonds <- function(bonddata, group, matrange="all",rse=FALSE, version = 0, SettlementCheck=FALSE, mythology = TRUE, Annulized=TRUE, Actual = 0,wmethod="rf",setdf=11,freque=2){
  solution_list <- vector(mode = "list", length = length(bonddata))
  for (i in 1:length(bonddata)){
    # print(i)
    cs_res <- estim_cs(bonddata[[i]], group=group, matrange=matrange,rse=rse, version = version, SettlementCheck=SettlementCheck, mythology = mythology, Annulized=Annulized, Actual = Actual,wmethod=wmethod,setdf=setdf,freque=freque)
    solution_list[[i]] <- cs_res
  }
  class(solution_list) <- "dyncouponbonds_cs"
  return(solution_list)
}

summary.dyncouponbonds_cs <- function(object, ...) {
  x <- object
  
  sumry <- list()
  
  f_p_mrsme <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_maabse <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_mrsme <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_maabse <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_smape <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_smape <- vector(mode = "list", length = length(x[[1]]$group))
  f_p_wrmse <- vector(mode = "list", length = length(x[[1]]$group))
  f_y_wrmse <- vector(mode = "list", length = length(x[[1]]$group))
  
  for (z in 1:length(x)){
    subset <- x[[z]]
    RMSE_p <- mapply(function(i) rmse(subset$p[[i]],subset$phat[[i]]),seq(subset$n_group))
    AABSE_p <- mapply(function(i) aabse(subset$p[[i]],subset$phat[[i]]) ,seq(subset$n_group))
    RMSE_y <- mapply(function(i) rmse(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100) ,seq(subset$n_group))
    AABSE_y <- mapply(function(i) aabse(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100) ,seq(subset$n_group))
    
    sMAPE_p <- mapply(function(i) smape(subset$p[[i]],subset$phat[[i]]),seq(subset$n_group))
    sMAPE_y <- mapply(function(i) smape(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100),seq(subset$n_group))
    
    wRMSE_p <- mapply(function(i) rmsewe(subset$p[[i]],subset$phat[[i]],subset[["duration"]][[i]][,3]),seq(subset$n_group))
    wRMSE_y <- mapply(function(i) rmsewe(subset$y[[i]][,2]*100,subset$yhat[[i]][,2]*100,subset[["duration"]][[i]][,3]),seq(subset$n_group))
    
    
    for (u in 1:length(x[[1]]$group)){
      f_p_mrsme[[u]] <- c(f_p_mrsme[[u]],RMSE_p[[u]])
      f_p_maabse[[u]] <- c(f_p_maabse[[u]],AABSE_p[[u]])
      f_y_mrsme[[u]] <- c(f_y_mrsme[[u]],RMSE_y[[u]])
      f_y_maabse[[u]] <- c(f_y_maabse[[u]],AABSE_y[[u]])
      f_p_smape[[u]] <- c(f_p_smape[[u]],sMAPE_p[[u]])
      f_y_smape[[u]]<- c(f_y_smape[[u]],sMAPE_y[[u]])
      f_p_wrmse[[u]] <- c(f_p_wrmse[[u]],wRMSE_p[[u]])
      f_y_wrmse[[u]] <- c(f_y_wrmse[[u]],wRMSE_y[[u]])
      
    }
  }
  
  p_mrsme <- mapply(function(k) mean(f_p_mrsme[[k]]), 1:length(x[[1]]$group))
  p_maabse <- mapply(function(k) mean(f_p_maabse[[k]]), 1:length(x[[1]]$group))
  p_smape <- mapply(function(k) mean(f_y_mrsme[[k]]), 1:length(x[[1]]$group))
  y_mrsme <- mapply(function(k) mean(f_y_maabse[[k]]), 1:length(x[[1]]$group))
  y_maabse <- mapply(function(k) mean(f_p_smape[[k]]), 1:length(x[[1]]$group))
  y_smape <- mapply(function(k) mean(f_y_smape[[k]]), 1:length(x[[1]]$group))
  p_wRMSE <- mapply(function(k) mean(f_p_wrmse[[k]]), 1:length(x[[1]]$group))
  y_wRMSE <- mapply(function(k) mean(f_y_wrmse[[k]]), 1:length(x[[1]]$group))
  
  sumry$gof <- rbind(p_mrsme,p_maabse,p_smape,p_wRMSE,y_mrsme,y_maabse,y_smape,y_wRMSE)
  colnames(sumry$gof) <- x[[1]]$group
  rownames(sumry$gof) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices","Weighted RMSE-Prices" ,"RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMSE-Yields (in %)")
  
  class(sumry) <- "summary.dyntermstrc_cs"
  sumry
}

print.summary.dyntermstrc_cs <- function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Goodness of fit:\n")
  cat("---------------------------------------------------\n")
  
  print.default(format(x$gof,digits=6,scientific=FALSE),quote=FALSE)
}

#########################################################
### Extracting Number of Bonds in specific          #####
### Maturity Brackets                               #####
#########################################################

restlaufzeitklassen <- function(bonddata,group,include_price=FALSE,version = 3,SettlementCheck=FALSE,freque=2,mythology=TRUE,wmethod="rf",Actual=0,Annulized=TRUE){
  # number of groups 
  n_group <- seq(length(bonddata))
  
  maturities <- mapply(function(j) create_maturities_matrix(bonddata[[j]][[group]],include_price=include_price,version = version,SettlementCheck=SettlementCheck),n_group,SIMPLIFY=FALSE)
  colmax <- function(m) apply(m,2,max)
  maturities_bonds <- mapply(function(j) colmax(maturities[[j]]),n_group,SIMPLIFY=FALSE)
  
  m_p <- mapply(function(j) create_maturities_matrix(bonddata[[j]][[group]],include_price=TRUE,version = version,SettlementCheck=SettlementCheck),n_group,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]][[group]],include_price=TRUE),n_group,SIMPLIFY=FALSE)
  m_p_2 <- mapply(function(k) maturities_matrix_1825(bonddata[[k]][[group]],include_price=TRUE,Actual = Actual,SettlementCheck=SettlementCheck),n_group,SIMPLIFY=FALSE)
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]],m_p_2[[k]], freque=freque, mythology=mythology, Annulized=Annulized),n_group,SIMPLIFY=FALSE)
  
  duration <- mapply(function(k) duration(cf_p[[k]],m_p[[k]],y[[k]][,2],m_p_2[[k]],freque,mythology,wmethod),n_group,SIMPLIFY=FALSE)
  
  
  restlaufzeit <- matrix(0,nrow = length(bonddata),ncol = 17)
  for (i in 1:dim(restlaufzeit)[1]){
    restlaufzeit[i,1] <- length(unname(maturities_bonds[[i]]))
    restlaufzeit[i,2] <- sum(unname(maturities_bonds[[i]])<1, na.rm=TRUE)
    for (j in 3:11){
      restlaufzeit[i,j] <- sum(unname(maturities_bonds[[i]])>=j-1 & unname(maturities_bonds[[i]])<j, na.rm=TRUE)
    }
    for (j in 1:4){
      restlaufzeit[i,j+11] <- sum(unname(maturities_bonds[[i]])>=(50-(5-j)*10) & unname(maturities_bonds[[i]])<(50-(4-j)*10), na.rm=TRUE)
    }
    restlaufzeit[i,16] <- sum(unname(maturities_bonds[[i]])>50, na.rm=TRUE)
    restlaufzeit[i,17] <- mean(unname(duration[[i]][,2]), na.rm=TRUE)
  }
  
  
  
  colnames(restlaufzeit) <- c("Gesamt","RLZ < 1","1<= RLZ < 2","2<= RLZ < 3","3<= RLZ < 4","4<= RLZ < 5","5<= RLZ < 6","6<= RLZ < 7","7<= RLZ < 8","8<= RLZ < 9","9<= RLZ < 10","10<= RLZ < 20","20<= RLZ < 30","30<= RLZ < 40","40<= RLZ < 50","RLZ > 50","Average Tenor")
  werktage <- bizseq(bonddata[[1]][[group]][["TODAY"]], bonddata[[length(bonddata)]][[group]][["TODAY"]], "QuantLib/UnitedKingdom/Exchange")
  rownames(restlaufzeit) <- as.character(werktage)
  
  restlaufzeit<- as.data.frame(restlaufzeit)
  write.xlsx(restlaufzeit,file ="Restlaufzeitklassen.xlsx",overwrite = TRUE,asTable = TRUE,rowNames = TRUE)
  
  return(restlaufzeit)
  
}

#########################################################
### Extracting yields to predefinded maturity dates #####
#########################################################

extractingyields <- function(evalmodel,country,range=c(0.25,0.5,1,2,3,5,7,10,15,30)){
  yieldsforpca <- matrix(0, nrow=length(evalmodel),ncol = length(range))
  for (i in 1:length(evalmodel)){
    for (j in 1:length(range)){
      yieldsforpca[i,j] <- evalmodel[[i]][["spot"]][[country]][which(round(evalmodel[[i]][["spot"]][[country]][,1],2) %in% round(range[j],2)),2]
    }
  }
  
  namvector <- character(length(range))
  for (i in 1:length(range)){
    if (range[i]<1){
      namvector[i] <- paste(as.character(range[i]*12),"M",sep = "")
    } else {
      namvector[i] <- paste(as.character(range[i]),"Y",sep = "")
    }
  }
  colnames(yieldsforpca) <- namvector
  
  return(yieldsforpca)
}

#########################################################
### Prinicpal component analysis                    #####
#########################################################

princomp <- function(pcayields){
  
  diffs = diff(pcayields, 1)
  result_pca = prcomp(diffs, scale = TRUE, center = TRUE)
  
  return(result_pca)
}

pve <- function(resultpca,cumulative=TRUE,limitpca = 0.98){
  # Compute variance
  my_pca.var <- resultpca$sdev ^ 2
  # Proportion of variance for a scree plot
  propve <- my_pca.var / sum(my_pca.var)
  if (!cumulative){
    plot(propve, xlab = "principal component",
         ylab = "Proportion of Variance Explained",
         ylim = c(0, 1), type = "b",
         main = "Scree Plot")
  } else {
    plot(cumsum(propve),
         xlab = "Principal Component",
         ylab = "Cumulative Proportion of Variance Explained",
         ylim = c(0, 1), type = "b")
  }
  pcaan <- which(cumsum(propve) >= limitpca)[1]
  cat("Number of Prinicpal Components needed to reach the threshold of Limitpca:\n")
  print(pcaan)
  cat("---------------------------------------------------\n")
  cat("Cumulative Proportion of Variance Explained \n")
  print(cumsum(resultpca$sdev^2)/sum(resultpca$sdev^2))
}

plot.prcomp <- function(resultpca,pcanumber="3",Flipp=TRUE,range=c(0.25,0.5,1,2,3,5,7,10,15,30),...) {
  if (Flipp){ # first principal component sign flip so it's an upward shift in yield
    resultpca$x[,1] <- -resultpca$x[,1]
    # also have to flip sign for eigenvector
    resultpca$rotation[,1] <- -resultpca$rotation[,1]}
  
  x1 <- seq(from = 1, to = length(resultpca$rotation[,1]))
  #x1<- range
  cdata <- mapply(function(i) resultpca$rotation[,i],1:as.numeric(pcanumber) ,SIMPLIFY = FALSE)
  
  minvalue <- round(min(mapply(function(i) min(cdata[[i]]), 1:as.numeric(pcanumber))),2)-0.5
  maxvalue <- ceiling(max(mapply(function(i) max(cdata[[i]]), 1:as.numeric(pcanumber))))+0.2
  
  cdata2 <- as.data.frame(cdata)
  write.xlsx(cdata2,file ="cdata.xlsx",overwrite = TRUE,asTable = TRUE,rowNames = TRUE)
  
  clegend <- switch(pcanumber, "1" = c("1. Loading"),
                    "2" = c("1. Loading","2. Loading"),
                    "3" = c("1. Loading","2. Loading","3. Loading"),
                    "4" = c("1. Loading","2. Loading","3. Loading","4. Loading"),
                    "5" = c("1. Loading","2. Loading","3. Loading","4. Loading","5. Loading"))
  ccol <- switch(pcanumber, "1" = c("red"),
                 "2" = c("red","steelblue"),
                 "3" = c("red","steelblue","green"),
                 "4" = c("red","steelblue","green","coral"),
                 "5" = c("red","steelblue","green","coral","deeppink"))
  
  plot(x1 ,unname(cdata[[1]]), type="l",xlab="Range",ylab="Loading", 
       col="red",lty=1,ylim = c(minvalue, maxvalue),bty='L')
  for (i in 2:as.numeric(pcanumber)){
    lines(x1,unname(cdata[[i]]), col=ccol[i])
  }
  legend("bottom",legend=clegend,col=ccol,lty = 1:1,cex=0.4)
}

datarebuildpca <- function(resultpca,yieldsforpca1,pcanumber=3){
  # data re-built using just first two principal components
  ypc <- resultpca$x[,1:pcanumber] %*% t(resultpca$rotation[,1:pcanumber])
  # scale the re-built series to have same mean and stdev as original series y
  ypc <- t(apply(ypc,1,function(x) x * resultpca$scale + resultpca$center))
  #Since we fit the difference, we need to adjust back
  newyieldsforpca <- matrix(data=0,nrow = dim(yieldsforpca1)[1],ncol = dim(yieldsforpca1)[2])
  newyieldsforpca[1,] <- unname(yieldsforpca1[1,])
  colnames(newyieldsforpca) <- colnames(yieldsforpca1)
  for (i in 2:dim(yieldsforpca1)[1]){
    for (j in 1:dim(yieldsforpca1)[2]){
      newyieldsforpca[i,j] <- newyieldsforpca[i-1,j]+ypc[i-1,j]
    }
  }
  return(newyieldsforpca)
}

plotyields <- function(newyieldsforpca1,yieldsforpca1,Multiple=TRUE,cplot="PCA",range=c(0.25,0.5,1,2,3,5,7,10,15,30)){
  
  if(Multiple&(length(range)>6)){
    cat("---------------------------------------------------\n")
    cat("If Multiple is TRUE range should only contain a maximum of 5 elements \n")
    cat("---------------------------------------------------\n")
    cat("Options to chose from: \n")
    print(colnames(newyieldsforpca1)[which(colnames(newyieldsforpca1)==colnames(yieldsforpca1))])
  } else {
    
    namvector <- character(length(range))
    for (i in 1:length(range)){
      if (range[i]<1){
        namvector[i] <- paste(as.character(range[i]*12),"M",sep = "")
      } else {
        namvector[i] <- paste(as.character(range[i]),"Y",sep = "")
      }
    }
    
    x1 <- seq(from = 1, to = dim(yieldsforpca1)[1])
    cdata <- switch(cplot, "PCA" = newyieldsforpca1,
                    "Model" = yieldsforpca1)
    valuerange = as.character(length(range))
    
    ccol1 <- switch(valuerange, "1" = c("brown1"),
                    "2"  = c("brown1","bisque2"),
                    "3"  = c("brown1","bisque2","darkseagreen2"),
                    "4"  = c("brown1","bisque2","darkseagreen2","coral"),
                    "5"  = c("brown1","bisque2","darkseagreen2","coral","deeppink"),
                    "6"  = c("brown1","bisque2","darkseagreen2","coral","deeppink","deepskyblue"),
                    "7"  = c("brown1","bisque2","darkseagreen2","coral","deeppink","deepskyblue","darkgoldenrod1"),
                    "8"  = c("brown1","bisque2","darkseagreen2","coral","deeppink","deepskyblue","darkgoldenrod1","cyan"),
                    "9"  = c("brown1","bisque2","darkseagreen2","coral","deeppink","deepskyblue","darkgoldenrod1","cyan","cornsilk"),
                    "10" = c("brown1","bisque2","darkseagreen2","coral","deeppink","deepskyblue","darkgoldenrod1","cyan","cornsilk","burlywood1"))
    
    ccol2 <- switch(valuerange, "1" = c("brown4"),
                    "2" = c("brown4","antiquewhite4"),
                    "3" = c("brown4","antiquewhite4","darkseagreen4"),
                    "4" = c("brown4","antiquewhite4","darkseagreen4","coral3"),
                    "5" = c("brown4","antiquewhite4","darkseagreen4","coral3","deeppink4"))
    
    ccol <- c(ccol1,ccol2)
    
    ylimit = c((min(min(newyieldsforpca1),min(yieldsforpca1)))*100-0.5,(max(max(newyieldsforpca1),max(yieldsforpca1)))*100+0.25)
    
    
    if(Multiple){
      plot(x1 ,unname(newyieldsforpca1[,which(colnames(newyieldsforpca1) %in% namvector)[1]])*100, type="l",xlab="Dates",ylab="Yield", 
           col=ccol1[1],lty=1,ylim = ylimit)
      for (i in 2:length(which(colnames(newyieldsforpca1) %in% namvector))) {
        lines(x1,unname(newyieldsforpca1[,which(colnames(newyieldsforpca1) %in% namvector)[i]])*100, col=ccol1[i])
      }
      for (i in 1:length(which(colnames(yieldsforpca1) %in% namvector))) {
        lines(x1,unname(yieldsforpca1[,which(colnames(yieldsforpca1) %in% namvector)[i]])*100, col=ccol2[i])
      }
      namvector <- c(paste("PCA",namvector),paste("Model",namvector))
      legend("bottom",legend=namvector,col=ccol,lty = 1:1,cex=0.4,ncol = 2)
      title("Comparison between PCA and Model")
      
    } else {
      plot(x1 ,unname(cdata[,which(colnames(cdata) %in% namvector)[1]])*100, type="l",xlab="Dates",ylab="Yield", 
           col=ccol1[1],lty=1,ylim = ylimit)
      for (i in 2:length(which(colnames(cdata) %in% namvector))){
        lines(x1,unname(cdata[,which(colnames(cdata) %in% namvector)[i]])*100, col=ccol1[i])
      }
      legend("bottom",legend=namvector,col=ccol,lty = 1:1,cex=0.4)
      title(cplot)
    }
  }
}

#########################################################
### LOOC                                            #####
#########################################################

# group <- c("Great Britian")
# matrange="all"
# method="ns"
# version = 0             #   for Daycountversion
# rse = FALSE
# SettlementCheck = FALSE   #   for SettlementDay Usage
# freque=2                  #   for Frequency of Coupon Payment
# mythology=TRUE            #   for Different Yield Method
# Annulized=TRUE            #   for Yield Method 2 (Annulized/Semiannulized)
# Actual = 0                #   for Daycount convention Yield Method 2
# wmethod="rf"              #   for Weight Method
# startparam=NULL           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
# # otherwise globally optimal parameters are searched automatically
# lambda=0.0609*12          # yearly lambda-value for "Diebold/Li" estimation
# tauconstr = NULL         # constraints for tau parameter grid
# constrOptimOptions = list(control = list(maxit = 2000), outer.iterations = 200, outer.eps = 1e-07)
# setdf = 11

looc<- function(dataset,Model,group,matrange="all",method="ns",version = 0,rse = FALSE, SettlementCheck = FALSE,freque=2,mythology=TRUE,Annulized=TRUE,Actual = 0,wmethod="rf",startparam=NULL,lambda=0.0609*12,tauconstr = NULL,setdf= 11,constrOptimOptions = list(control = list(maxit = 2000), outer.iterations = 200, outer.eps = 1e-03),...)UseMethod("looc") 

looc.dyncouponbonds <- function(dataset,
                                Model,                  # dataset (Dynamic)
                                group,                     # names of countries for estimation c("Country 1", "Country 2", ...)
                                matrange="all",            # maturity range in years c(min, max)
                                method="ns",
                                version = 0,               #   for Daycountversion
                                rse = FALSE, 
                                SettlementCheck = FALSE,   #   for SettlementDay Usage
                                freque=2,                  #   for Frequency of Coupon Payment
                                mythology=TRUE,            #   for Different Yield Method
                                Annulized=TRUE,            #   for Yield Method 2 (Annulized/Semiannulized)
                                Actual = 0,                #   for Daycount convention Yield Method 2
                                wmethod="rf",              #   for Weight Method
                                startparam=NULL,           # startparameter matrix with columns c("beta0","beta1","beta2","tau1","beta3","tau2")
                                # otherwise globally optimal parameters are searched automatically
                                lambda=0.0609*12,          # yearly lambda-value for "Diebold/Li" estimation
                                tauconstr = NULL,          # constraints for tau parameter grid
                                setdf= 11,
                                constrOptimOptions = list(control = list(maxit = 1000), outer.iterations = 200, outer.eps = 1e-03),...
){
  
  if (length(dataset)<20){
    cat("---------------------------------------------------\n")
    cat("The Dynamic Set needs to include at least 5 Buisnessdays \n")
    cat("---------------------------------------------------\n")
  } else {
    ## Start with a random selection of the Testgroup
    samplesize <- sort(sample(1:length(dataset), floor(length(dataset)/20), replace=FALSE))
    subset.dataset <- dataset[samplesize]
    
    maxvalueISIN <- max(mapply(function(k) length(subset.dataset[[k]][[group]]$ISIN), 1:length(subset.dataset)))
    
    phat <- matrix(0,nrow = length(subset.dataset), ncol = maxvalueISIN-2)
    p <- matrix(0,nrow = length(subset.dataset), ncol = maxvalueISIN-2)
    y <- matrix(0,nrow = length(subset.dataset), ncol = maxvalueISIN-2)
    yhat <- matrix(0,nrow = length(subset.dataset), ncol = maxvalueISIN-2)
    dura <- matrix(0,nrow = length(subset.dataset), ncol = maxvalueISIN-2)
    
    if (Model == "McCulloch"){
      
      for (i in 1:length(subset.dataset)){
        subsubset.dataset <- subset.dataset[[i]]
        for (j in 2:(length(subsubset.dataset[[group]]$ISIN)-1)){
          loocset <- rm_bond.couponbonds(subsubset.dataset, group, subsubset.dataset[[group]]$ISIN[[j]])
          mc_result <- estim_cs_mc(loocset, group, matrange=matrange,rse=rse, version = version, SettlementCheck=SettlementCheck, mythology = mythology, Annulized=Annulized, Actual = Actual,wmethod=wmethod)
          
          cf <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=FALSE)[,j]
          m <- create_maturities_matrix(subsubset.dataset[[group]],include_price=FALSE, version=version, SettlementCheck=SettlementCheck)[,j]
          p[i,j-1] <- subsubset.dataset[[group]][["PRICE"]][j] + subsubset.dataset[[group]][["ACCRUED"]][j]
          
          dt <- rep(1,length(m))
          for(sidx in 1:mc_result[["s"]][[group]]){
            dt <- dt + mc_result[["alpha"]][[group]][sidx]* mapply(function(p) gi(m[p],mc_result[["knotpoints"]][[1]],sidx,mc_result[["s"]][[group]]),1:length(m))
          }
          
          
          cf_p <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=TRUE)
          m_p <- create_maturities_matrix(subsubset.dataset[[group]],include_price=TRUE, version=version, SettlementCheck=SettlementCheck)
          m_p_2 <- maturities_matrix_1825(subsubset.dataset[[group]],include_price=TRUE,Actual = Actual,SettlementCheck=SettlementCheck)
          yall <- bond_yields(cf_p,m_p,m_p_2, freque=freque, mythology=mythology, Annulized=Annulized)
          
          # WEIGHTS
          dura[i,j-1]  <- duration(cf_p,m_p,yall[,2],m_p_2,freque,mythology,wmethod)[j,3]
          
          # calculate yield
          y[i,j-1] <-  yall[j,2]
          
          # calculate estimated price
          phat[i,j-1] <-  sum(cf*dt)
          
          # calculate estimated yield
          yhat[i,j-1] <-   bond_yields(c(-phat[i,j-1],cf),m_p[,j],m_p_2[,j], freque=freque, mythology=mythology, Annulized=Annulized)[1,2]
          
        }
        
      }
      
    } else if (Model == "RSplines"){
      
      
      for (i in 1:length(subset.dataset)){
        subsubset.dataset <- subset.dataset[[i]]
        
        for (j in 2:(length(subsubset.dataset[[group]]$ISIN)-1)){
          loocset <- rm_bond.couponbonds(subsubset.dataset, group, subsubset.dataset[[group]]$ISIN[[j]])
          mc_result <- estim_cs(loocset, group, matrange=matrange,rse=rse, version = version, SettlementCheck=SettlementCheck, mythology = mythology, Annulized=Annulized, Actual = Actual,wmethod=wmethod,setdf=setdf)
          
          cf <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=FALSE)[,j]
          m <- create_maturities_matrix(subsubset.dataset[[group]],include_price=FALSE, version=version, SettlementCheck=SettlementCheck)[,j]
          p[i,j-1] <- subsubset.dataset[[group]][["PRICE"]][j] + subsubset.dataset[["Great Britian"]][["ACCRUED"]][j]
          
          cf_p <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=TRUE)
          m_p <- create_maturities_matrix(subsubset.dataset[[group]],include_price=TRUE, version=version, SettlementCheck=SettlementCheck)
          m_p_2 <- maturities_matrix_1825(subsubset.dataset[[group]],include_price=TRUE,Actual = Actual,SettlementCheck=SettlementCheck)
          yall <- bond_yields(cf_p,m_p,m_p_2, freque=freque, mythology=mythology, Annulized=Annulized)
          
          # WEIGHTS
          dura[i,j-1]  <- duration(cf_p,m_p,yall[,2],m_p_2,freque,mythology,wmethod)[j,3]
          
          # calculate yield
          y[i,j-1] <-  yall[j,2]
          
          # calculate estimated price
          phat[i,j-1] <- cf%*%(predict(bs(mc_result[["x"]][[group]],df =setdf, degree = 3,intercept = TRUE),m)%*%(unname(mc_result[["alpha"]][[group]])))
          
          # calculate estimated yield
          yhat[i,j-1] <-   bond_yields(c(-phat[i,j-1],cf),m_p[,j],m_p_2[,j], freque=freque, mythology=mythology, Annulized=Annulized)[1,2]
          
        }
        
      }
      
    } else if (Model == "Exponential"){
      for (i in 1:length(subset.dataset)){
        subsubset.dataset <- subset.dataset[[i]]
        
        for (j in 2:(length(subsubset.dataset[[group]]$ISIN)-1)){
          loocset <- rm_bond.couponbonds(subsubset.dataset, group, subsubset.dataset[[group]]$ISIN[[j]])
          mc_result <- estim_nss(loocset, group,method=method, matrange=matrange, version = version, SettlementCheck=SettlementCheck, mythology = mythology, Annulized=Annulized, Actual = Actual,wmethod=wmethod,startparam=startparam,lambda=lambda,tauconstr = tauconstr,constrOptimOptions = constrOptimOptions,optimtype = "allglobal")
          
          cf <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=FALSE)[,j]
          m <- create_maturities_matrix(subsubset.dataset[[group]],include_price=FALSE, version=version, SettlementCheck=SettlementCheck)[,j]
          p[i,j-1] <- subsubset.dataset[[group]][["PRICE"]][j] + subsubset.dataset[[group]][["ACCRUED"]][j]
          
          cf_p <- create_cashflows_matrix(subsubset.dataset[[group]],include_price=TRUE)
          m_p <- create_maturities_matrix(subsubset.dataset[[group]],include_price=TRUE, version=version, SettlementCheck=SettlementCheck)
          m_p_2 <- maturities_matrix_1825(subsubset.dataset[[group]],include_price=TRUE,Actual = Actual,SettlementCheck=SettlementCheck)
          yall <- bond_yields(cf_p,m_p,m_p_2, freque=freque, mythology=mythology, Annulized=Annulized)
          
          # WEIGHTS
          dura[i,j-1]  <- duration(cf_p,m_p,yall[,2],m_p_2,freque,mythology,wmethod)[j,3]
          
          # calculate yield
          y[i,j-1] <-  yall[j,2]
          
          # calculate estimated price
          spot_rates <- spotrates(method,mc_result[["opt_result"]][[group]]$par,m,lambda)/100
          spot_rates[is.nan(spot_rates)] <- 0
          
          if (mythology){
            discount_factors <- exp(-m*spot_rates)
          } else{
            discount_factors <- 1/((1+spot_rates)^m)
          }
          
          phat[i,j-1] <- sum(cf*discount_factors)  
          
          # calculate estimated yield
          yhat[i,j-1] <-   bond_yields(c(-phat[i,j-1],cf),m_p[,j],m_p_2[,j], freque=freque, mythology=mythology, Annulized=Annulized)[1,2]
          
        }
      }
    }
  }
  returnvalue <- list(p=p,phat=phat,y=y,yhat=yhat,dura=dura,samplesize=samplesize)
  class(returnvalue) <- "dynlooc"
  return(returnvalue)
}

exportqualitylooc <- function(returnvalue,bondlist,country) UseMethod("exportqualitylooc")

exportqualitylooc.dynlooc <- function(returnvalue,bondlist,country){
  
  values <- matrix(0, nrow=length(returnvalue$samplesize), ncol=8)
  
  #RMSE Price
  for (i in 1:length(returnvalue$samplesize)){
    values[i,1] <- rmse(returnvalue$p[i,],returnvalue$phat[i,])
  }
  #MAE Price
  for (i in 1:length(returnvalue$samplesize)){
    values[i,2] <- aabse(returnvalue$p[i,],returnvalue$phat[i,])
  }
  #SMAPE Price
  for (i in 1:length(returnvalue$samplesize)){
    values[i,3] <- smape(returnvalue$p[i,],returnvalue$phat[i,])
  }
  #wRMSE Price
  for (i in 1:length(returnvalue$samplesize)){
    values[i,4] <- rmsewe(returnvalue$p[i,],returnvalue$phat[i,],returnvalue$dura[i,])
  }
  
  #RMSE Yield
  for (i in 1:length(returnvalue$samplesize)){
    values[i,5] <- rmse(returnvalue$y[i,]*100,returnvalue$yhat[i,]*100)
  }
  #MAE Yield
  for (i in 1:length(returnvalue$samplesize)){
    values[i,6] <- aabse(returnvalue$y[i,]*100,returnvalue$yhat[i,]*100)
  }
  #SMAPE Yield
  for (i in 1:length(returnvalue$samplesize)){
    values[i,7] <- smape(returnvalue$y[i,]*100,returnvalue$yhat[i,]*100)
  }
  #wRMSE Yields
  for (i in 1:length(returnvalue$samplesize)){
    values[i,8] <- rmsewe(returnvalue$y[i,]*100,returnvalue$yhat[i,]*100,returnvalue$dura[i,])
  }
  
  
  werktage <- rep(as.Date("2022-08-05"), times=length(returnvalue$samplesize))
  
  for (i in 1:length(returnvalue$samplesize)){
    werktage[i] <- bondlist[[returnvalue$samplesize[i]]][[country]][["TODAY"]]
  }
  
  colnames(values) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices", "Weighted RMSE-Prices", "RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMSE-Yields (in %)")
  rownames(values) <- as.character(werktage)
  values <- as.data.frame(values)
  write.xlsx(values,file ="QualMeasuresLOOC.xlsx",overwrite = TRUE,asTable = TRUE,rowNames = TRUE)
  class(values) <- "exportqual"
  return(values)
  
}

print.dynlooc <- function(x,...){
  values <- matrix(0, nrow=length(x$samplesize), ncol=8)
  
  #RMSE Price
  for (i in 1:length(x$samplesize)){
    values[i,1] <- rmse(x$p[i,],x$phat[i,])
  }
  #MAE Price
  for (i in 1:length(x$samplesize)){
    values[i,2] <- aabse(x$p[i,],x$phat[i,])
  }
  #SMAPE Price
  for (i in 1:length(x$samplesize)){
    values[i,3] <- smape(x$p[i,],x$phat[i,])
  }
  #wRMSE Price
  for (i in 1:length(x$samplesize)){
    values[i,4] <- rmsewe(x$p[i,],x$phat[i,],x$dura[i,])
  }
  
  #RMSE Yield
  for (i in 1:length(x$samplesize)){
    values[i,5] <- rmse(x$y[i,]*100,x$yhat[i,]*100)
  }
  #MAE Yield
  for (i in 1:length(x$samplesize)){
    values[i,6] <- aabse(x$y[i,]*100,x$yhat[i,]*100)
  }
  #SMAPE Yield
  for (i in 1:length(x$samplesize)){
    values[i,7] <- smape(x$y[i,]*100,x$yhat[i,]*100)
  }
  #wRMSE Yields
  for (i in 1:length(x$samplesize)){
    values[i,8] <- rmsewe(x$y[i,]*100,x$yhat[i,]*100,x$dura[i,])
  }
  
  RMSE_p <- mean(values[,1])
  AABSE_p <- mean(values[,2])
  sMAPE_p <- mean(values[,3])
  wRMSE_p <- mean(values[,4])
  RMSE_y <- mean(values[,5])
  AABSE_y <- mean(values[,6])
  sMAPE_y <- mean(values[,7])
  wRMSE_y <- mean(values[,8])
  
  gof <- rbind(RMSE_p,AABSE_p,sMAPE_p,wRMSE_p,RMSE_y,AABSE_y,sMAPE_y,wRMSE_y)
  colnames(gof) <- names("Country")
  rownames(gof) <- c("RMSE-Prices", "MAE-Prices","SMAPE-Prices","Weighted RMS-Prices", "RMSE-Yields (in %)", "MAE-Yields (in %)","SMAPE-Yields (in %)","Weighted RMS-Yields (in %)")
  sumry <- list(gof)
  
  cat("---------------------------------------------------\n")
  cat("Goodness of fit (OOS):\n")
  cat("---------------------------------------------------\n")
  
  print.default(format(sumry[[1]],digits=6,scientific=FALSE),quote=FALSE)
  
}


#########################################################
### Fama-Bliss without filters                      #####
#########################################################

fama_bliss <- function(bonddata,group,version = 0,SettlementCheck = FALSE,searchint=c(-1,9), tol=1e-10, Memorysaver = FALSE) UseMethod("fama_bliss")
fama_bliss.couponbonds <- function(bonddata,group,version = 0,SettlementCheck = FALSE,searchint=c(-1,9), tol=1e-10, Memorysaver = FALSE){
  
  
  m_p <- create_maturities_matrix(bonddata[[group]],include_price=TRUE,version = version,SettlementCheck=SettlementCheck)
  cf_p <- create_cashflows_matrix(bonddata[[group]],include_price=TRUE)
  
  n_group <- length(bonddata[[group]][["ISIN"]])
  forwardrate <- rep(0,n_group)
  
  colmax <- function(m) apply(m,2,max)
  maturities_bonds <- colmax(m_p)
  
  ##for the first cashflow it is equal to the YTM
  pvcashflows <- function(y) {
    t(cf_p[,1])%*%exp(-m_p[,1]*y)
  }
  # calculate bond yields
  forwardrate[1] <- uniroot(pvcashflows, searchint, tol = tol,maxiter=3000)$root 
  
  for (i in 2:n_group){
    cashnew <- cf_p[,i]
    mat <- m_p[,i]
    maturities_bonds2 <- maturities_bonds[(i-1)]
    index <- indexfinder(mat,maturities_bonds2)
    if (length(index)>1){
      index <- index[2]
    }
    discont = rep(0,index)
    discont[1]  <- exp(discont[1])
    # length((which(mat!=0))) Anzahl der nicht 0 Eintraege 
    # Till Index all Discount rates known
    
    if (index >1){
      for (j in 2:index){
        discont[j] <- exp(-discontfunction(forwardrate,m_p[j,i],maturities_bonds))
      } 
    }
    cf_mat_new <- cf_p[(index:length(cf_p[,i])),i]
    cf_mat_new[1] <- t(cf_p[(1:index),i])%*%discont
    #Hier muss ich noch cf_mat_new durch exp(-summe k = 1 bis Bigk-1 F^k*(tau^k-tau^(k-1))) abziehen
    dividendvalue <- exp(-discontfunction(forwardrate,maturities_bonds2,maturities_bonds))
    cf_mat_new[1] <- cf_mat_new[1]/dividendvalue
    
    mat_new <- mat[(index:length(cf_p[,i]))]
    mat_new[1] <- 0
    
    for (k in 2:length(mat_new)){
      if (mat_new[k] != 0){
        mat_new[k] <- mat_new[k]-maturities_bonds2
      }
    }
    
    if (!is.matrix(cf_mat_new))
      cf_mat_new <- as.matrix(cf_mat_new)
    if (!is.matrix(mat_new))
      mat_new <- as.matrix(mat_new)
    
    pvcashflows2 <- function(y) {
      t(cf_mat_new[,1])%*%exp(-mat_new[,1]*y)
    }
    
    # i loop
    # Here we need to solve for the ith forward rate
    forwardrate[i] <- uniroot(pvcashflows2, searchint, tol = tol,maxiter=3000)$root 
    
  }
  
  # maturity interval
  if(min(unname(maturities_bonds))>0.24){
    minvalue <- 0.20
  } else {
    minvalue<- round(min(unname(maturities_bonds)),2)
  }
  maxvalue<- round(max(unname(maturities_bonds)),digits = 2)-0.02
  
  
  t <- seq(minvalue,maxvalue,0.01)
  
  isin_discount <- rep(0,n_group)
  isin_spotrate <- rep(0,n_group)
  for (l in 1:n_group){
    maturities_bonds2 <- maturities_bonds[l]
    isin_discount[l] <- exp(-discontfunction(forwardrate,maturities_bonds2,maturities_bonds))
    isin_spotrate[l] <- -1/maturities_bonds2*log(isin_discount[l])
  }
  
  isin_discount<- cbind(maturities_bonds,isin_discount)
  isin_spotrate<- cbind(maturities_bonds,isin_spotrate)
  
  phat <- rep(0,n_group)
  bonderror <- rep(0,n_group)
  for (l in 1:n_group){
    mat2 <- unname(m_p[,l])
    if(length(which(mat2==0))!=0){
      mat2 <- mat2[-which(mat2==0)]
    }
    for (j in 1:length(mat2)){
      phat[l] <- phat[l] + t(cf_p[j+1,l])%*%as.vector(exp(-discontfunction(forwardrate,mat2[j],maturities_bonds)))
    }
    bonderror[l] <- cf_p[1,l] + phat[l]
  }
  
  if (!Memorysaver){  
    discount <- rep(0,length(t))
    spotrate <- rep(0,length(t))
    for (l in 1:length(t)){
      discount[l] <- exp(-discontfunction(forwardrate,t[l],maturities_bonds))
      spotrate[l] <- -1/t[l]*log(discount[l])
    }
    
    discount<- cbind(t,discount)
    spotrate<- cbind(t,spotrate)
    
    final_discount <- list()
    final_discount[[group]] <- discount
    for (k in 1:length(final_discount)) class(final_discount[[k]]) <- "ir_curve"
    class(final_discount) <- "df_curves"
    
    final_spotrate <- list()
    final_spotrate[[group]] <- spotrate
    for (k in 1:length(final_spotrate)) class(final_spotrate[[k]]) <- "ir_curve"
    class(final_spotrate) <- "spot_curves"
    
    result <- list(forwardrate = forwardrate,
                   isin_spotrate = isin_spotrate,
                   isin_discount = isin_discount,
                   discount = final_discount,
                   spot = final_spotrate,
                   phat =phat,
                   bonderror = bonderror)
    
    return(result)} else {
      result <- list(forwardrate = forwardrate,
                     isin_spotrate = isin_spotrate,
                     isin_discount = isin_discount)
      
      return(result)
    }
  
}

indexfinder <- function(mat,maturities_bonds){
  mat <- unname(mat)
  maturities_bonds <- unname(maturities_bonds)
  if(length(which(mat==0))!=0){
    mat <- mat[-which(mat==0)]
  }
  combined <- sort(c(mat,maturities_bonds))
  indexx <- c()
  for (i in 1:length(maturities_bonds)){
    index1 <- which(combined==maturities_bonds[i])
    indexx <- c(indexx,index1)
  }
  return(indexx)
}

indexfinder2 <- function(mat,maturities_bonds,val1,val2){
  mat <- unname(mat)
  maturities_bonds <- unname(maturities_bonds)
  if(length(which(mat==0))!=0){
    mat <- mat[-which(mat==0)]
  }
  combined <- sort(c(mat,maturities_bonds))
  indexx <- c()
  combined <- combined[(val1:val2)]
  for (i in 1:length(combined)){
    index1 <- which(mat==combined[i])
    indexx <- c(indexx,index1)
  }
  return(indexx)
}

discontfunction <- function(forwardr,mat,bondmat){
  mat <- unname(mat) #Mat will be a single maturity of a cashflow
  bondmat <- unname(bondmat) #Maturities of bonds
  bondmat <- c(0,bondmat)
  combined <- sort(c(mat,bondmat))
  bigK <- which(mat==combined) - 1 
  if (length(bigK)>1){
    bigK <- bigK[1]
  }
  if (bigK > 1){
    sumvalues <- forwardr[(bigK)]*(mat-bondmat[(bigK)])
    for (i in 1:(bigK-1)){
      sumvalues <- sumvalues + forwardr[i]*(bondmat[(i+1)]-bondmat[i])
    }
  } else {
    sumvalues <- forwardr[1]*mat
  }
  
  return(sumvalues)
}

fama_bliss.dyncouponbonds <- function(bonddata,group,version = 0,SettlementCheck = FALSE,searchint=c(-1,9), tol=1e-10, Memorysaver = FALSE){
  solution_list <- vector(mode = "list", length = length(bonddata))
  for (i in 1:length(bonddata)){
    # print(i)
    fama_res <- fama_bliss.couponbonds(bonddata[[i]], group=group, version = version,SettlementCheck = SettlementCheck,searchint=searchint, tol=tol,Memorysaver=Memorysaver)
    solution_list[[i]] <- fama_res
  }
  class(solution_list) <- "fama_res_dyn"
  return(solution_list)
}

fama_res_dyn_fun <- function(x,inputnodes=30,Yields=FALSE,zeit=FALSE,...){
  if (Yields & !zeit){
    set.seed(10)
    matr <- matrix(0,nrow = inputnodes,ncol = length(x))
    for (j in 1:length(x)){
      if (length(x[[j]][["forwardrate"]])<inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes-length(x[[j]][["forwardrate"]]), replace=FALSE))
        res1 <- c(unname(x[[j]][["isin_spotrate"]][,2]), unname(x[[j]][["isin_spotrate"]][samplesize,2]))
        matr[1:inputnodes,j] <- res1
        
      } else if (length(x[[j]][["forwardrate"]])>inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes, replace=FALSE))
        res1 <- unname(x[[j]][["isin_spotrate"]][samplesize,2])
        matr[1:inputnodes,j] <- res1
        
      } else if (length(x[[j]][["forwardrate"]])==inputnodes){
        samplesize <- c(0)
        res1 <- c(unname(x[[j]][["isin_spotrate"]][,2]), unname(x[[j]][["isin_spotrate"]][samplesize,2]))
        matr[1:inputnodes,j] <- res1
      }
    }
    return(matr)
  } else if (!Yields & zeit){
    set.seed(10)
    matr <- matrix(0,nrow = inputnodes,ncol = length(x))
    for (j in 1:length(x)){
      if (length(x[[j]][["forwardrate"]])<inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes-length(x[[j]][["forwardrate"]]), replace=FALSE))
        res2 <- c(unname(x[[j]][["isin_spotrate"]][,1]), unname(x[[j]][["isin_spotrate"]][samplesize,1]))
        matr[1:inputnodes,j] <- res2
        
      } else if (length(x[[j]][["forwardrate"]])>inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes, replace=FALSE))
        res2 <- unname(x[[j]][["isin_spotrate"]][samplesize,1])
        matr[1:inputnodes,j] <- res2
        
        
      } else if (length(x[[j]][["forwardrate"]])==inputnodes){
        samplesize <- c(0)
        res2 <- c(unname(x[[j]][["isin_spotrate"]][,1]), unname(x[[j]][["isin_spotrate"]][samplesize,1]))
        matr[1:inputnodes,j] <- res2
        
      }
    }
    return(matr)
  } else {
    set.seed(10)
    matr <- matrix(0,nrow = 2*inputnodes,ncol = length(x))
    for (j in 1:length(x)){
      if (length(x[[j]][["forwardrate"]])<inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes-length(x[[j]][["forwardrate"]]), replace=FALSE))
        res1 <- c(unname(x[[j]][["isin_spotrate"]][,2]), unname(x[[j]][["isin_spotrate"]][samplesize,2]))
        res2 <- c(unname(x[[j]][["isin_spotrate"]][,1]), unname(x[[j]][["isin_spotrate"]][samplesize,1]))
        matr[1:inputnodes,j] <- res1
        matr[(inputnodes+1):(2*inputnodes),j] <- res2
        
      } else if (length(x[[j]][["forwardrate"]])>inputnodes){
        samplesize <- sort(sample(1:length(x[[j]][["forwardrate"]]), inputnodes, replace=FALSE))
        res1 <- unname(x[[j]][["isin_spotrate"]][samplesize,2])
        res2 <- unname(x[[j]][["isin_spotrate"]][samplesize,1])
        matr[1:inputnodes,j] <- res1
        matr[(inputnodes+1):(2*inputnodes),j] <- res2
        
      } else if (length(x[[j]][["forwardrate"]])==inputnodes){
        samplesize <- c(0)
        res1 <- c(unname(x[[j]][["isin_spotrate"]][,2]), unname(x[[j]][["isin_spotrate"]][samplesize,2]))
        res2 <- c(unname(x[[j]][["isin_spotrate"]][,1]), unname(x[[j]][["isin_spotrate"]][samplesize,1]))
        matr[1:inputnodes,j] <- res1
        matr[(inputnodes+1):(2*inputnodes),j] <- res2
        
      }
    }
    return(matr)
  }
  
}


### AE Python
analysis_aepython <- function(export_yields,rawmodel,group,new_matrix,num_rows,version=0,Actual=0,SettlementCheck=FALSE,Annulized=TRUE,mythology=TRUE){
  export_discont <- matrix(0,nrow=dim(new_matrix)[1],ncol = dim(new_matrix)[2])
  export_discont <- exp(-export_yields*new_matrix)
  
  result <- list(list())
  
  for (i in 1:length(DynDataGB)){
    cashflows <- create_cashflows_matrix(DynDataGB[[i]][[group]],include_price=TRUE)
    m <- create_maturities_matrix(DynDataGB[[i]][[group]],include_price=TRUE, version = version, SettlementCheck=SettlementCheck)
    m2 <- maturities_matrix_1825(rawmodel[[i]][[group]],include_price=TRUE, Actual = Actual, SettlementCheck=SettlementCheck)
    y <- bond_yields(cashflows,m,m2,mythology = mythology,Annulized = Annulized)
    
    discontmatrix <- matrix(0,nrow=(dim(m)[1]-1),ncol = dim(m)[2] )
    p_hat <- rep(0,length(DynDataGB[[i]][[group]][["ISIN"]]))
    p <- -1*unname(cashflows[1,])
    if(i==1){
      a = 1
      b = sum(num_rows[1])
    } else{
      a = sum(num_rows[1:(i-1)])+1
      b = sum(num_rows[1:i])
    }
    for (j in 1:length(DynDataGB[[i]][[group]][["ISIN"]])){
      if(length(which(m[-1,j]==0))!=0){
        m_new <- m[-1,j][-which(m[-1,j]==0)]
      } else {
        m_new <- m[-1,j]
      }
      for (k in 1:length(m_new)){
        wert<-which(m_new[k] == c(as.matrix(new_matrix)[a:b,]))[1]
        discontmatrix[k,j]<-c(as.matrix(export_discont)[a:b,])[wert]
      }
      
      p_hat[j] <- t(cashflows[-1,j])%*%(as.matrix(discontmatrix[,j]))
    }
    
    p_error <- p_hat+p
    
    cashflowsnew <- cashflows
    cashflowsnew[1,] <- -p_hat
    y_hat <- bond_yields(cashflowsnew,m,m2,mythology = mythology,Annulized = Annulized)
    y_error <- cbind(y_hat[,1], y_hat[,2] - y[,2])
    
    result[[i]] <- list(p = p,
                        phat=p_hat,
                        perror=p_error,
                        y=y,
                        yhat = y_hat,
                        yerror = y_error)
    
  }
  return(result)
}

valuesfinal_aepython <- function(result){
  values <- matrix(0, nrow=length(result), ncol=6)
  
  #RMSE Price
  for (i in 1:length(result)){
    values[i,1] <- rmse(result[[i]]$p,result[[i]]$phat)
  }
  #MAE Price
  for (i in 1:length(result)){
    values[i,2] <- aabse(result[[i]]$p,result[[i]]$phat)
  }
  #SMAPE Price
  for (i in 1:length(result)){
    values[i,3] <- smape(result[[i]]$p,result[[i]]$phat)
  }
  #RMSE Yield
  for (i in 1:length(result)){
    values[i,4] <- rmse(result[[i]]$y*100,result[[i]]$yhat*100)
  }
  #MAE Yield
  for (i in 1:length(result)){
    values[i,5] <- aabse(result[[i]]$y*100,result[[i]]$yhat*100)
  }
  #SMAPE Yield
  for (i in 1:length(result)){
    values[i,6] <- smape(result[[i]]$y*100,result[[i]]$yhat*100)
  }
  
}
