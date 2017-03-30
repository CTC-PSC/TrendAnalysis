rm(list=ls()) #Clear the workspace

library(MASS); # Loads miscellaneous functions (ginv, etc.)
library(gplots)
#library(foreign)

#This program is an adaptation of trend analysis code from "A better way to estimate population trends"
#by Humbert et al. (2009), Oikos 118:1940:1946; the principal adaptation is that is meant to loop over
#many stocks simultaneously, deal with gappy/short time series, etc.


#-------------------------------------------------------------------------------
#   1. Set specs, paths, etc. and read in the abundance dataset(formatted as in template)
#-------------------------------------------------------------------------------
this_is_the_place<-file.path(choose.dir()) #Get file path string for writing output file 
MasterDataFrame<-read.table(file.choose(),header=TRUE) #direct it to the text file
stocklist<-names(MasterDataFrame)[3:length(names(MasterDataFrame))] #vector of pop'n names
#how many pops are involved
numstks<-length(stocklist)
#-------------------------------------------------------------------------------




#-------------------------------------------------------------------------------
#   2. Create output data frame
#-------------------------------------------------------------------------------
nm<-c('Stock','YearEnd',
      'OE.mu','OE.ssq','OE.tsq','OE.xo',
      'PN.mu','PN.ssq','PN.tsq','PN.xo',
      'SSml.mu','SSml.ssq','SSml.tsq','SSml.xo',
      'SSreml.mu','SSreml.ssq','SSreml.tsq','SSreml.xo',
      'SSreml.mUCB','SSreml.mLCB')

# Output header name prefixes and defs
#OE = expon. growth observer error model, log(N) vs. t
#PN = expon. growth process noise model, traditional Dennis model
#SS = state space model, estimation of trend (mu) while accounting for PN & OE
#SS computed with ML and REML variants
#mu = trend parameter (log(lambda))
#ssq = process error
#tsq = observer error
#x0 = initial values for SS models (average of OE and PN)
#mUCB, mLCB = upper and lower 95% CBs on mu 

#Create an output data frame (set number of rows to number of pops)
OutputDFLT <- as.data.frame(
  matrix(nrow = numstks, ncol = 20,
         dimnames = list(NULL, nm)))
#-------------------------------------------------------------------------------




#----------------------------------------------------------------------
#   3. DEFINE FUNCTIONS FOR COMPUTING ML & REML LOG-LIKELIHOODS
#----------------------------------------------------------------------
# ML objective function «negloglike.ml» is negative of log-likelihood;
# the Nelder-Mead optimization routine in R, «optim», is a minimization
# routine. The ML objective function uses equations 24-26 from Dennis et
# al. (2006). The three function arguments are: theta, vector of
# parameters (transformed to the real line), yt, vector of time series
# observations, and tt, vector of observation times.
negloglike.ml=function(theta,yt,tt)
{
  muu=theta[1];
  sigmasq=exp(theta[2]); # Constrains ssq > 0.
  tausq=exp(theta[3]); # Constrains tsq > 0.
  xzero=theta[4];
  q=length(yt)-1;
  qp1=q+1;
  yt=matrix(yt,nrow=qp1,ncol=1);
  vx=matrix(0,qp1,qp1);
  for (ti in 1:q)
  {
    vx[(ti+1):qp1,(ti+1):qp1]=matrix(1,1,(qp1-ti))*tt[ti+1];
  }
  Sigma.mat=sigmasq*vx;
  Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
  diag(Itausq)=rep(tausq,q+1);
  V=Sigma.mat+Itausq;
  mu=matrix((xzero+muu*tt),nrow=qp1,ncol=1);
  ofn=((qp1)/2)*log(2*pi)+(0.5*log(det(V)))+
    (0.5*(t(yt-mu)%*%ginv(V)%*%(yt-mu)));
  return(ofn);
}
# REML objective function „negloglike.reml“ is negative of log-likelihood
# for second differences of the log-scale observations. The REML objective
# function uses equations A18-A25 from Humbert et al. (2009). The three
# function arguments are: theta, vector of parameters (transformed to the
# real line), yt, vector of time series observations (log scale), and
# tt, vector of observation times. Function performs the differencing.
negloglike.reml=function(theta,yt,tt)
{
  sigsq=exp(theta[1]); # Constrains ssq > 0.
  tausq=exp(theta[2]); # Constrains tsq > 0.
  q=length(yt)-1;
  qp1=q+1;
  vx=matrix(0,qp1,qp1);
  for (ti in 1:q)
  {
    vx[(ti+1):qp1,(ti+1):qp1]=matrix(1,1,(qp1-ti))*tt[ti+1];
  }
  Sigma.mat=sigsq*vx;
  Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
  diag(Itausq)=rep(tausq,q+1);
  V=Sigma.mat+Itausq;
  ss=tt[2:qp1]-tt[1:q];
  D1mat=cbind(-diag(1/ss),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/ss));
  D2mat=cbind(-diag(1,q-1),matrix(0,q-1,1))+
    cbind(matrix(0,q-1,1),diag(1,q-1));
  Phi.mat=D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat);
  wt=(yt[2:qp1]-yt[1:q])/ss;
  ut=wt[2:q]-wt[1:q-1];
  ofn=(q/2)*log(2*pi)+(0.5*log(det(Phi.mat)))+
    (0.5*(ut%*%ginv(Phi.mat)%*%ut));
  return(ofn);
}
#----------------------------------------------------------------------




#----------------------------------------------------------------------
#   4. The main data processing and estimation loop
#----------------------------------------------------------------------
### Will work the entire dataset, pop'n by pop'n

sets<-1 #counting variable for output file (moves a single row with each stock/set combo in output file)
i=1

while(i <= numstks)
{
      my_stock<-MasterDataFrame[c(1,2,i+2)]
      stockname<-names(MasterDataFrame[i+2])
      

      #*************************************************************************************
      #Uncomment and specify a particular year if subsetting for some reason, e.g., pre-post agreement periods
      #my_stock<-subset(my_stock,my_stock$Year>=1999)   
      #*************************************************************************************
      
    
      #-----------------------------------------------------------------------
      #This next chunky block (between "<<<<<" rows) is ugly but little more than some fussy subsetting
      #and series modification is needed to deal with datasets of different shapes/sizes
      #-----------------------------------------------------------------------
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      for(k in 1:length(my_stock[,1])) #Skip to first record for stock in ragged array
      {
            if(is.na(my_stock[k,3])==TRUE)
              {my_stock$true[k]=1} else
              {
                my_stock$true[k]=0 #flag for not NA
                if(my_stock[k,3]==0)
                {my_stock[k,3]=1} #***If there's a 0 count in the series, set it to 1 b/c ln(0) is undef'd 
              }
            
      }
      
      #now dump trailing NAs
      LastYear<-NA
      z<-length(my_stock[,1])
      while(is.na(LastYear)==TRUE) #Skip to first record for stock in ragged array
      {
            if(is.na(my_stock[z,3])==FALSE)
            {LastYear<-my_stock$Year[z]}
            z<-z-1
      }
      
      my_stock<-subset(my_stock,Year<=LastYear) #clip ragged end
      
      
      gaps<-FALSE
      gaptest<-rep(0,length(my_stock[,1])-1)
    
      for(k in 1:(length(my_stock[,1])-1))
      {
            gaptest[k]<-my_stock$true[k]-my_stock$true[k+1]
      }
      if(min(gaptest)<0)
      {
            gaps<-TRUE #Set boolean to true if there's a gap in the series
      } 
      #-----------------------------------------------------------------------

      #Dump NAs and rename to my_data for compatibility with original Humbert code
      my_data <-subset(my_stock,my_stock$true==0)
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
      #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
      
        
      #FINALLY, Begin Humbert Et Al. Code
        
      ##----------------------------------------------------------------------
      # a. USER INPUT SECTION
      #----------------------------------------------------------------------
      Observed.t <- my_data[,3]
      Time.t <- my_data$Year
      EndYear<-max(Time.t)
      #print.table(cbind(Observed.t,Time.t)) #uncomment if you want to see it
        
      #----------------------------------------------------------------------
      # b. PROGRAM INITIALIZATION SECTION
      #----------------------------------------------------------------------
      T.t=Time.t-Time.t[1]; # Time starts at zero.
      Y.t=log(Observed.t); # Log-transform the observations.
      q=length(Y.t)-1; # Number of time series transitions, q.
      qp1=q+1; # q+1 gets used a lot, too.
      S.t=T.t[2:qp1]-T.t[1:q]; # Time intervals.
      m=rep(1,qp1); # Will contain Kalman means for Kalman calculations.
      v=rep(1,qp1); # Will contain variances for Kalman calculations.
        
    
      #----------------------------------------------------------------------
      # c. SECTION FOR CALCULATING EGOE AND EGPN ESTIMATES
      # (FOR USE AS INITIAL VALUES)
      #----------------------------------------------------------------------
      # The EGOE (Exponential Growth Observation Error) estimates
      Ybar=mean(Y.t);
      Tbar=mean(T.t);
      mu.egoe=sum((T.t-Tbar)*(Y.t-Ybar))/sum((T.t-Tbar)*(T.t-Tbar));
      x0.egoe=Ybar-mu.egoe*Tbar;
      ssq.egoe=0;
      Yhat.egoe=x0.egoe+mu.egoe*T.t;
      tsq.egoe=sum((Y.t-Yhat.egoe)*(Y.t-Yhat.egoe))/(q-1);

      # The EGPN (Exponential Growth Process Noise) estimates
      Ttr=sqrt(S.t);
      Ytr=(Y.t[2:qp1]-Y.t[1:q])/Ttr;
      mu.egpn=sum(Ttr*Ytr)/sum(Ttr*Ttr);
      Ytrhat=mu.egpn*Ttr;
      ssq.egpn=sum((Ytr-Ytrhat)*(Ytr-Ytrhat))/(q-1);
      tsq.egpn=0;
      x0.egpn=Y.t[1];
      # Initial values for EGSS are averages of EGOE and EGPN values
      mu0=(mu.egoe+mu.egpn)/2; # For ML only
      ssq0=ssq.egpn/2; # For ML and REML
      tsq0=tsq.egoe/2; # For ML and REML
      x00=x0.egoe; # For ML only
      
      # To set different initial values for iterations, enter manually a value
      # after the equal sign of the concern parameter instead of the
      # automatically generated value. Then run again the line and the program
      # section 5 below.
      # Initial values near the EGOE and EGPN models are good for exploring
      # possible alternative local maxima. The values which produce the largest
      # log-likelihood should be used. To see the log-likelihood for the REML
      # estimates type:
      # EGSSreml$value[1];
      # See Dennis et al. 2006 for more details.

      # collect for file writing
      parms.egoe=c(mu.egoe,ssq.egoe,tsq.egoe,x0.egoe); # --
      parms.egpn=c(mu.egpn,ssq.egpn,tsq.egpn,x0.egpn); # --
        
      #----------------------------------------------------------------------
      # d. SECTION FOR CALCULATING ML & REML PARAMETER ESTIMATES (STATE SPACE MODEL)
      #----------------------------------------------------------------------
        
      #This section will be skipped if there are any gaps in your time series.
      #If your series has only one gap, you should consider interpolating in 
      #your data file and reloading it... 
         
      if(gaps==FALSE)
          {
              # The ML estimates.
              EGSSml=optim(par=c(mu0,log(ssq0),log(tsq0),x00),
                           negloglike.ml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
              params.ml=c(EGSSml$par[1],exp(EGSSml$par[2]),exp(EGSSml$par[3]),
                          EGSSml$par[4]);
              lnlike.ml=-EGSSml$value[1];
              AIC.egss=-2*lnlike.ml+2*length(params.ml);
              mu.ml=params.ml[1]; # These are the ML estimates.
              ssq.ml=params.ml[2]; # --
              tsq.ml=params.ml[3]; # --
              x0.ml=params.ml[4]; # --
  
              # The REML estimates.
              EGSSreml=optim(par=c(log(ssq0),log(tsq0)),
                             negloglike.reml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
              params.reml=c(exp(EGSSreml$par[1]),exp(EGSSreml$par[2]))
              ssq.reml=params.reml[1]; # These are the REML estimates.
              tsq.reml=params.reml[2]; # --
              vx=matrix(0,qp1,qp1);
              for (ti in 1:q)
                {
                  vx[(ti+1):qp1,(ti+1):qp1]=matrix(1,1,(qp1-ti))*T.t[ti+1];
                }
              Sigma.mat=ssq.reml*vx;
              Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
              diag(Itausq)=rep(tsq.reml,q+1);
              V=Sigma.mat+Itausq;
              D1mat=cbind(-diag(1/S.t),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/S.t));
              V1mat=D1mat%*%V%*%t(D1mat);
              W.t=(Y.t[2:qp1]-Y.t[1:q])/S.t;
              j1=matrix(1,q,1);
              V1inv=ginv(V1mat);
              mu.reml=(t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1);
              j=matrix(1,qp1,1);
              Vinv=ginv(V);
              x0.reml=(t(j)%*%Vinv%*%(Y.t-mu.reml*T.t))/(t(j)%*%Vinv%*%j);
              Var_mu.reml=1/(t(j1)%*%V1inv%*%j1); # Variance of mu
              mu_hi.reml=mu.reml+1.96*sqrt(Var_mu.reml); # 95% CI for mu
              mu_lo.reml=mu.reml-1.96*sqrt(Var_mu.reml); # --
  
              # Calculate estimated population sizes for EGSS model
              # with Kalman filter, for plotting.
              #
              ###Choose ML or REML estimates here for calculating model values ###
              ###for plotting (by commenting out the unwanted, default is REML)###
              mu=mu.ml; ssq=ssq.ml; tsq=tsq.ml; x0=x0.ml;
              # mu=mu.reml; ssq=ssq.reml; tsq=tsq.reml; x0=x0.reml;
              m[1]=x0; # Initial mean of Y(t).
              v[1]=tsq; # Initial variance of Y(t).
              for (ti in 1:q) # Loop to generate estimated population abundances
                { # using Kalman filter (see equations 6 & 7,
                  # Dennis et al. (2006)).
                  m[ti+1]=mu+(m[ti]+((v[ti]-tsq)/v[ti])*(Y.t[ti]-m[ti]));
                  v[ti+1]=tsq*((v[ti]-tsq)/v[ti])+ssq+tsq;
                }
              #collect for file writing
              parms.reml=c(mu.reml,ssq.reml,tsq.reml,x0.reml); # --
              parms.ml=c(mu.ml,ssq.ml,tsq.ml,x0.ml); # --
          } else
          {
              #NAs for gappy data
              parms.reml=c(NA,NA,NA,NA); # --
              parms.ml=c(NA,NA,NA,NA); # --
              mu_hi.reml<-NA
              mu_lo.reml<-NA        
          }
      
    
        
      #Write the results to dataframe:
      OutputDFLT$Stock[sets]=stockname
      OutputDFLT$YearEnd[sets]<-EndYear
      OutputDFLT$OE.mu[sets]<-parms.egoe[1]
      OutputDFLT$OE.ssq[sets]<-parms.egoe[2]
      OutputDFLT$OE.tsq[sets]<-parms.egoe[3]
      OutputDFLT$OE.xo[sets]<-parms.egoe[4]
      OutputDFLT$PN.mu[sets]<-parms.egpn[1]
      OutputDFLT$PN.ssq[sets]<-parms.egpn[2]
      OutputDFLT$PN.tsq[sets]<-parms.egpn[3]
      OutputDFLT$PN.xo[sets]<-parms.egpn[4]
      OutputDFLT$SSml.mu[sets]<-parms.ml[1]
      OutputDFLT$SSml.ssq[sets]<-parms.ml[2]
      OutputDFLT$SSml.tsq[sets]<-parms.ml[3]
      OutputDFLT$SSml.xo[sets]<-parms.ml[4]
      OutputDFLT$SSreml.mu[sets]<-parms.reml[1]
      OutputDFLT$SSreml.ssq[sets]<-parms.reml[2]
      OutputDFLT$SSreml.tsq[sets]<-parms.reml[3]
      OutputDFLT$SSreml.xo[sets]<-parms.reml[4]
      OutputDFLT$SSreml.mUCB[sets]<-mu_hi.reml
      OutputDFLT$SSreml.mLCB[sets]<-mu_lo.reml
        
      sets<-sets+1
      i<-i+1

}


#----------------------------------------------------------------------
#   5. DONE! Now, write the output file!
#----------------------------------------------------------------------
stamp<-format(Sys.time(), "%a %b %d") #DateStampForOutputFile
write.csv(OutputDFLT,file = paste(this_is_the_place,"\\",stamp," TrendOutput.csv",sep=""),row.names=FALSE)
#----------------------------------------------------------------------











