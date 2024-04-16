# # DEFINING THE DATA:
# myData <-  matrix (
#   c(64.0, 62.3,   NA,   NA, 64.8, 57.5,   NA, 70.2, 63.9, 71.1, 
#     66.5, 68.1,   NA, 75.1, 64.6, 69.2, 68.1,   NA, 63.2,   NA, 
#     64.1, 71.5, 76.0, 69.7, 73.3, 61.7, 66.4, 65.7, 68.3, 66.9,
#     136.4,215.1,173.6,117.3,123.3, 96.5,178.3,191.1,158.0,193.9, 
#     127.1,147.9,119.0,204.4,143.4,124.4,140.9,164.7,139.8,110.2, 
#     134.1,193.6,180.0,155.0,188.2,187.4,139.2,147.9,178.6,111.1) ,
#   nrow=30  )
# colnames(myData) <- c("height","weight")
# myData <- as.data.frame(myData)

# this index will help setup priors and let us look at posterior values for missing x's
mIdx <- ifelse( is.na(data$weight) , 1 , 0)
mIdx <- sapply( 1:length(mIdx), 
                function(n) mIdx[n]*sum(mIdx[1:n]))
# result: mIdx = 
#              0, 0, 1, 2, 0, 0, 3, 0, 0, 0, 
#              0, 0, 4, 0, 0, 0, 0, 5, 0, 6, 
#              0, 0, 0, 0, 0, 0, 0, 0, 0, 0

# add missing index to myData
data$mIdx <- mIdx



# DATA PREP:
x = data[,"weight"]
# x = myData[,"height"]
# meanY = mean(y,na.rm=TRUE) #mean weight
meanX = mean(x,na.rm=TRUE) #mean weight
# sdY = sd(y,na.rm=TRUE)     
sdX = sd(x,na.rm=TRUE)
Ntotal = length(x)
zx <- NULL
for ( i in 1:Ntotal ) {
  zx[i] <- ifelse ( mIdx[i]==0, ( x[i] - meanX ) / sdX , x[i] ) # skips NA's
  # zy[i] <- ( y[i] - meanY ) / sdY
}
# Specify the data list for JAGS
dataList = list(
  zx = zx ,
  # zy = zy ,
  # meanY = meanY ,
  meanX = meanX ,
  # sdY = sdY ,
  sdX = sdX ,
  Ntotal = Ntotal
)



# THE MODEL

model_string <- "
# Specify the model for standardized data:
model {
  for ( i in 1:Ntotal ) {
    zy[i] ~ dt( zbeta0 + zbeta1 * zx[i] , 1/zsigma^2 , nu )
  }
  # prior for imputing missing zx's
  for (i in 1:Ntotal){
    zx[i] ~ dnorm( 0 , 1 )
    # Estimated covariates
    ximp[i] <- zx[i]*sdX + meanX
  }
  # Priors vague on standardized scale:
  zbeta0 ~ dnorm( 0 , 1/(10)^2 )  
  zbeta1 ~ dnorm( 0 , 1/(10)^2 )
  zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
  nu ~ dexp(1/30.0)
  # Transform back to original scale:
  beta1 <- zbeta1 * sdY / sdX  
  beta0 <- zbeta0 * sdY  + meanY - zbeta1 * meanX * sdY / sdX 
  sigma <- zsigma * sdY
}
"

# INITIALIZE VALUES
# values hardcoded for simplicity
zbeta0 = 0
zbeta1 = 0.5
zsigma = 1
nu = 30
# initial values for missing x data:
xInit = rep( NA , length(x) )
xInit[3] <- 68 ; xInit[4]<- 64 ; xInit[7] <- 68
xInit[13] <- 64 ; xInit[18] <- 68 ; xInit[20] <- 64
initsList = list( zbeta0=zbeta0 , zbeta1=zbeta1 , 
                  zsigma=zsigma , nu = nu , x=xInit )

nChains=2

library(rjags)
jagsModel = jags.model( textConnection(model_string) , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )

update(jagsModel,5000)

# Running the model
samp=coda.samples(jagsModel, variable.names=c("beta0", "beta1", "sigma", "nu", "ximp"), 
                  n.iter=2000)

plot(samp)