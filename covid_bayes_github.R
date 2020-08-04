#################
# COVID-19 spatial Bayesian analysis
# Citation: Goldstein ND, Wheeler DC, Gustafson P, Burstyn I. A Bayesian Approach to Improving Spatial Estimates After Accounting for Misclassification Bias in Surveillance Data for COVID-19 in Philadelphia, PA. Manuscript in preparation.
# 4/30/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library("R2jags") #JAGS interface
library("MCMCvis") #visualizations
library("spdep") #spatial modeling
library("RColorBrewer") #color palette
library("cartogram") #mapping cartograms 
library("Rfast") #colMedians
library("MASS") #mvrnodm


### READ DATA ###

load("covid_spatial.RData")


### DEFINE BAYESIAN MODELS ###

jags_model1 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)                       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)              #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    logit(r[i]) <- eta0 + eta1*adi[i]                 #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model2 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)                       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)              #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i]          #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model3 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)              #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    logit(r[i]) <- eta0 + eta1*adi[i] + v[i]          #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model4 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)              #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]   #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

#sensitivity analyses on SN
jags_model1_sn = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    #sn[i] ~ dbeta(14.022, 9.681)             #sensitivity of surveillance process to identify case
    sn[i] ~ dunif(0.25, 0.75)                 #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    logit(r[i]) <- eta0 + eta1*adi[i]                 #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model2_sn = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    #sn[i] ~ dbeta(14.022, 9.681)             #sensitivity of surveillance process to identify case
    sn[i] ~ dunif(0.25, 0.75)                 #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i]          #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model3_sn = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    #sn[i] ~ dbeta(14.022, 9.681)             #sensitivity of surveillance process to identify case
    sn[i] ~ dunif(0.25, 0.75)                 #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    logit(r[i]) <- eta0 + eta1*adi[i] + v[i]          #model 3: spatial structured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]  #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}

jags_model4_sn = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #adi = area deprivation index of the ZIP code
  #zip = number of ZIP codes
  #d = binary spatial weights matrix
  #w = diagonal matrix of number of neighbors per ZIP Code
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  #sn is defined by ZIP code
  sp ~ dbeta(100, 3.02)       #specificity of surveillance process to identify case
  
  #spatial effects
  eta0 ~ dnorm(-3.5, 1/(0.5^2))               #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))              #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                      #precision of random component u
  tau_v ~ dgamma(1, 0.1)                      #precision of random component v
  phi ~ dunif(0, 1)                           #spatial correlation
  v ~ dmnorm(rep(0,zip), tau_v * (d - phi*w)) #spatial structuted random effect (note: random spatial prior specified outside for loop)
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    #sn[i] ~ dbeta(14.022, 9.681)             #sensitivity of surveillance process to identify case
    sn[i] ~ dunif(0.25, 0.75)                 #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                    #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence as function of ADI and random effect (candidate models)
    #logit(r[i]) <- eta0 + eta1*adi[i]                #model 1: base model
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    #logit(r[i]) <- eta0 + eta1*adi[i] + v[i]         #model 3: spatial structured model
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i] + v[i]   #model 4: convolution model
    
    #true infected in population
    y.t[i] ~ dbinom(r[i], t[i])
    
    #observed infected if infected
    y.o1[i] ~ dbinom(sn[i], y.t[i])     
    
    #observed infected if not infected
    y.o2[i] ~ dbinom(1-sp, t[i]-y.t[i])   
    
    #observed counts
    y.o[i] ~ sum(y.o1[i],y.o2[i])          
    
    #error in count between observed and expected
    err[i] <- y.t[i] - y.o[i]
  }
}


### RUN JAGS ###

#define neighbors based on shared boundary
sa.nb = poly2nb(philly_sf, queen=T)

#create a binary spatial weights matrix from neighbors
sa.wt = nb2mat(neighbours=sa.nb, style="B")

line_data = list("y.o"=philly_data$Positive, "t"=philly_data$Population, "adi"=as.numeric(philly_data$census_NDI), "zip"=nrow(philly_data), "d"=diag(rowSums(sa.wt)), "w"=sa.wt)

#run JAGS with benchmarking
time1 = Sys.time()
model1_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model1, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model2_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model2, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model3_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model3, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model4_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model4, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model1_sn_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model1_sn, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model2_sn_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model2_sn, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model3_sn_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model3_sn, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model4_sn_posterior = jags.parallel(data=line_data, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","phi","u","v"), model.file=jags_model4_sn, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
time2 = Sys.time()
time2-time1
rm(time1, time2)

#save
save.image("bayes_posterior.RData")


### POSTERIOR VISUALIZATIONS and DIAGNOSTICS ###

#read saved data
load("bayes_posterior.RData")

#posterior for inference
covid_posterior = model4_posterior

#obtain Gelman-Rubin convergence diagnostic (Rhat: potential scale reduction factor)
options(max.print=9999)
print(covid_posterior)

#diagnostics
MCMCtrace(model1_posterior, params=c("sn","sp","eta0","eta1","u","v","r","deviance"), wd="/", filename="Appendix 4.pdf")
MCMCtrace(model2_posterior, params=c("sn","sp","eta0","eta1","u","v","r","deviance"), wd="/", filename="Appendix 5.pdf")
MCMCtrace(model3_posterior, params=c("sn","sp","eta0","eta1","u","v","r","deviance"), wd="/", filename="Appendix 6.pdf")
MCMCtrace(model4_posterior, params=c("sn","sp","eta0","eta1","u","v","r","deviance"), wd="/", filename="Appendix 7.pdf")

#save model syntax file
writeLines(capture.output(getAnywhere(jags_model)), file("bugs_model.txt"))

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="posteriors.pdf", type="pdf", width=11, height=8.5)

#posterior prevalence plot (r) comparing 4 models
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(2,0.4))
par(new=T)
MCMCplot(model2_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(3,0.4))
par(new=T)
MCMCplot(model3_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(4,0.4))
par(new=T)
MCMCplot(model4_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(6,0.4))
points(x=1:nrow(philly_data),y=(philly_data$Positive/philly_data$Population), col="black",pch=16)
legend("topright",legend=c("Model 1", "Model 2", "Model 3", "Model 4", "Observed"), pch=16, col=c(2:4,6,1), y.intersp=0.3)

#posterior count plot (y.t)
#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","err","eta0","eta1"), labels=philly_data$ZCTA, sz_labels=1, ylab="Expected cases (95% CrI)", main="Posterior ZCTA Expected Cases (y.t)", ylim=c(0,5000))
#points(x=1:nrow(philly_data),y=philly_data$Positive, col="red",pch=16)
#legend("topright",legend=c("Expected","Observed"), pch=16, col=c("black","red"))

#posterior error plot (err)
MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Expected cases (95% CrI)", main="Posterior ZCTA Missed Cases (err)", ylim=c(0,600))

#posterior SN, SP plot
MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","phi","tau"), ylab="Accuracy (95% CrI)", main="Posterior Sensitivity (sn) and Specificity (sp)", ylim=c(0.4,1))

#posterior logit prevalence (eta0, eta1)
MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","phi","tau"), ylab="Logit true prevalance", main="Average True Prevalence (eta0) and ADI Effect (eta)", ylim=c(-4.1,1))

#posterior spatial (phi, tau)
MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","eta0","eta1"), main="Spatial Correlation (phi) and Precision (tau)", ylim=c(0,20))

#posterior non-spatial (tau)
#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","eta0","eta1","phi"), main="Precision (tau)", ylim=c(-10,10))

dev.off()


### EXTRACT POSTERIOR TRUE PREVALENCE (r) ESTIMATES ##

r_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
  chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
  
  #add zip code specific posteriors to dataframe
  r_posterior = data.frame(matrix(NA,ncol=zip,nrow=nrow(chain1)*5))
  col_mapping = colnames(r_posterior)[order(names(r_posterior))]
  for (i in 1:ncol(r_posterior)) {
    
    #determine correct column (offset is based on first r column in mcmc output)
    col_index = which(col_mapping==colnames(r_posterior)[i]) + 51
    
    #add data
    r_posterior[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
  
  return(r_posterior)
}

r_model1 = r_posterior(model1_posterior)
r_model2 = r_posterior(model2_posterior)
r_model3 = r_posterior(model3_posterior)
r_model4 = r_posterior(model4_posterior)

r_model1_sn = r_posterior(model1_sn_posterior)
r_model2_sn = r_posterior(model2_sn_posterior)
r_model3_sn = r_posterior(model3_sn_posterior)
r_model4_sn = r_posterior(model4_sn_posterior)


### POSTERIOR ESTIMATES for MANUSCRIPT ###

#SN
sn_post = data.frame(covid_posterior$BUGSoutput$summary[grep("sn",row.names(covid_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
mean(sn_post$X50)
mean(sn_post$X2.5.)
mean(sn_post$X97.5.)

sn_post[sn_post$X50==min(sn_post$X50), ]
philly_data$ZCTA[which(sn_post$X50==min(sn_post$X50))]
sn_post[sn_post$X50==max(sn_post$X50), ]
philly_data$ZCTA[which(sn_post$X50==max(sn_post$X50))]

MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","phi","sp","u","v"), labels=philly_data$ZCTA, sz_labels=1, ylab="Sensitivity (95% CrI)", ylim=c(0.2,0.9))

#SP
#sp_post = data.frame(covid_posterior$BUGSoutput$summary[grep("sp",row.names(covid_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
#mean(sp_post$X50)
#mean(sp_post$X2.5.)
#mean(sp_post$X97.5.)

data.frame(covid_posterior$BUGSoutput$summary[grep("sp",row.names(covid_posterior$BUGSoutput$summary)),], stringsAsFactors=F)

#sp_post[sp_post$X50==min(sp_post$X50), ]
#sp_post[sp_post$X50==max(sp_post$X50), ]
#philly_data$ZCTA[which(sn_post$X50==min(sn_post$X50))]
#philly_data$ZCTA[which(sn_post$X50==max(sn_post$X50))]

#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","phi","tau","sn"), labels=philly_data$ZCTA, sz_labels=1, ylab="Accuracy (95% CrI)", ylim=c(0.9,1))

#prevalence
prev_post = data.frame(covid_posterior$BUGSoutput$summary[grep("r",substr(row.names(covid_posterior$BUGSoutput$summary),1,1)),], stringsAsFactors=F)
#mean(prev_post$X50)
#mean(prev_post$X2.5.)
#mean(prev_post$X97.5.)

#citywide
r_model4$prev = (r_model1$X1*philly_data$Population[1] +
                   r_model1$X2*philly_data$Population[2] +
                   r_model1$X3*philly_data$Population[3] +
                   r_model1$X4*philly_data$Population[4] +
                   r_model1$X5*philly_data$Population[5] +
                   r_model1$X6*philly_data$Population[6] +
                   r_model1$X7*philly_data$Population[7] +
                   r_model1$X8*philly_data$Population[8] +
                   r_model1$X9*philly_data$Population[9] +
                   r_model1$X10*philly_data$Population[10] +
                   r_model1$X11*philly_data$Population[11] +
                   r_model1$X12*philly_data$Population[12] +
                   r_model1$X13*philly_data$Population[13] +
                   r_model1$X14*philly_data$Population[14] +
                   r_model1$X15*philly_data$Population[15] +
                   r_model1$X16*philly_data$Population[16] +
                   r_model1$X17*philly_data$Population[17] +
                   r_model1$X18*philly_data$Population[18] +
                   r_model1$X19*philly_data$Population[19] +
                   r_model1$X20*philly_data$Population[20] +
                   r_model1$X21*philly_data$Population[21] +
                   r_model1$X22*philly_data$Population[22] +
                   r_model1$X23*philly_data$Population[23] +
                   r_model1$X24*philly_data$Population[24] +
                   r_model1$X25*philly_data$Population[25] +
                   r_model1$X26*philly_data$Population[26] +
                   r_model1$X27*philly_data$Population[27] +
                   r_model1$X28*philly_data$Population[28] +
                   r_model1$X29*philly_data$Population[29] +
                   r_model1$X30*philly_data$Population[30] +
                   r_model1$X31*philly_data$Population[31] +
                   r_model1$X32*philly_data$Population[32] +
                   r_model1$X33*philly_data$Population[33] +
                   r_model1$X34*philly_data$Population[34] +
                   r_model1$X35*philly_data$Population[35] +
                   r_model1$X36*philly_data$Population[36] +
                   r_model1$X37*philly_data$Population[37] +
                   r_model1$X38*philly_data$Population[38] +
                   r_model1$X39*philly_data$Population[39] +
                   r_model1$X40*philly_data$Population[40] +
                   r_model1$X41*philly_data$Population[41] +
                   r_model1$X42*philly_data$Population[42] +
                   r_model1$X43*philly_data$Population[43] +
                   r_model1$X44*philly_data$Population[44] +
                   r_model1$X45*philly_data$Population[45] +
                   r_model1$X46*philly_data$Population[46] +
                   r_model1$X47*philly_data$Population[47]) / sum(philly_data$Population)
quantile(r_model4$prev, probs=c(0.025,0.5,0.975))

prev_post[prev_post$X50==min(prev_post$X50), ]
philly_data$ZCTA[which(prev_post$X50==min(prev_post$X50))]
prev_post[prev_post$X50==max(prev_post$X50), ]
philly_data$ZCTA[which(prev_post$X50==max(prev_post$X50))]

#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05))
#points(x=1:nrow(philly_data),y=(philly_data$Positive/philly_data$Population), pch=1)
#legend("topright",legend=c("Expected","Observed"), pch=c(16,1))

#measurement error
err_post = data.frame(covid_posterior$BUGSoutput$summary[grep("err",row.names(covid_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","eta0","eta1","phi","u","v"), labels=philly_data$ZCTA, sz_labels=1, ylab="Missed cases (95% CrI)", ylim=c(-100,1500))
err_post[err_post$X50==max(err_post$X50), ]
philly_data$ZCTA[which(err_post$X50==max(err_post$X50))]
prev_post[philly_data$ZCTA==19120, ]


### POSTERIOR PREVALENCE PLOT ###

#create the boxplot dataframe
boxplot_df = data.frame("r"=c(c(r_model1$X1, r_model2$X1, r_model3$X1, r_model4$X1),
                              c(r_model1$X2, r_model2$X2, r_model3$X2, r_model4$X2),
                              c(r_model1$X3, r_model2$X3, r_model3$X3, r_model4$X3),
                              c(r_model1$X4, r_model2$X4, r_model3$X4, r_model4$X4),
                              c(r_model1$X5, r_model2$X5, r_model3$X5, r_model4$X5),
                              c(r_model1$X6, r_model2$X6, r_model3$X6, r_model4$X6),
                              c(r_model1$X7, r_model2$X7, r_model3$X7, r_model4$X7),
                              c(r_model1$X8, r_model2$X8, r_model3$X8, r_model4$X8),
                              c(r_model1$X9, r_model2$X9, r_model3$X9, r_model4$X9),
                              c(r_model1$X10, r_model2$X10, r_model3$X10, r_model4$X10),
                              c(r_model1$X11, r_model2$X11, r_model3$X11, r_model4$X11),
                              c(r_model1$X12, r_model2$X12, r_model3$X12, r_model4$X12),
                              c(r_model1$X13, r_model2$X13, r_model3$X13, r_model4$X13),
                              c(r_model1$X14, r_model2$X14, r_model3$X14, r_model4$X14),
                              c(r_model1$X15, r_model2$X15, r_model3$X15, r_model4$X15),
                              c(r_model1$X16, r_model2$X16, r_model3$X16, r_model4$X16),
                              c(r_model1$X17, r_model2$X17, r_model3$X17, r_model4$X17),
                              c(r_model1$X18, r_model2$X18, r_model3$X18, r_model4$X18),
                              c(r_model1$X19, r_model2$X19, r_model3$X19, r_model4$X19),
                              c(r_model1$X20, r_model2$X20, r_model3$X20, r_model4$X20),
                              c(r_model1$X21, r_model2$X21, r_model3$X21, r_model4$X21),
                              c(r_model1$X22, r_model2$X22, r_model3$X22, r_model4$X22),
                              c(r_model1$X23, r_model2$X23, r_model3$X23, r_model4$X23),
                              c(r_model1$X24, r_model2$X24, r_model3$X24, r_model4$X24),
                              c(r_model1$X25, r_model2$X25, r_model3$X25, r_model4$X25),
                              c(r_model1$X26, r_model2$X26, r_model3$X26, r_model4$X26),
                              c(r_model1$X27, r_model2$X27, r_model3$X27, r_model4$X27),
                              c(r_model1$X28, r_model2$X28, r_model3$X28, r_model4$X28),
                              c(r_model1$X29, r_model2$X29, r_model3$X29, r_model4$X29),
                              c(r_model1$X30, r_model2$X30, r_model3$X30, r_model4$X30),
                              c(r_model1$X31, r_model2$X31, r_model3$X31, r_model4$X31),
                              c(r_model1$X32, r_model2$X32, r_model3$X32, r_model4$X32),
                              c(r_model1$X33, r_model2$X33, r_model3$X33, r_model4$X33),
                              c(r_model1$X34, r_model2$X34, r_model3$X34, r_model4$X34),
                              c(r_model1$X35, r_model2$X35, r_model3$X35, r_model4$X35),
                              c(r_model1$X36, r_model2$X36, r_model3$X36, r_model4$X36),
                              c(r_model1$X37, r_model2$X37, r_model3$X37, r_model4$X37),
                              c(r_model1$X38, r_model2$X38, r_model3$X38, r_model4$X38),
                              c(r_model1$X39, r_model2$X39, r_model3$X39, r_model4$X39),
                              c(r_model1$X40, r_model2$X40, r_model3$X40, r_model4$X40),
                              c(r_model1$X41, r_model2$X41, r_model3$X41, r_model4$X41),
                              c(r_model1$X42, r_model2$X42, r_model3$X42, r_model4$X42),
                              c(r_model1$X43, r_model2$X43, r_model3$X43, r_model4$X43),
                              c(r_model1$X44, r_model2$X44, r_model3$X44, r_model4$X44),
                              c(r_model1$X45, r_model2$X45, r_model3$X45, r_model4$X45),
                              c(r_model1$X46, r_model2$X46, r_model3$X46, r_model4$X46),
                              c(r_model1$X47, r_model2$X47, r_model3$X47, r_model4$X47)),
                        "zcta"=rep(philly_data$ZCTA[1:zip], each=nrow(r_model1)*4),
                        "model"=rep(rep(c("Model 1","Model 2","Model 3", "Model 4"), each=nrow(r_model1)), zip),
                        stringsAsFactors=F)

#x-axis zip codes
plot_labels = rep(philly_data$ZCTA, each=4)
plot_labels[seq(2,188, by=4)] = ""
plot_labels[seq(3,188, by=4)] = ""
plot_labels[seq(4,188, by=4)] = ""

#export @ 1600 x 600
boxplot(r ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", xlab="", ylab="Prevalence", col=c("#F0E442", "#56B4E9", "#009E73","#D55E00"))
axis(side=1, labels=philly_data$ZCTA, at=seq(2.5,188, by=4), las=2, tick=F)
rect(seq(0.5,188,by=8),rep(-1,24),seq(4.5,188.5,by=8),1,col = rgb(0.5,0.5,0.5,1/5), border=NA)
points(x=seq(2.5,188,by=4),y=(philly_data$Positive/philly_data$Population), col="black", pch=8, cex=1.5)
legend("topleft", legend=c("Model 1","Model 2", "Model 3", "Model 4", "Observed"), fill=c(c("#F0E442", "#56B4E9", "#009E73","#D55E00"), NA), border=c("black","black","black","black",NA), y.intersp=0.6, cex=0.8, pch=c(NA,NA,NA,NA,8))


### PRIOR and POSTERIOR PLOTS ##

par(mfrow=c(3,2))

#phi
plot(NULL, main="A) phi", xlab="Value", ylab="Density", yaxt='n', xlim=c(-0.1,1.1), ylim=c(0,1.5))

#obtain posterior
phi = c(as.data.frame(as.mcmc(covid_posterior)[[1]])$phi,
        as.data.frame(as.mcmc(covid_posterior)[[2]])$phi,
        as.data.frame(as.mcmc(covid_posterior)[[3]])$phi,
        as.data.frame(as.mcmc(covid_posterior)[[1]])$phi,
        as.data.frame(as.mcmc(covid_posterior)[[1]])$phi)

#add lines
lines(density(phi), lwd=2, col="gray")
lines(density(runif(covid_posterior$BUGSoutput$n.sims, 0, 1)), lwd=2, col="black")
#legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
rm(phi)

#sp
plot(NULL, main="B) specificity", xlab="Value", ylab="Density", yaxt='n', xlim=c(0.95,1), ylim=c(0,500))

#obtain posterior
sp = c(as.data.frame(as.mcmc(covid_posterior)[[1]])$sp,
       as.data.frame(as.mcmc(covid_posterior)[[2]])$sp,
       as.data.frame(as.mcmc(covid_posterior)[[3]])$sp,
       as.data.frame(as.mcmc(covid_posterior)[[1]])$sp,
       as.data.frame(as.mcmc(covid_posterior)[[1]])$sp)

#add lines
lines(density(sp), lwd=2, col="gray")
lines(density(rbeta(covid_posterior$BUGSoutput$n.sims, 100, 3.02)), lwd=2, col="black")
#legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
rm(sp)

#sn
plot(NULL, main="C) sensitivity", xlab="Value", ylab="Density", yaxt='n', xlim=c(0.2,0.9), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first r column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + 98
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, xlim=c(0.95, 1), yaxt='n', col="gray")
}
rm(i, sn)

lines(density(rbeta(covid_posterior$BUGSoutput$n.sims, 14.022, 9.681)), lwd=2, col="black")
#legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))

#u
plot(NULL, main="D) u", xlab="Value", ylab="Density", yaxt='n', xlim=c(-20,20), ylim=c(0,2))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])

#add zip code specific posteriors to dataframe
u = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
col_mapping = colnames(u)[order(names(u))]
for (i in 1:ncol(u)) {
  
  #determine correct column (offset is based on first r column in mcmc output)
  col_index = which(col_mapping==colnames(u)[i]) + 146
  
  #add data
  u[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(u[,i]), lwd=1, xlim=c(0.95, 1), yaxt='n', col="gray")
}
rm(i, u)

lines(density(rnorm(covid_posterior$BUGSoutput$n.sims, 0, median(rgamma(covid_posterior$BUGSoutput$n.sims, 1, 0.1)))), lwd=2, col="black")
#legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))

#v
plot(NULL, main="E) v", xlab="Value", ylab="Density", yaxt='n', xlim=c(-20,20), ylim=c(0,2))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])

#add zip code specific posteriors to dataframe
v = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
col_mapping = colnames(v)[order(names(v))]
for (i in 1:ncol(v)) {
  
  #determine correct column (offset is based on first r column in mcmc output)
  col_index = which(col_mapping==colnames(v)[i]) + 193
  
  #add data
  v[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(v[,i]), lwd=1, xlim=c(0.95, 1), yaxt='n', col="gray")
}
rm(i, v)

lines(density(mvrnorm(covid_posterior$BUGSoutput$n.sims, mu=rep(0,zip), Sigma=median(rgamma(covid_posterior$BUGSoutput$n.sims, 1, 0.1)) * (d - median(runif(covid_posterior$BUGSoutput$n.sims, 0, 1))*w))), lwd=2, col="black")
#legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))


### SENSITIVITY ANALYSIS CORRELATION MATRIX ###

par(mar=rep(2,4), mfrow=c(2,2))

plot(colMedians(as.matrix(r_model1_sn)), colMedians(as.matrix(r_model1)), pch=16, xlim=c(0.01,0.06), ylim=c(0,0.06), xlab="", ylab="", main="A) Model 1")
text(0.0125,0.05, expression(paste(rho)))
text(0.016, 0.05, paste("=","",round(cor(colMedians(as.matrix(r_model1)), colMedians(as.matrix(r_model1_sn))), 2)))

plot(colMedians(as.matrix(r_model2_sn)), colMedians(as.matrix(r_model2)), pch=16, xlim=c(0.01,0.06), ylim=c(0,0.06), xlab="", ylab="", main="B) Model 2")
text(0.0125,0.05, expression(paste(rho)))
text(0.018, 0.05, paste("=","",round(cor(colMedians(as.matrix(r_model2)), colMedians(as.matrix(r_model2_sn))), 2)))

plot(colMedians(as.matrix(r_model3_sn)), colMedians(as.matrix(r_model3)), pch=16, xlim=c(0.01,0.06), ylim=c(0,0.06), xlab="", ylab="", main="C) Model 3")
text(0.0125,0.05, expression(paste(rho)))
text(0.018, 0.05, paste("=","",round(cor(colMedians(as.matrix(r_model3)), colMedians(as.matrix(r_model3_sn))), 2)))

plot(colMedians(as.matrix(r_model4_sn)), colMedians(as.matrix(r_model4)), pch=16, xlim=c(0.01,0.06), ylim=c(0,0.06), xlab="", ylab="", main="D) Model 4")
text(0.0125,0.05, expression(paste(rho)))
text(0.018, 0.05, paste("=","",round(cor(colMedians(as.matrix(r_model4)), colMedians(as.matrix(r_model4_sn))), 2)))


### SPATIAL ANALYSIS of POSTERIOR ###

#set CRS; https://epsg.io/3652
philly_sf_proj = st_transform(st_as_sf(philly_sf), crs=3652)

#merge case data with shapefile
philly_sf_joined = merge(x=philly_sf_proj, y=philly_data, by.x="ZCTA5CE10", by.y="ZCTA", all.x=T, duplicateGeoms=T)

#define neighbors based on shared boundary
sa.nb = poly2nb(philly_sf_joined, queen=T)

#create a representation of the binary weights
sa.wt = nb2listw(neighbours=sa.nb, style="B")

#moran's spatial autocorrelation on missed infections per capita
moran.mc(prev_post$X50., listw=sa.wt, nsim=10000)

#spatial correlogram, shows autocorrelation by order of neighbor (1=neighbor, 2=neighbor's neighbor, etc)
plot(sp.correlogram(neighbours=sa.nb,var=prev_post$X50.,order=4,method="I",style="B",zero.policy=T), main="")

#choropleth map: see https://edzer.github.io/sp/ for helpful tips
philly_sf_joined$prev_post = prev_post$X50.
sp_arrow = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(2730000,220000), scale = 10000)
sp_scale = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(2725000,215000), scale = 21120, fill=c("transparent","black"))
sp_scale_text1 = list("sp.text", c(2726000,211000), "0")
sp_scale_text2 = list("sp.text", c(2742000,211000), "4 mi")
spplot(as_Spatial(philly_sf_joined), "prev_post", cuts=7, col.regions=brewer.pal(8, "Reds"), sp.layout=list(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2))
rm(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2)

#cartogram map
philly_sf_joined$err = err_post$X50. / philly_sf_joined$Population
carto = cartogram_cont(philly_sf_joined, "err", itermax=10)

par(mar=rep(0.1,4), mfrow=c(1,2))
plot(philly_sf_joined$geometry)
text(t(sapply(slot(as_Spatial(philly_sf_joined), "polygons"), function(i) slot(i, "labpt"))), cex=0.6, labels=philly_sf_joined$ZCTA5CE10)
plot(carto$geometry)

