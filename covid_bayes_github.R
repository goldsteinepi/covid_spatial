#################
# COVID-19 spatial Bayesian analysis
# Citation: Goldstein ND, Wheeler DC, Gustafson P, Burstyn I. A Bayesian Approach to Improving Spatial Estimates of Prevalence of COVID-19 After Accounting for Misclassification Bias in Surveillance Data in Philadelphia, PA. Manuscript in preparation.
# 4/30/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library("R2jags") #JAGS interface
library("MCMCvis") #visualizations
library("spdep") #spatial modeling
library("RColorBrewer") #color palette
library("cartogram") #mapping cartograms 
library("Rfast") #colMedians
#library("MASS") #mvrnodm


### READ DATA ###

load("covid_spatial.RData")


### DEFINE BAYESIAN MODELS ###

#initial specifications
jags_model = function()
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

#candidate model 1
jags_model1 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #zip = number of ZIP codes
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  sp ~ dbeta(100, 3.02)                     #specificity of surveillance process to identify case
  #sn is defined by ZIP code
  
  #priors: these do not vary by ZIP code
  eta0 ~ dnorm(-3.5, 1/(0.5^2))            #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))           #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                   #precision of random component u
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)          #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i]    #candidate model 1: unstructured model
    #r[i] ~ dunif(0, 0.1)                       #candidate model 2: direct prior
    
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

#candidate model 2
jags_model2 = function()
{    
  ## data
  
  #y.o = observed count of COVID-19 cases in each ZIP code
  #t = population of the ZIP code
  #zip = number of ZIP codes
  
  ##priors: these do not vary by ZIP code
  
  #accuracy of surveillance process to identify case
  sp ~ dbeta(100, 3.02)                     #specificity of surveillance process to identify case
  #sn is defined by ZIP code
  
  #priors: these do not vary by ZIP code
  eta0 ~ dnorm(-3.5, 1/(0.5^2))            #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))           #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                   #precision of random component u
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)          #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence
    #logit(r[i]) <- eta0 + eta1*adi[i] + u[i]   #candidate model 1: unstructured model
    r[i] ~ dunif(0, 0.1)                        #candidate model 2: direct prior
    
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
#sa.nb = poly2nb(philly_sf, queen=T)

#create a binary spatial weights matrix from neighbors
#sa.wt = nb2mat(neighbours=sa.nb, style="B")

#line_data_observed = list("y.o"=philly_data$Positive, "t"=philly_data$Population, "adi"=as.numeric(philly_data$census_NDI), "zip"=nrow(philly_data), "d"=diag(rowSums(sa.wt)), "w"=sa.wt)
line_data_observed = list("y.o"=philly_data$Positive, "t"=philly_data$Population, "adi"=as.numeric(philly_data$census_NDI), "zip"=nrow(philly_data))

#run JAGS with benchmarking
time1 = Sys.time()
model1_posterior = jags.parallel(data=line_data_observed, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","u"), model.file=jags_model1, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
model2_posterior = jags.parallel(data=line_data_observed, parameters.to.save=c("r","y.t","sn","sp","err","eta0","eta1","u"), model.file=jags_model2, n.chains=5, n.thin=10, n.iter=50000, n.burnin=5000)
time2 = Sys.time()
time2-time1
rm(time1, time2)

#save
save.image("bayes_posterior.RData")


### POSTERIOR VISUALIZATIONS and DIAGNOSTICS ###

#read saved data
load("bayes_posterior.RData")

#obtain Gelman-Rubin convergence diagnostic (Rhat: potential scale reduction factor)
options(max.print=9999)
print(model1_posterior)
print(model2_posterior)

#diagnostics
MCMCtrace(model1_posterior, params=c("sn","sp","eta0","eta1","r","u","deviance"), wd="", filename="Appendix 3.pdf", open_pdf=F)
MCMCtrace(model2_posterior, params=c("sn","sp","eta0","eta1","r","u","deviance"), wd="", filename="Appendix 4.pdf", open_pdf=F)

#save model syntax file
writeLines(capture.output(getAnywhere(jags_model1)), file("bugs_model.txt"))
writeLines(capture.output(getAnywhere(jags_model2)), file("bugs_model.txt"))

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="posteriors.pdf", type="pdf", width=11, height=8.5)

#posterior prevalence plot (r) comparing 4 models
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(2,0.4))
par(new=T)
MCMCplot(model2_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05), col=adjustcolor(3,0.4))
par(new=T)
points(x=1:nrow(philly_data),y=(philly_data$Positive/philly_data$Population), col="black",pch=16)
legend("topright",legend=c("Model 1", "Model 2", "Observed"), pch=16, col=c(2,3,1), y.intersp=0.5)

#posterior count plot (y.t)
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","r","sn","sp","err","eta0","eta1","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Expected cases (95% CrI)", main="Posterior ZCTA Expected Cases (y.t)", ylim=c(0,2500))
points(x=1:nrow(philly_data),y=philly_data$Positive, col="red",pch=16)
legend("topright",legend=c("Expected","Observed"), pch=16, col=c("black","red"))

#posterior error plot (err)
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","eta0","eta1","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Expected cases (95% CrI)", main="Posterior ZCTA Missed Cases (err)", ylim=c(0,1000))

#posterior SN, SP plot
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","u"), ylab="Accuracy (95% CrI)", main="Posterior Sensitivity (sn) and Specificity (sp)", ylim=c(0.4,1))

#posterior logit prevalence (eta0, eta1)
MCMCplot(model1_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","u"), ylab="Logit true prevalance", main="Average True Prevalence (eta0) and ADI Effect (eta)", ylim=c(-4.1,1))

#posterior spatial (phi, tau)
#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","eta0","eta1"), main="Spatial Correlation (phi) and Precision (tau)", ylim=c(0,20))

#posterior non-spatial (tau)
#MCMCplot(covid_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","err","eta0","eta1","phi"), main="Precision (tau)", ylim=c(-10,10))

dev.off()


### EXTRACT POSTERIOR ESTIMATES ##

#prevalance
prev_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
  chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
  
  #add zip code specific posteriors to dataframe
  prev_posterior = data.frame(matrix(NA,ncol=zip,nrow=nrow(chain1)*5))
  col_mapping = colnames(prev_posterior)[order(names(prev_posterior))]
  for (i in 1:ncol(prev_posterior)) {
    
    #determine correct column (offset is based on first y.t column in mcmc output)
    col_index = which(col_mapping==colnames(prev_posterior)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
    
    #add data
    prev_posterior[,i] = c(chain1[,col_index]/philly_data$Population[i], chain2[,col_index]/philly_data$Population[i], chain3[,col_index]/philly_data$Population[i], chain4[,col_index]/philly_data$Population[i], chain5[,col_index]/philly_data$Population[i])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
  
  return(prev_posterior)
}

prev_model1 = prev_posterior(model1_posterior)
prev_model2 = prev_posterior(model2_posterior)

#sensitivity
sn_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
  chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
  
  #add zip code specific posteriors to dataframe
  sn_posterior = data.frame(matrix(NA,ncol=zip,nrow=nrow(chain1)*5))
  col_mapping = colnames(sn_posterior)[order(names(sn_posterior))]
  for (i in 1:ncol(sn_posterior)) {
    
    #determine correct column (offset is based on first sn column in mcmc output)
    col_index = which(col_mapping==colnames(sn_posterior)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
    
    #add data
    sn_posterior[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
  
  return(sn_posterior)
}

sn_model1 = sn_posterior(model1_posterior)
sn_model2 = sn_posterior(model2_posterior)

#missed counts
err_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
  chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
  
  #add zip code specific posteriors to dataframe
  err_posterior = data.frame(matrix(NA,ncol=zip,nrow=nrow(chain1)*5))
  col_mapping = colnames(err_posterior)[order(names(err_posterior))]
  for (i in 1:ncol(err_posterior)) {
    
    #determine correct column (offset is based on first sn column in mcmc output)
    col_index = which(col_mapping==colnames(err_posterior)[i]) + (which(colnames(chain1)=="err[1]") - 1)
    
    #add data
    err_posterior[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
  
  return(err_posterior)
}

err_model1 = err_posterior(model1_posterior)
err_model2 = err_posterior(model2_posterior)


### POSTERIOR MODEL COMPARISON PLOTS ###

##Prevalence

#create the boxplot dataframe
boxplot_df = data.frame("r"=c(c(prev_model1$X1, prev_model2$X1),
                              c(prev_model1$X2, prev_model2$X2),
                              c(prev_model1$X3, prev_model2$X3),
                              c(prev_model1$X4, prev_model2$X4),
                              c(prev_model1$X5, prev_model2$X5),
                              c(prev_model1$X6, prev_model2$X6),
                              c(prev_model1$X7, prev_model2$X7),
                              c(prev_model1$X8, prev_model2$X8),
                              c(prev_model1$X9, prev_model2$X9),
                              c(prev_model1$X10, prev_model2$X10),
                              c(prev_model1$X11, prev_model2$X11),
                              c(prev_model1$X12, prev_model2$X12),
                              c(prev_model1$X13, prev_model2$X13),
                              c(prev_model1$X14, prev_model2$X14),
                              c(prev_model1$X15, prev_model2$X15),
                              c(prev_model1$X16, prev_model2$X16),
                              c(prev_model1$X17, prev_model2$X17),
                              c(prev_model1$X18, prev_model2$X18),
                              c(prev_model1$X19, prev_model2$X19),
                              c(prev_model1$X20, prev_model2$X20),
                              c(prev_model1$X21, prev_model2$X21),
                              c(prev_model1$X22, prev_model2$X22),
                              c(prev_model1$X23, prev_model2$X23),
                              c(prev_model1$X24, prev_model2$X24),
                              c(prev_model1$X25, prev_model2$X25),
                              c(prev_model1$X26, prev_model2$X26),
                              c(prev_model1$X27, prev_model2$X27),
                              c(prev_model1$X28, prev_model2$X28),
                              c(prev_model1$X29, prev_model2$X29),
                              c(prev_model1$X30, prev_model2$X30),
                              c(prev_model1$X31, prev_model2$X31),
                              c(prev_model1$X32, prev_model2$X32),
                              c(prev_model1$X33, prev_model2$X33),
                              c(prev_model1$X34, prev_model2$X34),
                              c(prev_model1$X35, prev_model2$X35),
                              c(prev_model1$X36, prev_model2$X36),
                              c(prev_model1$X37, prev_model2$X37),
                              c(prev_model1$X38, prev_model2$X38),
                              c(prev_model1$X39, prev_model2$X39),
                              c(prev_model1$X40, prev_model2$X40),
                              c(prev_model1$X41, prev_model2$X41),
                              c(prev_model1$X42, prev_model2$X42),
                              c(prev_model1$X43, prev_model2$X43),
                              c(prev_model1$X44, prev_model2$X44),
                              c(prev_model1$X45, prev_model2$X45),
                              c(prev_model1$X46, prev_model2$X46),
                              c(prev_model1$X47, prev_model2$X47)),
                        "zcta"=rep(philly_data$ZCTA[1:zip], each=nrow(prev_model1)*2),
                        "model"=rep(rep(c("Model 1","Model 2"), each=nrow(prev_model1)), zip),
                        stringsAsFactors=F)

#x-axis zip codes
plot_labels = rep(philly_data$ZCTA, each=2)
plot_labels[seq(2,94, by=2)] = ""

#export @ 1600 x 600
boxplot(r ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", xlab="", ylab="Prevalence", col=c("#999999", "#FFFFFF"))
axis(side=1, labels=philly_data$ZCTA, at=seq(2,94, by=2), las=2, tick=F)
rect(seq(0.5,94,by=4), rep(-1,24), seq(2.5,94.5,by=4), 1, col=rgb(0.5,0.5,0.5,1/5), border=NA)
par(new=T)
boxplot(r ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", yaxt="n", xlab="", ylab="", col=c("#999999", "#FFFFFF"))
points(x=seq(1.5,94,by=2),y=(philly_data$Positive/philly_data$Population), col="black", pch=8, cex=1.5)
legend("topleft", legend=c("Model 1","Model 2", "Observed"), fill=c(c("#999999", "#FFFFFF"), NA), border=c("black","black",NA), y.intersp=0.6, cex=0.8, pch=c(NA,NA,8))

# #create the boxplot dataframe
# boxplot_df = data.frame("r"=c(unlist(prev_model3[, 1:zip])),
#                         "zcta"=rep(philly_data$ZCTA[1:zip], each=nrow(prev_model3)),
#                         stringsAsFactors=F)
# 
# boxplot(r ~ zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", xlab="", ylab="Prevalence")
# axis(side=1, labels=philly_data$ZCTA, at=1:zip, las=2, tick=F)
# points(x=1:zip,y=(philly_data$Positive/philly_data$Population), col="black", pch=8, cex=1)
# legend("topleft", legend=c("Posterior","Observed"), fill=c("#FFFFFF", NA), border=c("black",NA), y.intersp=0.6, cex=0.8, pch=c(NA,8))

##Sensitivity

#create the boxplot dataframe
boxplot_df = data.frame("sn"=c(c(sn_model1$X1, sn_model2$X1),
                               c(sn_model1$X2, sn_model2$X2),
                               c(sn_model1$X3, sn_model2$X3),
                               c(sn_model1$X4, sn_model2$X4),
                               c(sn_model1$X5, sn_model2$X5),
                               c(sn_model1$X6, sn_model2$X6),
                               c(sn_model1$X7, sn_model2$X7),
                               c(sn_model1$X8, sn_model2$X8),
                               c(sn_model1$X9, sn_model2$X9),
                               c(sn_model1$X10, sn_model2$X10),
                               c(sn_model1$X11, sn_model2$X11),
                               c(sn_model1$X12, sn_model2$X12),
                               c(sn_model1$X13, sn_model2$X13),
                               c(sn_model1$X14, sn_model2$X14),
                               c(sn_model1$X15, sn_model2$X15),
                               c(sn_model1$X16, sn_model2$X16),
                               c(sn_model1$X17, sn_model2$X17),
                               c(sn_model1$X18, sn_model2$X18),
                               c(sn_model1$X19, sn_model2$X19),
                               c(sn_model1$X20, sn_model2$X20),
                               c(sn_model1$X21, sn_model2$X21),
                               c(sn_model1$X22, sn_model2$X22),
                               c(sn_model1$X23, sn_model2$X23),
                               c(sn_model1$X24, sn_model2$X24),
                               c(sn_model1$X25, sn_model2$X25),
                               c(sn_model1$X26, sn_model2$X26),
                               c(sn_model1$X27, sn_model2$X27),
                               c(sn_model1$X28, sn_model2$X28),
                               c(sn_model1$X29, sn_model2$X29),
                               c(sn_model1$X30, sn_model2$X30),
                               c(sn_model1$X31, sn_model2$X31),
                               c(sn_model1$X32, sn_model2$X32),
                               c(sn_model1$X33, sn_model2$X33),
                               c(sn_model1$X34, sn_model2$X34),
                               c(sn_model1$X35, sn_model2$X35),
                               c(sn_model1$X36, sn_model2$X36),
                               c(sn_model1$X37, sn_model2$X37),
                               c(sn_model1$X38, sn_model2$X38),
                               c(sn_model1$X39, sn_model2$X39),
                               c(sn_model1$X40, sn_model2$X40),
                               c(sn_model1$X41, sn_model2$X41),
                               c(sn_model1$X42, sn_model2$X42),
                               c(sn_model1$X43, sn_model2$X43),
                               c(sn_model1$X44, sn_model2$X44),
                               c(sn_model1$X45, sn_model2$X45),
                               c(sn_model1$X46, sn_model2$X46),
                               c(sn_model1$X47, sn_model2$X47)),
                        "zcta"=rep(philly_data$ZCTA[1:zip], each=nrow(sn_model1)*2),
                        "model"=rep(rep(c("Model 1","Model 2"), each=nrow(sn_model1)), zip),
                        stringsAsFactors=F)

#x-axis zip codes
plot_labels = rep(philly_data$ZCTA, each=2)
plot_labels[seq(2,94, by=2)] = ""

#export @ 1600 x 600
boxplot(sn ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", xlab="", ylab="Sensitivity", col=c("#999999", "#FFFFFF"), ylim=c(0.1,1))
axis(side=1, labels=philly_data$ZCTA, at=seq(2,94, by=2), las=2, tick=F)
rect(seq(0.5,94,by=4), rep(-1,24), seq(2.5,94.5,by=4), 2, col=rgb(0.5,0.5,0.5,1/5), border=NA)
par(new=T)
boxplot(sn ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", yaxt="n", xlab="", ylab="", col=c("#999999", "#FFFFFF"), ylim=c(0.1,1))
legend("topleft", legend=c("Model 1","Model 2"), fill=c("#999999", "#FFFFFF"), border=c("black","black"), y.intersp=0.6, cex=0.8)

##Missed counts: err

#create the boxplot dataframe
boxplot_df = data.frame("err"=c(c(err_model1$X1, err_model2$X1),
                                c(err_model1$X2, err_model2$X2),
                                c(err_model1$X3, err_model2$X3),
                                c(err_model1$X4, err_model2$X4),
                                c(err_model1$X5, err_model2$X5),
                                c(err_model1$X6, err_model2$X6),
                                c(err_model1$X7, err_model2$X7),
                                c(err_model1$X8, err_model2$X8),
                                c(err_model1$X9, err_model2$X9),
                                c(err_model1$X10, err_model2$X10),
                                c(err_model1$X11, err_model2$X11),
                                c(err_model1$X12, err_model2$X12),
                                c(err_model1$X13, err_model2$X13),
                                c(err_model1$X14, err_model2$X14),
                                c(err_model1$X15, err_model2$X15),
                                c(err_model1$X16, err_model2$X16),
                                c(err_model1$X17, err_model2$X17),
                                c(err_model1$X18, err_model2$X18),
                                c(err_model1$X19, err_model2$X19),
                                c(err_model1$X20, err_model2$X20),
                                c(err_model1$X21, err_model2$X21),
                                c(err_model1$X22, err_model2$X22),
                                c(err_model1$X23, err_model2$X23),
                                c(err_model1$X24, err_model2$X24),
                                c(err_model1$X25, err_model2$X25),
                                c(err_model1$X26, err_model2$X26),
                                c(err_model1$X27, err_model2$X27),
                                c(err_model1$X28, err_model2$X28),
                                c(err_model1$X29, err_model2$X29),
                                c(err_model1$X30, err_model2$X30),
                                c(err_model1$X31, err_model2$X31),
                                c(err_model1$X32, err_model2$X32),
                                c(err_model1$X33, err_model2$X33),
                                c(err_model1$X34, err_model2$X34),
                                c(err_model1$X35, err_model2$X35),
                                c(err_model1$X36, err_model2$X36),
                                c(err_model1$X37, err_model2$X37),
                                c(err_model1$X38, err_model2$X38),
                                c(err_model1$X39, err_model2$X39),
                                c(err_model1$X40, err_model2$X40),
                                c(err_model1$X41, err_model2$X41),
                                c(err_model1$X42, err_model2$X42),
                                c(err_model1$X43, err_model2$X43),
                                c(err_model1$X44, err_model2$X44),
                                c(err_model1$X45, err_model2$X45),
                                c(err_model1$X46, err_model2$X46),
                                c(err_model1$X47, err_model2$X47)),
                        "zcta"=rep(philly_data$ZCTA[1:zip], each=nrow(err_model1)*2),
                        "model"=rep(rep(c("Model 1","Model 2"), each=nrow(err_model1)), zip),
                        stringsAsFactors=F)

#x-axis zip codes
plot_labels = rep(philly_data$ZCTA, each=2)
plot_labels[seq(2,94, by=2)] = ""

#export @ 1600 x 600
boxplot(err ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", xlab="", ylab="Missed cases", col=c("#999999", "#FFFFFF"), ylim=c(-250,1500))
axis(side=1, labels=philly_data$ZCTA, at=seq(2,94, by=2), las=2, tick=F)
rect(seq(0.5,94,by=4), rep(-1000,24), seq(2.5,94.5,by=4), 10000, col=rgb(0.5,0.5,0.5,1/5), border=NA)
par(new=T)
boxplot(err ~ model + zcta, data=boxplot_df, outline=F, whisklty=1, xaxt="n", yaxt="n", xlab="", ylab="", col=c("#999999", "#FFFFFF"), ylim=c(-250,1500))
legend("topleft", legend=c("Model 1","Model 2"), fill=c("#999999", "#FFFFFF"), border=c("black","black"), y.intersp=0.6, cex=0.8)


### POSTERIOR ESTIMATES for MANUSCRIPT ###

inference_posterior = model1_posterior
inference_prev = prev_model1
plot_title = "A)"

#prevalence
prev_post = data.frame(inference_posterior$BUGSoutput$summary[grep("r",substr(row.names(inference_posterior$BUGSoutput$summary),1,1)),], stringsAsFactors=F)
#mean(prev_post$X50)
#mean(prev_post$X2.5.)
#mean(prev_post$X97.5.)

#citywide
inference_prev$prev = (inference_prev$X1*philly_data$Population[1] +
                         inference_prev$X2*philly_data$Population[2] +
                         inference_prev$X3*philly_data$Population[3] +
                         inference_prev$X4*philly_data$Population[4] +
                         inference_prev$X5*philly_data$Population[5] +
                         inference_prev$X6*philly_data$Population[6] +
                         inference_prev$X7*philly_data$Population[7] +
                         inference_prev$X8*philly_data$Population[8] +
                         inference_prev$X9*philly_data$Population[9] +
                         inference_prev$X10*philly_data$Population[10] +
                         inference_prev$X11*philly_data$Population[11] +
                         inference_prev$X12*philly_data$Population[12] +
                         inference_prev$X13*philly_data$Population[13] +
                         inference_prev$X14*philly_data$Population[14] +
                         inference_prev$X15*philly_data$Population[15] +
                         inference_prev$X16*philly_data$Population[16] +
                         inference_prev$X17*philly_data$Population[17] +
                         inference_prev$X18*philly_data$Population[18] +
                         inference_prev$X19*philly_data$Population[19] +
                         inference_prev$X20*philly_data$Population[20] +
                         inference_prev$X21*philly_data$Population[21] +
                         inference_prev$X22*philly_data$Population[22] +
                         inference_prev$X23*philly_data$Population[23] +
                         inference_prev$X24*philly_data$Population[24] +
                         inference_prev$X25*philly_data$Population[25] +
                         inference_prev$X26*philly_data$Population[26] +
                         inference_prev$X27*philly_data$Population[27] +
                         inference_prev$X28*philly_data$Population[28] +
                         inference_prev$X29*philly_data$Population[29] +
                         inference_prev$X30*philly_data$Population[30] +
                         inference_prev$X31*philly_data$Population[31] +
                         inference_prev$X32*philly_data$Population[32] +
                         inference_prev$X33*philly_data$Population[33] +
                         inference_prev$X34*philly_data$Population[34] +
                         inference_prev$X35*philly_data$Population[35] +
                         inference_prev$X36*philly_data$Population[36] +
                         inference_prev$X37*philly_data$Population[37] +
                         inference_prev$X38*philly_data$Population[38] +
                         inference_prev$X39*philly_data$Population[39] +
                         inference_prev$X40*philly_data$Population[40] +
                         inference_prev$X41*philly_data$Population[41] +
                         inference_prev$X42*philly_data$Population[42] +
                         inference_prev$X43*philly_data$Population[43] +
                         inference_prev$X44*philly_data$Population[44] +
                         inference_prev$X45*philly_data$Population[45] +
                         inference_prev$X46*philly_data$Population[46] +
                         inference_prev$X47*philly_data$Population[47]) / sum(philly_data$Population)
quantile(inference_prev$prev, probs=c(0.025,0.5,0.975))

prev_post[prev_post$X50==min(prev_post$X50), ]
philly_data$ZCTA[which(prev_post$X50==min(prev_post$X50))]
prev_post[prev_post$X50==max(prev_post$X50), ]
philly_data$ZCTA[which(prev_post$X50==max(prev_post$X50))]

#MCMCplot(inference_posterior, horiz=F, excl=c("deviance","y.t","sn","sp","err","eta0","eta1","phi","tau"), labels=philly_data$ZCTA, sz_labels=1, ylab="Prevalence (95% CrI)", ylim=c(0,0.05))
#points(x=1:nrow(philly_data),y=(philly_data$Positive/philly_data$Population), pch=1)
#legend("topright",legend=c("Expected","Observed"), pch=c(16,1))

#SN
sn_post = data.frame(inference_posterior$BUGSoutput$summary[grep("sn",row.names(inference_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
mean(sn_post$X50)
mean(sn_post$X2.5.)
mean(sn_post$X97.5.)

sn_post[sn_post$X50==min(sn_post$X50), ]
philly_data$ZCTA[which(sn_post$X50==min(sn_post$X50))]
sn_post[sn_post$X50==max(sn_post$X50), ]
philly_data$ZCTA[which(sn_post$X50==max(sn_post$X50))]

#MCMCplot(inference_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","sp","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Sensitivity (95% CrI)", ylim=c(0.2,0.9), main=plot_title)

#SP
#sp_post = data.frame(inference_posterior$BUGSoutput$summary[grep("sp",row.names(inference_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
#mean(sp_post$X50)
#mean(sp_post$X2.5.)
#mean(sp_post$X97.5.)

data.frame(inference_posterior$BUGSoutput$summary[grep("sp",row.names(inference_posterior$BUGSoutput$summary)),], stringsAsFactors=F)

#sp_post[sp_post$X50==min(sp_post$X50), ]
#sp_post[sp_post$X50==max(sp_post$X50), ]
#philly_data$ZCTA[which(sn_post$X50==min(sn_post$X50))]
#philly_data$ZCTA[which(sn_post$X50==max(sn_post$X50))]

#MCMCplot(inference_posterior, horiz=F, excl=c("deviance","r","err","y.t","eta0","eta1","phi","tau","sn"), labels=philly_data$ZCTA, sz_labels=1, ylab="Accuracy (95% CrI)", ylim=c(0.9,1))

#measurement error
err_post = data.frame(inference_posterior$BUGSoutput$summary[grep("err",row.names(inference_posterior$BUGSoutput$summary)),], stringsAsFactors=F)
#MCMCplot(inference_posterior, horiz=F, excl=c("deviance","r","sn","sp","y.t","eta0","eta1","u"), labels=philly_data$ZCTA, sz_labels=1, ylab="Missed cases (95% CrI)", ylim=c(-100,1500), main=plot_title)

err_post[err_post$X50==max(err_post$X50), ]
philly_data$ZCTA[which(err_post$X50==max(err_post$X50))]
philly_data$Positive[which(err_post$X50==max(err_post$X50))]/philly_data$Population[which(err_post$X50==max(err_post$X50))]
prev_post[which(err_post$X50==max(err_post$X50)), ]


# ### PRIOR and POSTERIOR PLOTS ##
# 
# par(mfrow=c(3,2))
# 
# #phi
# plot(NULL, main="A) phi", xlab="Value", ylab="Density", yaxt='n', xlim=c(-0.1,1.1), ylim=c(0,1.5))
# 
# #obtain posterior
# phi = c(as.data.frame(as.mcmc(covid_posterior)[[1]])$phi,
#         as.data.frame(as.mcmc(covid_posterior)[[2]])$phi,
#         as.data.frame(as.mcmc(covid_posterior)[[3]])$phi,
#         as.data.frame(as.mcmc(covid_posterior)[[1]])$phi,
#         as.data.frame(as.mcmc(covid_posterior)[[1]])$phi)
# 
# #add lines
# lines(density(phi), lwd=2, col="gray")
# lines(density(runif(covid_posterior$BUGSoutput$n.sims, 0, 1)), lwd=2, col="black")
# #legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
# rm(phi)
# 
# #sp
# plot(NULL, main="B) specificity", xlab="Value", ylab="Density", yaxt='n', xlim=c(0.95,1), ylim=c(0,500))
# 
# #obtain posterior
# sp = c(as.data.frame(as.mcmc(covid_posterior)[[1]])$sp,
#         as.data.frame(as.mcmc(covid_posterior)[[2]])$sp,
#         as.data.frame(as.mcmc(covid_posterior)[[3]])$sp,
#         as.data.frame(as.mcmc(covid_posterior)[[1]])$sp,
#         as.data.frame(as.mcmc(covid_posterior)[[1]])$sp)
# 
# #add lines
# lines(density(sp), lwd=2, col="gray")
# lines(density(rbeta(covid_posterior$BUGSoutput$n.sims, 100, 3.02)), lwd=2, col="black")
# #legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
# rm(sp)
# 
# #sn
# plot(NULL, main="C) sensitivity", xlab="Value", ylab="Density", yaxt='n', xlim=c(0.2,0.9), ylim=c(0,5))
# 
# #obtain posterior and add lines
# chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
# chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
# chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
# chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
# chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
# 
# #add zip code specific posteriors to dataframe
# sn = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
# col_mapping = colnames(sn)[order(names(sn))]
# for (i in 1:ncol(sn)) {
#   
#   #determine correct column (offset is based on first sn column in mcmc output)
#   col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
#   
#   #add data
#   sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
# }
# rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
# 
# #add lines
# for (i in 1:zip) {
#   lines(density(sn[,i]), lwd=1, col="gray")
# }
# rm(i, sn)
# 
# lines(density(rbeta(covid_posterior$BUGSoutput$n.sims, 14.022, 9.681)), lwd=2, col="black")
# #legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
# 
# # #u
# # plot(NULL, main="D) u", xlab="Value", ylab="Density", yaxt='n', xlim=c(-20,20), ylim=c(0,2))
# # 
# # #obtain posterior and add lines
# # chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
# # chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
# # chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
# # chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
# # chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
# # 
# # #add zip code specific posteriors to dataframe
# # u = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
# # col_mapping = colnames(u)[order(names(u))]
# # for (i in 1:ncol(u)) {
# #   
# #   #determine correct column (offset is based on first u column in mcmc output)
# #   col_index = which(col_mapping==colnames(u)[i]) + (which(colnames(chain1)=="u[1]") - 1)
# #   
# #   #add data
# #   u[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
# # }
# # rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
# # 
# # #add lines
# # for (i in 1:zip) {
# #   lines(density(u[,i]), lwd=1, col="gray")
# # }
# # rm(i, u)
# # 
# # lines(density(rnorm(covid_posterior$BUGSoutput$n.sims, 0, median(rgamma(covid_posterior$BUGSoutput$n.sims, 1, 0.1)))), lwd=2, col="black")
# # #legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))
# 
# #v
# plot(NULL, main="E) v", xlab="Value", ylab="Density", yaxt='n', xlim=c(-20,20), ylim=c(0,2))
# 
# #obtain posterior and add lines
# chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
# chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
# chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
# chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
# chain5 = as.data.frame(as.mcmc(covid_posterior)[[5]])
# 
# #add zip code specific posteriors to dataframe
# v = data.frame(matrix(NA,ncol=zip,nrow=covid_posterior$BUGSoutput$n.sims))
# col_mapping = colnames(v)[order(names(v))]
# for (i in 1:ncol(v)) {
#   
#   #determine correct column (offset is based on first v column in mcmc output)
#   col_index = which(col_mapping==colnames(v)[i]) + (which(colnames(chain1)=="v[1]") - 1)
#   
#   #add data
#   v[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
# }
# rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)
# 
# #add lines
# for (i in 1:zip) {
#   lines(density(v[,i]), lwd=1, col="gray")
# }
# rm(i, v)
# 
# lines(density(mvrnorm(covid_posterior$BUGSoutput$n.sims, mu=rep(0,zip), Sigma=median(rgamma(covid_posterior$BUGSoutput$n.sims, 1, 0.1)) * (d - median(runif(covid_posterior$BUGSoutput$n.sims, 0, 1))*w))), lwd=2, col="black")
# #legend("topleft", legend=c("Prior","Posterior"), lty=1, lwd=2, col=c("black","gray"))


# ### SENSITIVITY ANALYSIS CORRELATION MATRIX ###
# 
# plot(colMedians(as.matrix(prev_model3_sn)), colMedians(as.matrix(prev_model3)), pch=16, xlim=c(0.01,0.06), ylim=c(0,0.06), xlab="", ylab="")
# text(0.0125,0.05, expression(paste(rho)))
# text(0.018, 0.05, paste("=","",round(cor(colMedians(as.matrix(prev_model3)), colMedians(as.matrix(prev_model3_sn))), 2)))


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
#moran.mc(prev_post$X50., listw=sa.wt, nsim=10000)

#spatial correlogram, shows autocorrelation by order of neighbor (1=neighbor, 2=neighbor's neighbor, etc)
#plot(sp.correlogram(neighbours=sa.nb,var=prev_post$X50.,order=4,method="I",style="B",zero.policy=T), main="")

#choropleth map: see https://edzer.github.io/sp/ for helpful tips
philly_sf_joined$prev_post = prev_post$X50.
sp_arrow = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(2730000,220000), scale = 10000)
sp_scale = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(2725000,215000), scale = 21120, fill=c("transparent","black"))
sp_scale_text1 = list("sp.text", c(2726000,211000), "0")
sp_scale_text2 = list("sp.text", c(2742000,211000), "4 mi")
spplot(as_Spatial(philly_sf_joined), "prev_post", cuts=7, col.regions=brewer.pal(8, "Reds"), sp.layout=list(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2), main=plot_title)
rm(sp_arrow,sp_scale,sp_scale_text1,sp_scale_text2)

# #cartogram map
# philly_sf_joined$err = err_post$X50. / philly_sf_joined$Population
# carto = cartogram_cont(philly_sf_joined, "err", itermax=10)
# 
# par(mar=c(0.1,0.1,1,0.1))
# 
# #referent
# plot(philly_sf_joined$geometry, main="A)")
# text(t(sapply(slot(as_Spatial(philly_sf_joined), "polygons"), function(i) slot(i, "labpt"))), cex=0.6, labels=philly_sf_joined$ZCTA5CE10)
# 
# plot(carto$geometry, main="B)")
# plot(carto$geometry, main="C)")
