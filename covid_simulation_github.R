#################
# COVID-19 Bayesian simulation study
# Citation: Goldstein ND, Wheeler DC, Gustafson P, Burstyn I. A Bayesian Approach to Improving Spatial Estimates of Prevalence of COVID-19 After Accounting for Misclassification Bias in Surveillance Data in Philadelphia, PA. Manuscript in preparation.
# 9/30/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library("R2jags") #JAGS interface
library("MCMCvis") #visualizations


### READ DATA ###

load("covid_spatial.RData")


### CREATE SIMULATION DATA ###

sim_data = philly_data

#true prevalence is conditioned on ADI, so specify logistic model

#intercept: base prevalance when ADI=0 on log scale: 1%, 3%, 5%
beta0_1per=rnorm(length(sim_data$census_NDI), -4.5, 0.1)
beta0_3per=rnorm(length(sim_data$census_NDI), -3.5, 0.1)
beta0_5per=rnorm(length(sim_data$census_NDI), -3, 0.1)

#effect of ADI on prev
beta1=rnorm(length(sim_data$census_NDI), 0.04, 0.04)

#true prevalence by ZIP code conditioned on ADI
true_prev_1per = 1/(1+exp(-(beta0_1per+beta1*sim_data$census_NDI)))
true_prev_3per = 1/(1+exp(-(beta0_3per+beta1*sim_data$census_NDI)))
true_prev_5per = 1/(1+exp(-(beta0_5per+beta1*sim_data$census_NDI)))
rm(beta0_1per,beta0_3per,beta0_5per,beta1)

#100*quantile(true_prev_1per, prob=c(0.025, .5, .975))
#plot(density(true_prev_1per), main="average 1%", xlim=c(0.01, 0.07), xlab="true prevalence per zip code given ADI")

#100*quantile(true_prev_3per, prob=c(0.025, .5, .975))
#plot(density(true_prev_3per), main="average 3%", xlim=c(0.01, 0.07), xlab="true prevalence per zip code given ADI")

#100*quantile(true_prev_5per, prob=c(0.025, .5, .975))
#plot(density(true_prev_5per), main="average 5%", xlim=c(0.01, 0.07), xlab="true prevalence per zip code given ADI")

#accuracy of surveillance
SN_30per = rbeta(nrow(sim_data), 3, 7)
SN_50per = rbeta(nrow(sim_data), 5, 5)
SN_70per = rbeta(nrow(sim_data), 7, 3)
SP = 0.99

#hist(SN_30per, breaks="fd")
#hist(SN_50per, breaks="fd")
#hist(SN_70per, breaks="fd")

#scaling factor to ensure Y/Y* sufficiently large
sim_data$Population = sim_data$Population * 1

#generate new ZIP code case counts with specific combinations of parameters
sim_data$Positive_0.3_0.01 = rbinom(nrow(sim_data), sim_data$Population, (SN_30per*true_prev_1per + (1 - SP)*(1-true_prev_1per)))
sim_data$Positive_0.5_0.01 = rbinom(nrow(sim_data), sim_data$Population, (SN_50per*true_prev_1per + (1 - SP)*(1-true_prev_1per)))
sim_data$Positive_0.7_0.01 = rbinom(nrow(sim_data), sim_data$Population, (SN_70per*true_prev_1per + (1 - SP)*(1-true_prev_1per)))
sim_data$Positive_0.3_0.03 = rbinom(nrow(sim_data), sim_data$Population, (SN_30per*true_prev_3per + (1 - SP)*(1-true_prev_3per)))
sim_data$Positive_0.5_0.03 = rbinom(nrow(sim_data), sim_data$Population, (SN_50per*true_prev_3per + (1 - SP)*(1-true_prev_3per)))
sim_data$Positive_0.7_0.03 = rbinom(nrow(sim_data), sim_data$Population, (SN_70per*true_prev_3per + (1 - SP)*(1-true_prev_3per)))
sim_data$Positive_0.3_0.05 = rbinom(nrow(sim_data), sim_data$Population, (SN_30per*true_prev_5per + (1 - SP)*(1-true_prev_5per)))
sim_data$Positive_0.5_0.05 = rbinom(nrow(sim_data), sim_data$Population, (SN_50per*true_prev_5per + (1 - SP)*(1-true_prev_5per)))
sim_data$Positive_0.7_0.05 = rbinom(nrow(sim_data), sim_data$Population, (SN_70per*true_prev_5per + (1 - SP)*(1-true_prev_5per)))

#verify data generation
summary(sim_data$Positive_0.3_0.01/sim_data$Population)
summary(sim_data$Positive_0.5_0.01/sim_data$Population)
summary(sim_data$Positive_0.7_0.01/sim_data$Population)
summary(sim_data$Positive_0.3_0.03/sim_data$Population)
summary(sim_data$Positive_0.5_0.03/sim_data$Population)
summary(sim_data$Positive_0.7_0.03/sim_data$Population)
summary(sim_data$Positive_0.3_0.05/sim_data$Population)
summary(sim_data$Positive_0.5_0.05/sim_data$Population)
summary(sim_data$Positive_0.7_0.05/sim_data$Population)

sum(sim_data$Population)

#save
save.image("bayes_simulation.RData")


### DEFINE BAYESIAN MODELS ###

jags_model = function()
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
  #eta0 ~ dnorm(-3.5, 1/(0.5^2))            #average true prevalence (logit scale)
  eta0 ~ dnorm(0,0.368)                    #average true prevalence (logit scale)
  eta1 ~ dnorm(0.04, 1/(0.04^2))           #effect of ADI on true prevalence (logit scale)
  tau_u ~ dgamma(1, 0.1)                   #precision of random component u
  
  #construct model for each ZIP code
  for (i in 1:zip) {
    
    ## priors: these vary by ZIP code
    
    sn[i] ~ dbeta(14.022, 9.681)            #sensitivity of surveillance process to identify case
    #sn[i] ~ dunif(0.25, 0.75)              #sensitivity of surveillance process to identify case
    u[i] ~ dnorm(0, tau_u)                 #unstructured (non-spatial) random effect
    
    ## models
    
    # true prevalence
    logit(r[i]) <- eta0 + eta1*adi[i] + u[i]         #model 2: unstructured model
    #r[i] ~ dunif(0, 0.1)
    
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

#read saved data
load("bayes_simulation.RData")

line_data_sim_0.3_0.01 = list("y.o"=sim_data$Positive_0.3_0.01, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.5_0.01 = list("y.o"=sim_data$Positive_0.5_0.01, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.7_0.01 = list("y.o"=sim_data$Positive_0.7_0.01, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.3_0.03 = list("y.o"=sim_data$Positive_0.3_0.03, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.5_0.03 = list("y.o"=sim_data$Positive_0.5_0.03, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.7_0.03 = list("y.o"=sim_data$Positive_0.7_0.03, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.3_0.05 = list("y.o"=sim_data$Positive_0.3_0.05, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.5_0.05 = list("y.o"=sim_data$Positive_0.5_0.05, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))
line_data_sim_0.7_0.05 = list("y.o"=sim_data$Positive_0.7_0.05, "t"=sim_data$Population, "adi"=as.numeric(sim_data$census_NDI), "zip"=nrow(sim_data))

#run JAGS with benchmarking
time1 = Sys.time()
model_posterior_sim_0.3_0.01 = jags.parallel(data=line_data_sim_0.3_0.01, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.5_0.01 = jags.parallel(data=line_data_sim_0.5_0.01, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.7_0.01 = jags.parallel(data=line_data_sim_0.7_0.01, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.3_0.03 = jags.parallel(data=line_data_sim_0.3_0.03, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.5_0.03 = jags.parallel(data=line_data_sim_0.5_0.03, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.7_0.03 = jags.parallel(data=line_data_sim_0.7_0.03, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.3_0.05 = jags.parallel(data=line_data_sim_0.3_0.05, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.5_0.05 = jags.parallel(data=line_data_sim_0.5_0.05, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
model_posterior_sim_0.7_0.05 = jags.parallel(data=line_data_sim_0.7_0.05, parameters.to.save=c("r","y.t","sn","sp","err","eta0"), model.file=jags_model, n.chains=5, n.thin=10, n.iter=10000, n.burnin=1000)
time2 = Sys.time()
time2-time1
rm(time1, time2)

#save
save.image("Sen_posterior.RData")


### POSTERIOR VISUALIZATIONS and DIAGNOSTICS ###

#obtain Gelman-Rubin convergence diagnostic (Rhat: potential scale reduction factor)
options(max.print=9999)
print(model_posterior_sim_0.3_0.01)

#diagnostics
MCMCtrace(model_posterior_sim_0.3_0.01, params=c("sn","sp","r","deviance"), wd="", filename="Sen_diagnostics.pdf", open_pdf=F)


### SENSITIVITY ANALYSIS SIMULATION PLOTS: SENSITIVITY ##

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="Sensitivity.pdf", type="pdf", width=11, height=8.5)

par(mfrow=c(3,3))

#model_posterior_sim_0.3_0.01
plot(NULL, main="Sensitivity=0.3, Prevalence=0.01", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.01$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_30per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.01
plot(NULL, main="Sensitivity=0.5, Prevalence=0.01", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.01$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_50per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.01
plot(NULL, main="Sensitivity=0.7, Prevalence=0.01", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.01$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_70per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.03
plot(NULL, main="Sensitivity=0.3, Prevalence=0.03", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.03$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_30per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.03
plot(NULL, main="Sensitivity=0.5, Prevalence=0.03", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.03$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_50per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.03
plot(NULL, main="Sensitivity=0.7, Prevalence=0.03", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.03$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_70per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.05
plot(NULL, main="Sensitivity=0.3, Prevalence=0.05", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.05$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_30per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.05
plot(NULL, main="Sensitivity=0.5, Prevalence=0.05", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.05$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_50per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.05
plot(NULL, main="Sensitivity=0.7, Prevalence=0.05", xlab="Sensitivity", ylab="Density", yaxt='n', xlim=c(0,1), ylim=c(0,5))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[5]])

#add zip code specific posteriors to dataframe
sn = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.05$BUGSoutput$n.sims))
col_mapping = colnames(sn)[order(names(sn))]
for (i in 1:ncol(sn)) {
  
  #determine correct column (offset is based on first sn column in mcmc output)
  col_index = which(col_mapping==colnames(sn)[i]) + (which(colnames(chain1)=="sn[1]") - 1)
  
  #add data
  sn[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index], chain5[,col_index])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(sn[,i]), lwd=1, col="gray")
}
rm(i, sn)

lines(density(SN_70per), lwd=2, lty=2)
lines(density(rbeta(10000, 14.022, 9.681)), lwd=2, col="black")
legend(x=0.7, y=5, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

dev.off()


### SENSITIVITY ANALYSIS SIMULATION PLOTS: SPECIFICITY ##

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="Specificity.pdf", type="pdf", width=11, height=8.5)

par(mfrow=c(3,3))

#model_posterior_sim_0.3_0.01
plot(NULL, main="Sensitivity=0.3, Prevalence=0.01", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,300))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.01
plot(NULL, main="Sensitivity=0.5, Prevalence=0.01", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,300))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.01
plot(NULL, main="Sensitivity=0.7, Prevalence=0.01", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,300))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.03
plot(NULL, main="Sensitivity=0.3, Prevalence=0.03", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,200))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.03
plot(NULL, main="Sensitivity=0.5, Prevalence=0.03", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,200))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.03
plot(NULL, main="Sensitivity=0.7, Prevalence=0.03", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,200))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.05
plot(NULL, main="Sensitivity=0.3, Prevalence=0.05", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.05
plot(NULL, main="Sensitivity=0.5, Prevalence=0.05", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.05
plot(NULL, main="Sensitivity=0.7, Prevalence=0.05", xlab="Specificity", ylab="Density", yaxt='n', xlim=c(0.5,1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[5]])

#add lines
lines(density(c(chain1$sp,chain2$sp,chain3$sp,chain4$sp,chain5$sp)), lwd=1, col="gray")
rm(chain1,chain2,chain3,chain4,chain5)

#add an empty plot to reset scaling
par(new=T)
plot(NULL,type='n',axes=FALSE,ann=FALSE, xlim=c(0.5,1), ylim=c(0,50))

abline(v=SP, lwd=2, lty=2)
lines(density(rbeta(10000, 100, 3.02)), lwd=2, col="black")
legend(x=0.5, y=45, legend=c("True","Prior", "Posterior"), lty=c(2,1,1), lwd=c(2,2,1), col=c("black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

dev.off()


### SENSITIVITY ANALYSIS SIMULATION PLOTS: PREVALENCE ##

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="Prevalence.pdf", type="pdf", width=11, height=8.5)

par(mfrow=c(3,3))

#model_posterior_sim_0.3_0.01
plot(NULL, main="Sensitivity=0.3, Prevalence=0.01", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,250))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_1per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.3_0.01/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=250, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.01
plot(NULL, main="Sensitivity=0.5, Prevalence=0.01", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,250))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_1per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.5_0.01/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=250, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.01
plot(NULL, main="Sensitivity=0.7, Prevalence=0.01", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,250))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_1per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.7_0.01/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=250, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.03
plot(NULL, main="Sensitivity=0.3, Prevalence=0.03", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_3per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.3_0.03/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.03
plot(NULL, main="Sensitivity=0.5, Prevalence=0.03", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_3per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.5_0.03/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.03
plot(NULL, main="Sensitivity=0.7, Prevalence=0.03", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_3per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.7_0.03/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.3_0.05
plot(NULL, main="Sensitivity=0.3, Prevalence=0.05", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_5per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.3_0.05/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.5_0.05
plot(NULL, main="Sensitivity=0.5, Prevalence=0.05", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_5per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.5_0.05/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

#model_posterior_sim_0.7_0.05
plot(NULL, main="Sensitivity=0.7, Prevalence=0.05", xlab="Prevalence", ylab="Density", yaxt='n', xlim=c(0,0.1), ylim=c(0,150))

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#add lines
for (i in 1:zip) {
  lines(density(prev[,i]), lwd=1, col="gray")
}
rm(i, prev)

lines(density(true_prev_5per), lwd=2, lty=2)
lines(density(sim_data$Positive_0.7_0.05/sim_data$Population), lwd=2, lty=3)
lines(density(runif(10000, 0, 0.1)), lwd=2, col="black")
legend(x=0.07, y=150, legend=c("True", "Observed", "Prior", "Posterior"), lty=c(2,3,1,1), lwd=c(2,2,2,1), col=c("black","black","black","gray"), cex=0.6, y.intersp=0.8, bty="n")

dev.off()


### SENSITIVITY ANALYSIS SIMULATION PLOTS: PREVALENCE by TRUE ##

#close any display devices if open
if (length(dev.list())) dev.off()
quartz(file="Prevalence by True.pdf", type="pdf", width=11, height=8.5)

par(mfrow=c(3,3))

#model_posterior_sim_0.3_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.3, Prevalence=0.01", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_1per),max(true_prev_1per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=(true_prev_1per[i]), y=(median(prev[,i])), pch=8)
  arrows(x0=(true_prev_1per[i]), y0=(quantile(prev[,i],probs=0.025)), x1=(true_prev_1per[i]), y1=(quantile(prev[,i],probs=0.975)), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.5_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.5, Prevalence=0.01", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_1per),max(true_prev_1per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_1per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_1per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_1per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.7_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.7, Prevalence=0.01", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_1per),max(true_prev_1per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_1per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_1per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_1per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.3_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.3, Prevalence=0.03", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_3per),max(true_prev_3per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_3per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_3per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_3per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.5_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.5, Prevalence=0.03", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_3per),max(true_prev_3per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_3per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_3per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_3per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.7_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.7, Prevalence=0.03", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_3per),max(true_prev_3per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_3per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_3per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_3per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.3_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.3, Prevalence=0.05", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_5per),max(true_prev_5per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_5per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_5per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_5per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.5_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.5, Prevalence=0.05", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_5per),max(true_prev_5per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_5per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_5per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_5per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

#model_posterior_sim_0.7_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#determine bounds for plot
ylim_min = 100
ylim_max = 0
for (i in 1:ncol(prev)) {
  if (ylim_min>quantile(prev[,i],probs=0.025)) ylim_min = quantile(prev[,i],probs=0.025)
  if (ylim_max<quantile(prev[,i],probs=0.975)) ylim_max = quantile(prev[,i],probs=0.975)
  
}
rm(i)

plot(NULL, main="Sensitivity=0.7, Prevalence=0.05", xlab="True Prevalence", ylab="Posterior Prevalence", xlim=c(min(true_prev_5per),max(true_prev_5per)), ylim=c(ylim_min,ylim_max))
rm(ylim_min, ylim_max)

#add median points and credible intervals
for (i in 1:zip) {
  
  points(x=true_prev_5per[i], y=median(prev[,i]), pch=8)
  arrows(x0=true_prev_5per[i], y0=quantile(prev[,i],probs=0.025), x1=true_prev_5per[i], y1=quantile(prev[,i],probs=0.975), angle=90, code=3, length=0.1, lwd=0.5)
}
rm(i, prev)

#add identity and prior line
abline(a=0, b=1, lty=1, lwd=2)
abline(h=0.05, lty=2, lwd=2)

dev.off()


### MEAN SQUARED ERROR of PREVALENCE ###

mse_allmodels = NA

#model_posterior_sim_0.3_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.3_0.01[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.5_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.5_0.01[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.7_0.01

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.01)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.01$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.7_0.01[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.3_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.3_0.03[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.5_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.5_0.03[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.7_0.03

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.03)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.03$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.7_0.03[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.3_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.3_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.3_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.3_0.05[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.5_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.5_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.5_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.5_0.05[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

#model_posterior_sim_0.7_0.05

#obtain posterior and add lines
chain1 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[1]])
chain2 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[2]])
chain3 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[3]])
chain4 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[4]])
chain5 = as.data.frame(as.mcmc(model_posterior_sim_0.7_0.05)[[5]])

#add zip code specific posteriors to dataframe
prev = data.frame(matrix(NA,ncol=zip,nrow=model_posterior_sim_0.7_0.05$BUGSoutput$n.sims))
col_mapping = colnames(prev)[order(names(prev))]
for (i in 1:ncol(prev)) {
  
  #determine correct column (offset is based on first y.t column in mcmc output)
  col_index = which(col_mapping==colnames(prev)[i]) + (which(colnames(chain1)=="y.t[1]") - 1)
  
  #add data
  prev[,i] = c(chain1[,col_index]/sim_data$Population[i], chain2[,col_index]/sim_data$Population[i], chain3[,col_index]/sim_data$Population[i], chain4[,col_index]/sim_data$Population[i], chain5[,col_index]/sim_data$Population[i])
}
rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4,chain5)

#compute mean square error by ZIP
mse = NA
for (i in 1:ncol(prev)) {
  mse = c(mse, ((sim_data$Positive_0.7_0.05[i]/sim_data$Population[i]) - as.numeric(quantile(prev[,i],probs=0.5)))^2)
}
rm(i)

#mse
mse_allmodels = c(mse_allmodels, sum(mse, na.rm=T))

mse_allmodels
