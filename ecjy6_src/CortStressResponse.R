#Appendix. Main R syntax document (PSM_CortStressResponse.R)

library(PSM)
library(nlme)

PKdata <- groupedData(DV ~ TAD | ID, read.csv("CortStressResponse.csv")[,-1])

#format data for PSM modeling functions
N <- length(unique(PKdata$ID))
PKdat <- vector(mode="list",length=N)
for(i in 1:N){
  PKdat[[i]]$Time <- with(PKdata[as.numeric(PKdata$ID) == i,], TAD)
  PKdat[[i]]$Y <- with(PKdata[as.numeric(PKdata$ID) == i,], matrix(DV, nrow=1))
  if(PKdata$FLAG[as.numeric(PKdata$ID) == i][1] == 2){
    PKdat[[i]]$Time <- c(PKdat[[i]]$Time,0,1)
    PKdat[[i]]$Y <- matrix(c(PKdat[[i]]$Y,NA,NA)[order(PKdat[[i]]$Time)], nrow=1)
    PKdat[[i]]$Time <- sort(PKdat[[i]]$Time)
  }
  if(PKdata$FLAG[as.numeric(PKdata$ID) == i][1] == 1){
    PKdat[[i]]$Time <- c(PKdat[[i]]$Time,1)
    PKdat[[i]]$Y <- matrix(c(PKdat[[i]]$Y,NA)[order(PKdat[[i]]$Time)], nrow=1)
    PKdat[[i]]$Time <- sort(PKdat[[i]]$Time)
  }
  PKdat[[i]]$U <- matrix(c(ifelse(PKdat[[i]]$Time != 0, 0, 1),
                           rep(1, length(PKdat[[i]]$Time))
  ), nrow=2, byrow=T)
  
  PKdat[[i]]$covar <- c(STUDY = PKdata$FLAG[as.numeric(PKdata$ID) == i][1] - 1,
                        FEMALE = PKdata$SEX[as.numeric(PKdata$ID) == i][1])
}

#load auxilliary modeling functions
#note: all parallel computing functions rely on the doMC package which is only available for Unix-based systems
#accordingly Windows-based systems require alternative pipelines for parallel computing (e.g. the snow package)
source("PSM_aux.R")
registerDoSEQ()

#determine the number of transit compartments (ODE models)
modt <- list()
for(i in 0:6){
  cat(paste("* Fitting deterministic model with fixed effects and",i,"transit compartments *\n"))
  pars <- c(stress = 45, ka = 0.10, kt= 0.3, ks=.60, ke = .10, init=0.001, scal=.8, sex=.7, S = 60, sigma=0.001)
  parA <- list(LB=pars*.2, Init=pars, UB=pars*2.5) #bounds + inits
  PK <- ODE_lm(nTransit = i) #determine number of transit compartment in structural model
  fitA <- PSM.estim(PK, PKdat, parA, trace=1, control=list(optimizer="ucminf"))
  modt[[i+1]] <- fitA
}
saveRDS(modt, "PSM_transits_ODE.RDS")

#assemble and inspect results of structural modeling
res <- Vectorize(function(x) round(modt[[x]]$THETA,3))(1:7) #parameter estimates for each fitted model
res <- rbind(res, "Mtt"= Vectorize(function(x) round((x)/modt[[x]]$THETA["kt"],3))(1:7)) #add mean transit time for each model
res <- rbind(res, "LL"= -Vectorize(function(x) modt[[x]]$NegLogL)(1:7))
res <- rbind(res, "AIC"= 2*c(6,rep(7,6)) +2*Vectorize(function(x) modt[[x]]$NegLogL)(1:7))
res <- rbind(res, "wAIC"= round(exp(-(res["AIC",] - min(res["AIC",]))/2) / sum(exp(-(res["AIC",] - min(res["AIC",]))/2)),4) )
print(res)

#determine the number of random effects (ODE population models)
modt <- readRDS("PSM_transits_ODE.RDS")
pars <- c(round(modt[[4]]$THETA,1), "OMEGA_stress"=1.5, "OMEGA_ke"=.05, "OMEGA_kt"=.15, "OMEGA_init"=.15)
pars[c("init","sigma")] <- c(.001,.001)
mods <- list()

for(i in 0:4){
  cat(paste("* Fitting deterministic model with",i,"random effects and 3 transit compartments *\n"))
  if(i == 1) pars["S"] <- 7.5
  parA <- list(LB=pars*.2, Init=pars, UB=pars*2.5) #bounds + inits
  PK <- ODE_lme(nTransit = 3, nOMEGAs = i)
  fitA <- PSM.estim(PK, PKdat, parA, trace=1, control=list(optimizer="ucminf"))
  mods[[i+1]] <- fitA
}

#add stochastic perturbations (SDE population model)
pars <- round(mods[[4]]$THETA,1)
pars[c("init","sigma")] <- c(.001,1)

parA <- list(LB=pars*.2, Init=pars, UB=pars*2.5) #bounds + inits
PK <- SDE_lme(nTransit = 3, nOMEGAs = 4)
fitA <- PSM.estim(PK, PKdat, parA, trace=1, control=list(optimizer="ucminf"))
mods[[6]] <- fitA

saveRDS(mods, "PSM_population_ODE.RDS")

#assemble results and generate parameter/fit table for population ODE/SDE models
i <- length(mods)
npars <- Vectorize(function(x) sum(parA$Init != round(mods[[x]]$THETA,3)))(1:i)

res <- Vectorize(function(x) round(mods[[x]]$THETA,3))(1:i) #parameter estimates for each fitted model
res <- rbind(res, "Mtt"= Vectorize(function(x) round(4/mods[[x]]$THETA["kt"],3))(1:i)) #add mean transit time for each model
res <- rbind(res, "LL"= -Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "AIC"= 2*npars +2*Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "wAIC"= round(exp(-(res["AIC",] - min(res["AIC",]))/2) / sum(exp(-(res["AIC",] - min(res["AIC",]))/2)),4) )
res <- rbind(res, "RMSEA"= sqrt((-2*res["LL",] +2*res["LL",i])/ npars - 1) / sqrt(dim(PKdata)[1]-1)  )
res <- rbind(res, "R2"= Vectorize(function(x) 1 - mods[[x]]$THETA["S"]/var(PKdata$DV))(1:i))
colnames(res) <- letters[1:i]
round(res,2)

#determine wisorized confidence intervals from the case-bootstrapped models a, e, f (250 replicates)
#bootstrapping requires massively parallel computation to obtain results in a reasonable amount of time
PK <- ODE_lme(nTransit = 3, nOMEGAs = 0)
CIs <- PSM.ci(mods[[1]], PK, PKdat, bootstrap=5, cores = 50)
saveRDS(CIs, "PSM_CIs.RDS")

PK <- ODE_lme(nTransit = 3, nOMEGAs = 4)
CIs <- PSM.ci(mods[[5]], PK, PKdat, bootstrap=5, cores = 50)
saveRDS(CIs, "PSM_CIs_ODE.RDS")

PK <- SDE_lme(nTransit = 3, nOMEGAs = 4)
CIs <- PSM.ci(mods[[6]], PK, PKdat, bootstrap=5, cores = 50)
saveRDS(CIs, "PSM_CIs_SDE.RDS")

CIs <- list()
CIs[["a"]] <- readRDS("PSM_CIs.RDS")
CIs[["e"]] <- readRDS("PSM_CIs_ODE.RDS")
CIs[["f"]] <- readRDS("PSM_CIs_SDE.RDS")

for(j in c("a","e","f")){
  for(i in names(res[1:14,j])){
    reps <- CIs[[j]]$replicates[,i]
    reps[reps < quantile(reps, .025)] <- quantile(reps, .025)
    reps[reps > quantile(reps, .975)] <- quantile(reps, .975)
    sdtemp <- sqrt(sum((reps - res[i,j])^2) / length(reps))
    CIs[[j]]$Bias[i] <- CIs[[j]]$SD[i] - sdtemp
    CIs[[j]]$SD[i] <- sdtemp 
    CIs[[j]]$CI[,i] <- qnorm(c(.025,.5,.975), mean = res[i,j], sd = sdtemp)
  }
  CIs[[j]]$dCI <- matrix(NA, ncol=3, nrow=3, dimnames = list(1:3, c("stress_f","init","Mtt")))
  
  reps <- CIs[[j]]$replicates[,"stress"] * exp(-CIs[[j]]$replicates[,"sex"])
  reps[reps < quantile(reps, .025)] <- quantile(reps, .025)
  reps[reps > quantile(reps, .975)] <- quantile(reps, .975)
  sdtemp <- sqrt(sum((reps - (res["stress",j]*exp(-res["sex",j])))^2) / length(reps))
  CIs[[j]]$dCI[,"stress_f"] <- qnorm(c(.025,.5,.975), mean = res["stress",j]*exp(-res["sex",j]), sd = sdtemp) 
  
  reps <- CIs[[j]]$replicates[,"ks"] / CIs[[j]]$replicates[,"ke"]
  reps[reps < quantile(reps, .025)] <- quantile(reps, .025)
  reps[reps > quantile(reps, .975)] <- quantile(reps, .975)
  sdtemp <- sqrt(sum((reps - (res["ks",j]/res["ke",j]))^2) / length(reps))
  CIs[[j]]$dCI[,"init"] <- qnorm(c(.025,.5,.975), mean = res["ks",j]/res["ke",j], sd = sdtemp)
  
  CIs[[j]]$dCI[,"Mtt"] <- rev(4 / CIs[[j]]$CI[,"kt"])
}

round(CIs[["a"]]$CI,2)
round(CIs[["a"]]$dCI,2)
round(CIs[["e"]]$CI,4)
round(CIs[["e"]]$dCI,2)
round(CIs[["f"]]$CI,2)
round(CIs[["f"]]$dCI,2)


#visualize and compare fits of structural and population models (Figure 5)
res <- Vectorize(function(x) round(modt[[x]]$THETA,3))(1:7) #parameter estimates for each fitted model
res <- rbind(res, "Mtt"= Vectorize(function(x) round((x)/modt[[x]]$THETA["kt"],3))(1:7)) #add mean transit time for each model
res <- rbind(res, "LL"= -Vectorize(function(x) modt[[x]]$NegLogL)(1:7))
res <- rbind(res, "AIC"= 2*c(6,rep(7,6)) +2*Vectorize(function(x) modt[[x]]$NegLogL)(1:7))
res <- rbind(res, "wAIC"= round(exp(-(res["AIC",] - min(res["AIC",]))/2) / sum(exp(-(res["AIC",] - min(res["AIC",]))/2)),4) )

par(mfrow=c(2,2), mar=c(3.5,3.5,1,1.5), mgp=c(2.2,0.8,0)) #10.00 x 10.00

#panel A
plot(0:6, res["wAIC",], type="b", col="black", pch=19, lty=2, lwd=2, ylab="Akaike weights", xlab="number of transit compartments", ylim=c(0,.35))
grid()
legend("bottomright", bty="n", lty=2, lwd=2, cex=.9, pch=19, legend=c("Akaike weight","log(likelihood)"), col=c("black","darkgrey") )
text(0.5, 0.325, "A", cex=2)
LL <- scale(res["LL",])/10 + .275
lines(0:6, LL, type="b", pch=19, lty=2, lwd=2, col="darkgrey")
mgp.axis(4, col="darkgrey", at=seq(0, .35, by = .05), cex.axis=.7, labels=round((seq(0, .35, by = .05)-.275)*10*attr(LL, "scaled:scale") + attr(LL, "scaled:center")), col.axis="darkgrey", col.ticks = "darkgrey", line=0)

#panel B
plot(c(-20,70),c(0,35), type="n", ylab="salivary cortisol (nM)", xlab="time relative to TSST onset (min)")
grid()
legend("topright", bty="n", pch=19, lty=2, lwd=2, cex=.9, col=c("tomato","darkblue","royalblue"), legend=c("Dresden ♀", "Dresden ♂‚", "Montreal ♂"))
text(-10,32.5, "B", cex=2)

obs <- with(PKdata[PKdata$FLAG == 1,], tapply(DV, TAD, mean))
obs_se <- with(PKdata[PKdata$FLAG == 1,], tapply(DV, TAD, sd))/sqrt(10)
errbar(x = unique(PKdata$TAD[PKdata$FLAG == 1]), y = obs, yplus = obs+1.96*obs_se, yminus = obs-1.96*obs_se, col="royalblue", errbar.col = "royalblue", add=T)
preds <- PSM.smooth(ODE_lm(nTransit = 3), PKdat[100], THETA = as.list(modt[[4]])$THETA, subsample = 20)
lines(preds[[1]]$Time, preds[[1]]$Ys, lwd=2, lty=2, col="royalblue")

obs <- with(PKdata[PKdata$FLAG == 2 & PKdata$SEX == 1,], tapply(DV, TAD, mean))
obs_se <-with(PKdata[PKdata$FLAG == 2 & PKdata$SEX == 1,], tapply(DV, TAD, sd))/sqrt(100)
errbar(x = unique(PKdata$TAD[PKdata$FLAG == 2])-.5, y = obs, yplus = obs+1.96*obs_se, yminus = obs-1.96*obs_se, col="tomato", errbar.col = "tomato", add=T)
dats <- PKdat[1] #select female participant, location Dresden
dats[[1]]$Time <- c(-20,dats[[1]]$Time)
dats[[1]]$Y <- cbind(NA,dats[[1]]$Y)
dats[[1]]$U <- cbind(c(0,1),dats[[1]]$U)
preds <- PSM.smooth(ODE_lm(nTransit = 3), dats, THETA = as.list(modt[[4]])$THETA, subsample = 20)
lines(preds[[1]]$Time, preds[[1]]$Ys, lwd=2, lty=2, col="tomato")

obs <- with(PKdata[PKdata$FLAG == 2 & PKdata$SEX == 0,], tapply(DV, TAD, mean))
obs_se <-with(PKdata[PKdata$FLAG == 2 & PKdata$SEX == 0,], tapply(DV, TAD, sd))/sqrt(100)
errbar(x = unique(PKdata$TAD[PKdata$FLAG == 2])+.5, y = obs, yplus = obs+1.96*obs_se, yminus = obs-1.96*obs_se, col="darkblue", errbar.col = "darkblue", add=T)
dats <- PKdat[6] #select male participant, location Dresden
dats[[1]]$Time <- c(-20,dats[[1]]$Time)
dats[[1]]$Y <- cbind(NA,dats[[1]]$Y)
dats[[1]]$U <- cbind(c(0,1),dats[[1]]$U)
preds <- PSM.smooth(ODE_lm(nTransit = 3), dats, THETA = as.list(modt[[4]])$THETA, subsample = 20)
lines(preds[[1]]$Time, preds[[1]]$Ys, lwd=2, lty=2, col="darkblue")

#panel C
fits <- PSM.smooth(ODE_lme(nTransit = 3, nOMEGAs = 4), PKdat, as.list(mods[[5]])$THETA, subsample = 20)
fitse <- PSM.smooth(SDE_lme(nTransit = 3, nOMEGAs = 4), PKdat, as.list(mods[[6]])$THETA, subsample = 20)
fitsb <- PSM.smooth(ODE_lme(nTransit = 3, nOMEGAs = 0), PKdat, as.list(mods[[1]])$THETA, subsample = 20)

pred_x <- predi_x <- predie_x <- obs_y <- NULL
for(i in 1:210){
  sel <- fitse[[i]]$Time %in% PKdat[[i]]$Time
  obs_y <- cbind(obs_y, PKdat[[i]]$Y)
  predi_x <- c(predi_x, fits[[i]]$Ys[sel])
  predie_x <- c(predie_x, fitse[[i]]$Ys[sel])
  pred_x <- c(pred_x, fitsb[[i]]$Ys[sel])
}
plot(obs_y, predi_x, col=rgb(0,0,0,.4), pch=19, ylab="predicted salivary cortisol (nM)", xlab="observed salivary cortisol (nM)", xlim=c(0,60), ylim=c(0,60))
grid()
points(predie_x, obs_y, col=rgb(135/255,206/255,235/255,.4), pch=19)
abline(a = 0, b = 1, lty=2, lwd=2, col="red")

legend("bottomright", bty="n", pch=c(19,19), cex=.9, col=c("black","skyblue"), legend=c("ODE","SDE"))
text(5,55, "C", cex=2)

#panel D
retain <- unlist(Vectorize(function(x) if(PKdat[[x]]$covar["FEMALE"] == 0 & PKdat[[x]]$covar["STUDY"] == 0){return(x)})(1:210))

plot(fits[[retain[9]]]$Time, fits[[retain[9]]]$Xs[4,]*mods[[5]]$THETA["ks"], lty=3, lwd=2, col="darkblue",  type="l", xlim=c(-20,70), ylim=c(0, 2.5), xlab="time relative to TSST onset (min)", ylab="cortisol secretion rate (nM/min)")
lines(fitse[[retain[9]]]$Time, fitse[[retain[9]]]$Xs[4,]*mods[[6]]$THETA["ks"], lty=1, lwd=2, col="darkblue")
lines(fits[[retain[8]]]$Time, fits[[retain[8]]]$Xs[4,]*mods[[5]]$THETA["ks"], lty=3, lwd=2, col="royalblue")
lines(fitse[[retain[8]]]$Time, fitse[[retain[8]]]$Xs[4,]*mods[[6]]$THETA["ks"], lty=1, lwd=2, col="royalblue")

grid()
legend("topright", bty="n", lty=rep(c(1,3), 2), lwd=2, cex=.9, col=rep(c("darkblue","royalblue"), each=2), legend=c(paste0("ID ", retain[9], ", SDE"), paste0("ID ", retain[9], ", ODE"), paste0("ID ", retain[8],", SDE"), paste0("ID ", retain[8], ", ODE")))
text(-10,2.25, "D", cex=2)

quartz.save("Figure 5.pdf", type = "pdf", dpi=300)

rm(dats, obs, obs_se, preds, LL)
rm(fits, fitse, fitsb, pred_x, predi_x, predie_x, obs_y)


#visual predictive checks (Figure 6; model e vs model f)
N <- 2000
simdat <- list()
for(i in 1:(N/2))  simdat[[i]] <- list("Time" = -10:75, "U" = matrix(c(ifelse(-10:75 == 1, 1, 0), ifelse(-10:75 == 100, 0, 1)), nrow=2, ncol=length(-10:75), byrow = T),
                                       "covar" = c("STUDY" = 1, "FEMALE" = 0))
for(i in (N/2+1):N)  simdat[[i]] <- list("Time" = -10:75, "U" = matrix(c(ifelse(-10:75 == 1, 1, 0), ifelse(-10:75 == 100, 0, 1)), nrow=2, ncol=length(-10:75), byrow = T),
                                         "covar" = c("STUDY" = 1, "FEMALE" = 1))
PKsimODE <- PSM.simulate(Model = ODE_lme(3,4), THETA = res[1:14,"e"], Data = simdat, deltaTime = .1)
PKsimSDE <- PSM.simulate(Model = SDE_lme(3,4), THETA = res[1:14,"f"], Data = simdat, deltaTime = .1)

predquant <- Vectorize(function(time, female, perc, simdata){
  if(female == 0) return(quantile(unlist(Vectorize(function(x) simdata[[x]]$Y[time+11])(1:(N/2))), perc))
  if(female == 1) return(quantile(unlist(Vectorize(function(x) simdata[[x]]$Y[time+11])((N/2+1):N)), perc))
}, vectorize.args = c("time"))

par(mfrow=c(1,2), mar=c(3.5,3.5,1,1.5), mgp=c(2.2,0.8,0))

plot(c(-10,75), c(0,70), type="n", xlab="time relative to TSST onset (min)", ylab="salivary cortisol (nM)", axes=F) #females
grid()
axis(1, at = seq(-10,70,10))
axis(2)
poly_y <- c(mA(predquant(time = -10:75, female = 1, perc = .025, PKsimODE)), mA(predquant(time = 75:-10, female = 1, perc = .975, PKsimODE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(.5,.5,.5,.5), border = NA) #ODE
lines(-10:75, mA(predquant(time = -10:75, female = 1, perc = .5, PKsimODE)), lwd=2, lty=2, col="darkgrey")
poly_y <- c(mA(predquant(time = -10:75, female = 1, perc = .025, PKsimSDE)), mA(predquant(time = 75:-10, female = 1, perc = .975, PKsimSDE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(135/255,206/255,235/255,.5), border = NA) #SDE
#with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], points(TAD, DV, col=rgb(0,0,0,.1), pch=19))
with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], vioplotm(DV[TAD == -5], DV[TAD == 11], DV[TAD == 20], DV[TAD == 30], DV[TAD == 40], DV[TAD == 55], DV[TAD == 70], at=c(-5,11,20,30,40,55,70), range = 0, wex=4, add=T, axes=F))
text(0, 60, "A", cex=2.5)
legend("topright", bty = "n", col=c("darkgrey","skyblue"), lwd=10, lty=1, legend=c("ODE","SDE"))

plot(c(-10,75), c(0,70), type="n", xlab="time relative to TSST onset (min)", ylab="salivary cortisol (nM)", axes=F) #females
grid()
axis(1, at = seq(-10,70,10))
axis(2)
poly_y <-c(mA(predquant(time = -10:75, female = 0, perc = .05, PKsimODE)), mA(predquant(time = 75:-10, female = 0, perc = .95, PKsimODE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(.5,.5,.5,.5), border = NA) #ODE
lines(-10:75, mA(predquant(time = -10:75, female = 0, perc = .5, PKsimODE)), lwd=2, lty=2, col="darkgrey")
poly_y <- c(mA(predquant(time = -10:75, female = 0, perc = .05, PKsimSDE)), mA(predquant(time = 75:-10, female = 0, perc = .95, PKsimSDE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(135/255,206/255,235/255,.5), border = NA) #SDE
#with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], points(TAD, DV, col=rgb(0,0,0,.1), pch=19))
with(PKdata[PKdata$SEX == 0 & PKdata$FLAG == 2,], vioplotm(DV[TAD == -5], DV[TAD == 11], DV[TAD == 20], DV[TAD == 30], DV[TAD == 40], DV[TAD == 55], DV[TAD == 70], at=c(-5,11,20,30,40,55,70), range = 0, wex=4, add=T, axes=F))
text(0, 60, "B", cex=2.5)

quartz.save("Figure 6.pdf", type = "pdf", dpi=300)


#compile parameter estimates and fit of population ODE/SDE models
modt <- readRDS("PSM_transits_ODE.RDS")
mods <- readRDS("PSM_population_ODE.RDS")
i <- length(mods)
pars <- c(round(modt[[4]]$THETA,1), "OMEGA_stress"=1.5, "OMEGA_ke"=.05, "OMEGA_kt"=.15, "OMEGA_init"=.15)
pars[c("init","sigma")] <- c(.001,.001)
parA <- list(LB=pars*.2, Init=pars, UB=pars*2.5) #bounds + inits
npars <- Vectorize(function(x) sum(parA$Init != round(mods[[x]]$THETA,3)))(1:i)

res <- Vectorize(function(x) round(mods[[x]]$THETA,3))(1:i) #parameter estimates for each fitted model
res <- rbind(res, "Mtt"= Vectorize(function(x) round(4/mods[[x]]$THETA["kt"],3))(1:i)) #add mean transit time for each model
res <- rbind(res, "LL"= -Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "AIC"= 2*npars +2*Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "R2"= Vectorize(function(x) 1 - mods[[x]]$THETA["S"]/var(PKdata$DV))(1:i))
colnames(res) <- letters[1:i]
print(res)


#simulate artificial cortisol data using the population SDE model (model f)
N <- 10000
samps <- -20:80
simdat <- list()
for(i in 1:(N/2))  simdat[[i]] <- list("Time" = samps, "U" = matrix(c(ifelse(samps == 1, 1, 0), ifelse(samps == 100, 0, 1)), nrow=2, ncol=length(samps), byrow = T),
                                       "covar" = c("STUDY" = 0, "FEMALE" = 0))
for(i in (N/2+1):N)  simdat[[i]] <- list("Time" = samps, "U" = matrix(c(ifelse(samps == 1, 1, 0), ifelse(samps == 100, 0, 1)), nrow=2, ncol=length(samps), byrow = T),
                                         "covar" = c("STUDY" = 0, "FEMALE" = 1))

set.seed(1234)
PKsim <- PSM.simulate(Model = SDE_lme(3,4), THETA = res[1:14,"f"], Data = simdat, deltaTime = 1)
pars <- matrix(Vectorize(function(x) c(res["stress","f"], res["ke","f"], res["kt","f"], res["ks","f"]/res["ke","f"]) * exp(PKsim[[x]]$eta[1:4]) * exp(c(-res["sex","f"]*simdat[[x]]$covar["FEMALE"],0,0,0)))(1:N), ncol=4, byrow=T)


#calculate non-compartmental parameters from artificial cortisol data and run principal component analysis
noncomps <- matrix(NA, ncol=12, nrow=N)
colnames(noncomps) <- c("Cmin","Cmax","Tmax","MaxMin","React%","Recov%","React","Recov","AUCg","AUCi","Cinit","Tmin")
freq <- 15

for(i in 1:N){
  noncomps[i,1] <- min(PKsim[[i]]$Y[PKsim[[i]]$Time %in% seq(0, 60, by=freq)]) #minimum concentration
  noncomps[i,2] <- max(PKsim[[i]]$Y[PKsim[[i]]$Time %in% seq(0, 60, by=freq)]) #maximum (peak) concentration
  noncomps[i,3] <- PKsim[[i]]$Time[which(PKsim[[i]]$Y == noncomps[i,2])] #time of peak concentration
  noncomps[i,4] <- noncomps[i,2] - noncomps[i,1] #maximum-minimum (max increase)
  noncomps[i,5] <- noncomps[i,2] / PKsim[[i]]$Y[PKsim[[i]]$Time == 0] #percent reactivity
  noncomps[i,6] <- PKsim[[i]]$Y[PKsim[[i]]$Time == 60] / noncomps[i,2] #percent recovery
  noncomps[i,7] <- noncomps[i,2] - PKsim[[i]]$Y[PKsim[[i]]$Time == 0] #abs reactivity
  noncomps[i,8] <- PKsim[[i]]$Y[PKsim[[i]]$Time == 60] - noncomps[i,2] #abs recovery
  secretion <- approxfun(seq(0, 60, by=freq), PKsim[[i]]$Y[PKsim[[i]]$Time %in% seq(0, 60, by=freq)])
  noncomps[i,9] <- as.numeric(integrate(secretion, lower=0, upper=60, subdivisions = 500L)[1])
  noncomps[i,10] <- noncomps[i,9] - PKsim[[i]]$Y[PKsim[[i]]$Time == 0]*60
  noncomps[i,11] <- PKsim[[i]]$Y[PKsim[[i]]$Time == -20] #initial concentration
  noncomps[i,12] <- PKsim[[i]]$Time[which(PKsim[[i]]$Y == noncomps[i,1])] #time of min concentration
}

tcomps <- round(cor(noncomps, pars, method = "spearman"),2)
round(apply(noncomps, 2, quantile, probs = c(.25,.5,.75)),2)

library(psych)
pc <- principal(cor(cbind(noncomps[,1:11],pars), method="spearman"), n.obs = N, nfactors = 4, rotate="varimax", scores=T)
print(pc)

#visualize results of simulation (Figure 7)
#Panel A
concs <- matrix(NA, ncol=length(PKsim[[1]]$Time), nrow=N)
for(i in 1:N) concs[i,] <- round(PKsim[[i]]$Y, 2)

tcors <- round(cor(concs, pars, method = "spearman"),2)
tcors <- apply(tcors, 2, mA)

par(mar = c(5, 4, 1, 2) + 0.1)
layout(matrix(c(1,1,2,3), ncol=2, nrow=2, byrow=T))

matplot(samps, tcors, type="l", ylim=c(-1,1), lty=c(1,2,2,1), lwd=2, ylab=paste("rank correlation with cortisol"), xlab="time relative to TSST onset (t - 20 min)", col = c("black","darkgrey"))
grid()
legend("topright", lty=c(1,1,2,2), lwd=2, col = c("black","darkgrey"), legend=c(expression(C[A](0)),expression(C[B](0)),expression(k[T]),expression(k[E])), bty="n")

text(-17.5, -.75, "A", cex = 2.5)

#Panel B
plot(pc$loadings[,1], pc$loadings[,2], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.75, col=c(rep(rgb(0,0,0,.25),11),rep(rgb(0,0,1,.25),4)), xlab="component 1 (secretory magnitude)", ylab="component 2 (steady state)")
xs <- c(.08,-.08,.08,0,-.1,-.1,.04,-.08,-.09,0, -.08,-.1,-.07,.07,-.1)
ys <- c(0,0,.05,-.08,0,0,-0,0,0,.1, -.05,0,0,-.05,.05)
for(i in 1:15) text(x = pc$loadings[i,1]+xs[i], y = pc$loadings[i,2]+ys[i], adj = c(.5,.5), c(rownames(tcomps),expression(C[A](0)),expression(k[E]),expression(k[T]),expression(C[B](0)))[i], cex=.75, col=c(rep("black",11),rep("royalblue",4))[i] )
abline(h=0, lty=2)
abline(v=0, lty=2)
grid()
text(-.8, .8, "B", cex = 2.5)

#Panel C
plot(pc$loadings[,4], pc$loadings[,3], xlim=c(-1,1), ylim=c(-1,1), pch=19, cex=.75, col=c(rep(rgb(0,0,0,.25),11),rep(rgb(0,0,1,.25),4)), xlab="component 3 (secretory delay)", ylab="component 4 (initial concentration)")
xs <- c(.08,.05,.05,0,0,0,.1,0,0,.1, .08,0,-.07,.07,-.1)
ys <- c(.05,-.05,.1,-.1,-.1,.07,0,-.07,0,0, 0,.1,.05,-.05,0)
for(i in 1:15) text(x = pc$loadings[i,4]+xs[i], y = pc$loadings[i,3]+ys[i], adj = c(.5,.5), c(rownames(tcomps),expression(C[A](0)),expression(k[E]),expression(k[T]),expression(C[B](0)))[i], cex=.75, col=c(rep("black",11),rep("royalblue",4))[i] )
abline(h=0, lty=2)
abline(v=0, lty=2)
grid()
text(-.8, .8, "C", cex = 2.5)

quartz.save("Figure 7.pdf", type = "pdf", dpi=300)


#power analyses (Figure 8)
IDs <- as.numeric(which(Vectorize(function(x) PKdat[[x]]$covar["STUDY"])(1:210) == 1)) #select Dresden sample

PKdatas <- data.frame(cbind(IDs,t(Vectorize(function(x) PKdat[[x]]$Y[-3:-2])(IDs))))
names(PKdatas) <- c("ID",PKdat[[1]]$Time[!is.na(PKdat[[1]]$Y)])

PKdatas$Cmin <- apply(PKdatas[,2:8], 1, min) #minimum concentration
PKdatas$Cmax <- apply(PKdatas[,2:8], 1, max) #maximum (peak) concentration
PKdatas$Tmax <- apply(PKdatas[,2:10], 1, function(x) as.numeric(names(x)[which(x[1:7] == x[9])])) #time of peak concentration
PKdatas$MaxMin <- apply(PKdatas[,2:10], 1, function(x) x[9] - x[8]) #maximum-minimum (max increase)
PKdatas$"React%" <- apply(PKdatas[,2:10], 1, function(x) x[9] / x[1]) #percent reactivity
PKdatas$"Recov%" <- apply(PKdatas[,2:10], 1, function(x) x[7] / x[9]) #percent recovery
PKdatas$React <- apply(PKdatas[,2:10], 1, function(x) x[9] - x[1]) #abs reactivity
PKdatas$Recov <- apply(PKdatas[,2:10], 1, function(x) x[7] - x[9]) #abs recovery
secretion <- function(x) approxfun(as.numeric(names(PKdatas)[2:8]), x)
PKdatas$AUCg <- apply(PKdatas[,2:8], 1, function(x) try(as.numeric(integrate(secretion(x), lower=-5, upper=70, subdivisions = 500L)[1])))
PKdatas$AUCi <- PKdatas$AUCg - PKdatas[,2]*75
PKdatas$Cinit <- PKdatas[,2] #initial concentration

PKdatas$SEX <- as.numeric(Vectorize(function(x) PKdat[[x]]$covar["FEMALE"])(IDs))
etas <- PSM.smooth(SDE_lme(3,4), PKdat, mods[[6]]$THETA, trace = 1)
etas <- t(Vectorize(function(x) etas[[x]]$eta)(IDs))
PKdatas <- cbind(PKdatas,etas)
names(PKdatas)[21:24] <- c("R0","kE","kT","C0")
PKdatas$R0 <-  PKdatas$R0 - mods[[6]]$THETA["sex"] * PKdatas$SEX

round(cor(PKdatas[,c(-8:-1,-24:-21)], PKdatas[,c(20:24)], method="spearman"),2)

N <- c(5,seq(10,100,by=10))
reps <- 100000
store_AUCg <- store_AUCi <- store_Tmax <- store_MaxMin <- store_Cmin <- store_Cinit <- store_R0 <- store_Cmax <-  matrix(NA, nrow=reps, ncol=length(N))
store_C0 <- store_kE <- store_kT <-  matrix(NA, nrow=reps, ncol=length(N))

pb <- txtProgressBar(min = 0, max = length(N), style = 3)
for(j in 1:length(N)){
  for(i in 1:reps){
    temp <- PKdatas[PKdatas$ID %in% sample(PKdatas$ID[PKdatas$SEX == 0], size = N[j], replace = T) | PKdatas$ID %in% sample(PKdatas$ID[PKdatas$SEX == 1], size = N[j], replace = T),]
    store_R0[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$R0, method="spearman", exact = T)$p.value)
    store_C0[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$C0, method="spearman", exact = T)$p.value)
    store_kE[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$kE, method="spearman", exact = T)$p.value)
    store_kT[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$kT, method="spearman", exact = T)$p.value)
    store_AUCg[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$AUCg, method="spearman", exact = T)$p.value)
    store_AUCi[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$AUCi, method="spearman", exact = T)$p.value)
    store_Tmax[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$Tmax, method="pearson", exact = T)$p.value)
    store_MaxMin[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$MaxMin, method="spearman", exact = T)$p.value)
    store_Cmin[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$Cmin, method="spearman", exact = T)$p.value)
    store_Cmax[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$Cmax, method="spearman", exact = T)$p.value)
    store_Cinit[i,j] <- suppressWarnings(cor.test(temp$SEX, temp$Cinit, method="spearman", exact = T)$p.value)
  }
  setTxtProgressBar(pb, j)
}

plot(spline(N, apply(store_MaxMin, 2, function(x) sum(x < .05) / reps)), lwd=2, xlab="sample size per group", ylab="statistical power", type="l", ylim=c(0,1), xlim=c(5,100))
grid()
lines(spline(N, apply(store_AUCi, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=3)
lines(spline(N, apply(store_AUCg, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=4)
lines(spline(N, apply(store_Cmax, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=2)
lines(spline(N, apply(store_Tmax, 2, function(x) sum(x < .05, na.rm=T) / (reps-sum(is.na(x))))), lwd=2, lty=1, col="darkgrey")
lines(spline(N, apply(store_Cmin, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=2, col="darkgrey")
lines(spline(N, apply(store_Cinit, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=3, col="darkgrey")
lines(spline(N, apply(store_R0, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=1, col="royalblue")
legend(80,.75, bty="n", lty=c(1,1:4,1:3), lwd=2, col=c("royalblue",rep("black",4),rep("darkgrey",3)), cex=1,
       legend=c("R(0)","MaxMin","Cmax","AUCi","AUCg","Tmax","Cmin","Cinit"))

quartz.save("Figure 8.pdf", type = "pdf", dpi=300)

plot(spline(N, apply(store_R0, 2, function(x) sum(x < .05) / reps)), lwd=2, col="royalblue", xlab="sample size per group", ylab="statistical power", type="l", ylim=c(0,1), xlim=c(5,100))
grid()
lines(spline(N, apply(store_C0, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=2, col="royalblue")
lines(spline(N, apply(store_kE, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=3, col="royalblue")
lines(spline(N, apply(store_kT, 2, function(x) sum(x < .05) / reps)), lwd=2, lty=4, col="royalblue")
legend(80,.75, bty="n", lty=1:4, lwd=2, col="royalblue", legend=c("R(0)","C(0)","kE","kT"))
