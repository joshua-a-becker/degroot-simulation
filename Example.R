### THIS CODE RUNS THE MODEL DESCRIBED IN
### Becker, J., Brackbill, D., & Centola, D. (2017). Network dynamics of social influence in the wisdom of crowds. 
### Proceedings of the national academy of sciences, 114(26), E5070-E5076.

rm(list=ls());gc()
require(ggplot2)
require(igraph)
require(mvtnorm)
require(dplyr)

### SET YOUR WORKING DIRECTORY
### THEN LOAD THE CORE OF THE MODEL
source('DeGroot_Function.R', echo=F)

generatePopulation = function(g, corr, truth=0){
  sigma <- matrix(c(1, abs(corr), abs(corr), 1), nrow = 2)
  x <- rmvnorm(vcount(g), mean = c(0,0), sigma = sigma, method = "chol")
  V(g)$guess=x[,1]
  V(g)$alpha = pchisq((x[,2]-truth)^2,1)
  
  # invert for negative correlation
  if( (corr/abs(corr))<0 ) { V(g)$alpha = 1-V(g)$alpha}
  
  g$corr = cor(abs(V(g)$guess-truth), V(g)$alpha)
  g
}

N=40
trial = 0
results = data.frame(  change_err_mu = numeric()
                     , change_ind_err = numeric()
                     , corr=numeric()
                     , trial=numeric()
                     , truth=numeric()
                     , centralization=numeric()
                     , net=character()
                     , stringsAsFactors=F
                     )



### 1. TEST THE EFFECT OF CORRELATION BTW ACCURACY/ALPHA 
###    IN A DECENTRALIZED POPULATION
replications=50
for(truth in c(0,1)) {
  for(corr in seq(-0.99, 0.99, by=0.2)) {
    print(trial)
    for(i in 1:replications) {
      trial=trial+1
      g=degree.sequence.game(rep(4,N))
      pop=generatePopulation(g, corr, truth=truth)
      output = simDegroot(pop, 5, truth)
      
      output = output %>% mutate(
        err_mu = abs(truth-mean)
      )
      
      change_err_mu = output$err_mu[output$round==5]-output$err_mu[output$round==0]
      change_ind_err = output$mean_ind_err[output$round==5]-output$mean_ind_err[output$round==0]
      results[nrow(results)+1,]=data.frame(
        change_err_mu = change_err_mu
        , change_ind_err = change_ind_err
        , corr = corr
        , trial=trial
        , truth = unique(output$truth)
        , centralization=unique(output$deg.cent)
        , net=g$name
        , stringsAsFactors = F
      )
    }
  }
}

### 2. TEST THE EFFECT OF CENTRALIZATION ON GROUP ACCURACY
###    WITH CORR=0, ERR=0

replications=50
for(truth in c(0,1)) {
  for(pow in seq(0,3,by=0.1)) {
    corr = 0.001
    truth = 0
      print(trial)
      for(i in 1:replications) {
        trial=trial+1
        g=barabasi.game(n=N, pow=pow, m=2, directed=F)
        pop=generatePopulation(g, corr, truth=truth)
        output = simDegroot(pop, 5, truth)
        
        output = output %>% mutate(
          err_mu = abs(truth-mean)
        )
        
        change_err_mu = output$err_mu[output$round==5]-output$err_mu[output$round==0]
        change_ind_err = output$mean_ind_err[output$round==5]-output$mean_ind_err[output$round==0]
        results[nrow(results)+1,]=data.frame(
          change_err_mu = change_err_mu
          , change_ind_err = change_ind_err
          , corr = pop$corr
          , trial=trial
          , truth = unique(output$truth)
          , centralization=unique(output$deg.cent)
          , net=g$name
          , stringsAsFactors = F
        )
      }
  }
}

### PLOT EFFECT OF ACCURACY CORELATION
results$truth_label = paste0("Expected Error = ", results$truth)
ggplot(subset(results,net=="Degree sequence random graph"),
       aes(x=round(corr,1), y=change_err_mu) ) +
  stat_summary(fun.y="mean", geom="point") +
  geom_hline(yintercept=0) + 
  facet_grid(.~truth_label) + 
  labs(x="Correlation", y="Change in Error of Mean")

ggplot(subset(results,net=="Degree sequence random graph"),
       aes(x=round(corr,1), y=change_ind_err) ) +
  stat_summary(fun.y="mean", geom="point") +
  geom_hline(yintercept=0) + 
  facet_grid(.~truth_label) + 
  labs(x="Correlation", y="Change in Individual Error")


### PLOT EFFECT OF ACCURACY CENTRALIZATION
ggplot(subset(results,net=="Barabasi graph"),
       aes(x=round(centralization,1), y=change_err_mu) ) +
  stat_summary(fun.y="mean", geom="point") +
  geom_hline(yintercept=0) + 
  facet_grid(.~truth_label) + 
  labs(x="Centralization", y="Change in Error of Mean")


ggplot(subset(results,net=="Barabasi graph"),
       aes(x=round(centralization,1), y=change_ind_err) ) +
  stat_summary(fun.y="mean", geom="point") +
  geom_hline(yintercept=0) + 
  facet_grid(.~truth_label) + 
  labs(x="Centralization", y="Change in Individual Error")

