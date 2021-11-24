# install.packages('coda')   
library(coda)
library(statnet)
data(package='ergm')
data("faux.desert.high")
fmh <- faux.desert.high

#  a few data summaries
fmh
plot(fmh, displayisolates = FALSE)
table(component.dist(fmh)$csize)

plot(fmh, displayisolates = FALSE, vertex.col = "grade")
fmh.degreedist <- table(sna::degree(fmh, cmode = "indegree"))
fmh.degreedist
summary(fmh ~ odegree(0:8))
summary(fmh ~ odegree(0:8, "sex"))
summary(fmh ~ triangle)
summary(fmh ~ edges + triangle)   
mixingmatrix(fmh, "grade")   # cross-classification of ego and alter grades
gr <- fmh %v% "grade"   
table(gr)  # distribution of ego grade

#  Time to model
# set a directory to store some output graphs
setwd('/Users/davidschaefer/Desktop')
setwd('/Users/drschaef/Desktop')
setwd('C:/temp/')


model1 <- ergm(fmh ~ edges)
summary(model1)     # Show output
names(model1)       # List all items that accompany output
#  No need to evaluate degneracy of model becuase it is a dyad inedependence model (not an ERGM)
#  (We will evaluate the fit of model1 below below)

# control for limits on degree (max outdegree = 11)
model1b <- ergm(fmh ~ edges, constraints=~bd(maxout=11))
summary(model1b)

# control for mutuality
model1c <- ergm(fmh ~ edges + mutual, constraints=~bd(maxout=11))
summary(model1c)

# control for homophily, sociality 
model1d <- ergm(fmh ~ edges + nodematch('grade') + nodeocov('grade') + nodeicov('grade'),
	constraints=~bd(maxout=11))
summary(model1d)

# control for homophily, sociality, mutuality 
model1e <- ergm(fmh ~ edges + mutual + nodematch('grade') + nodeocov('grade') + nodeicov('grade'),
	constraints=~bd(maxout=11))
summary(model1e)

# compare coefficients across models
(sum1 <- round(cbind(coef(model1e), coef(model1d)[names(coef(model1e))], coef(model1c)[names(coef(model1e))]),3))


model1f <- ergm(fmh ~ edges + mutual + nodematch('grade') + nodeocov('grade') + nodeicov('grade') +
	nodematch('race') + nodematch('sex'),
	constraints=~bd(maxout=11))
summary(model1f)

effs <- names(coef(model1f))
(sum2 <- round(cbind(coef(model1f), coef(model1e)[effs], coef(model1d)[effs], coef(model1c)[effs]), 3))



#  Estimate a model that tests transitivity, includes lower-order structures
model.tt <- ergm(fmh ~ edges + ostar(2) + istar(2) + m2star + triangle)
#  This model is so degenerate that it won't even estimate



#  Estimate models with better (i.e., gw) degree and triad closure terms 
#  - gwesp instead of triangles
model2 <- ergm(fmh ~ edges + gwesp(0, fixed = TRUE), 
  control = control.ergm(MCMLE.steplength=.25, MCMC.samplesize = 10000, MCMLE.maxit=5),
  verbose = TRUE)  # to get less output, set verbose to FALSE
summary(model2)
# gives output, but is it a good model

### EVALUATE DEGENERACY
#  To look for degeneracy, observe the model fitting process
pdf('temp2.pdf')
mcmc.diagnostics(model2)
dev.off()
#    Graphs: lefthand show how statistics changed in last iteration (relative to all)
#    * numbers of statistics should fluctuate around 0 (smooth line), not depart
#    * trend away from 0 indicates degeneracy
#    Text: Sample statistics should NOT differ from observed
#          Examine correlations at Lag 100+; these should be less than .4
#          Non-significant GEWEKE statistics: compares means from begin and end of Markov chain
#  Conclusion: degenerate model

#  Recommended solutions to degeneracy (assuming model is not mis-specified)
#  Don't try these now - they take a long time 
#    1. Increase the number of networks sampled, default is 10,000
#    2. Increase interval b/w sampled networks, default is 100 (give them more time to change)
#    3. Increase number of iterations from default of 20
#    4. Decrease the step-length from default of .5
#    5. Increase MCMC chain burnin from default of 10,000 


#  In addition, to eliminate degeneracy, use appropriate set of parameters
#  NOTE: The following models take a long time to run
#  I have decreased maxit to 5 in the next two models - it really should be 20
#  I have also decreased the number of samples to 10,000 (1e+4) - it should be 100,000 (1e+5)
#  Estimate a model using GWESP, alpha=0 (~alternating k-triangles)
model4.take1 <- ergm(fmh ~ edges + nodematch("grade") + nodematch("race") + nodematch("sex") +
     mutual + gwesp(0, fixed = TRUE), MCMCsamplesize = 1e+5, maxit = 5, verbose = TRUE)
summary(model4.take1)
#  Examine model output to look for degeneracy
#    decreasing log-likelihood improvements indicate convergence on a solution
#  Simulate a model and compare number of triangles (sim4.take1 is a network object)
sim4.take1 <- simulate(model4.take1, burnin = 1e+6)  
c(fmh = summary(fmh ~ triangle), sim4.take1 = summary(sim4.take1 ~ triangle))

#  Observe the model fitting process
#    we save as pdf because there are too many graphs to fit on screen
pdf("model4take1diagnostics.pdf") 
mcmc.diagnostics(model4.take1)
dev.off()
#  No sign of degeneracy - proceed to evaluating model fit

# Note: if the following doesn't work switch to model2
###  EVALUATE FIT OF MODEL - observed vs. simulated network
#  Begin by evaluating the fit of model1b (which was degenerate)
#  Older way - useful to understand what is being done
#    simulate one network using estimated parameters (could simulate as many networks as we want)
sim1b <- simulate(model1b, burnin = 1e+6)  
#  plot degree distribution of network, followed by simulated network 
plot(summary(sim1b ~ odegree(0:30)), type = "l", lty = 2, lwd = 3, xlab = "Degree", ylab = "Count")
lines(summary(fmh ~ odegree(0:30)), lty = 1, lwd = 2)
legend("topright", legend = c("Observed", "Simulated"), lwd = 3, lty = 1:2)
title("Outdegree Distribution, model1b")
#  too many actors with low degree in the simulated network
#    compare number of triangles in observed vs. simulated networks (too many simulated)
c(fmh = summary(fmh ~ triangle), sim1b = summary(sim1b ~ triangle))
#  not enough triangles - as we saw above, this model is degenerate

#  Newer way to evaluate fit - expect this step to take a while
#  Simulate several graphs and record the statistics listed
#    (note: the numbers below are smaller than one should use in practice)
gof.model1b <- gof(model1b)
#  second, print and then graph the observed and simulated networks
gof.model1b
pdf("gof_model1b.pdf") 
plot(gof.model1b)
dev.off()
#  open the "gof_model1b.pdf" file to view the fit
#    compare the observed distribution (dark line) to the simulated distribution (boxplots)
#      confidence intervals are indicated by light gray lines
#    model1b produced too few low degree nodes, too much & too little esp and too low geodesic distances

# Evaluate fit of model 4, which was not degenerate 
gof.model4.take1 <- gof(model4.take1)
#  second, print and then graph the observed and simulated networks
gof.model4.take1
pdf("gof_model4.take1.pdf") 
plot(gof.model4.take1)
dev.off()
#  open the "gof_model4.take1.pdf" file to view the fit
#  much better for all statistics, esp still a bit off


#  Improving the fit of model 4
#  Try model 4 with different alpha scaling paramters for the GWESP effect (0.5)
model4.take2 <- ergm(fmh ~ edges + nodematch("grade") + nodematch("sex") +
     mutual + gwesp(0.5, fixed = TRUE), MCMCsamplesize = 1e+4, maxit = 5, verbose = TRUE)
pdf("model4take2diagnostics.pdf")  
mcmc.diagnostics(model4.take2)
dev.off()
summary(model4.take2)
gof.model4.take2 <- gof(model4.take2)
gof.model4.take2
pdf("gof_model4.take2.pdf") 
plot(gof.model4.take2)
dev.off()

#  Try model 4 with terms for isolates and in/outdegree distributions
model4.take3 <- ergm(fmh ~ edges + nodematch("grade") + nodematch("sex") + 
     mutual + gwesp(0, fixed = TRUE) + isolates + gwodegree(.2, fixed = TRUE) +
     gwidegree(.2, fixed = TRUE), MCMCsamplesize = 1e+4, maxit = 15, verbose = TRUE)
pdf("model4take3diagnostics.pdf") 
mcmc.diagnostics(model4.take3)
dev.off()
summary(model4.take3)
gof.model4.take3 <- gof(model4.take3)
gof.model4.take3
pdf("gof_model4.take3.pdf") 
plot(gof.model4.take3)
dev.off()

#  Compare fit across the models
#    these are the estimated log-likelihood - values closer to zero indicate better fit
c(model4.take1$mle.lik, model4.take2$mle.lik, model4.take3$mle.lik)

# Learn more about ergm package
?ergm  # lots of information here; click on 'Index' link at bottom for all functions

# What terms are possible with ERGM
# - from the index page, click on 'ergm.terms'
# - or from ?ergm, navigate to 'ergm-terms' page
#  - you can also search for specific keywords in the ergm terms descriptions

search.ergmTerms('triangle')
search.ergmTerms('homophily')




