# A Bayesian framework to model variance of grain yield response to plant density
 
This repository contains the codes and algorithms developed for the article "A Bayesian framework to model variance of grain yield response to plant density"  to be published in Plant Methods.

Code Authors:
Nicolas Giordano, Dustin Hayes

Manuscript authors:
Nicolas Giordano, Dustin Hayes, Trevor J. Hefley, Josefina Lacasa, Brian Beres, Lucas Haag, Romulo P. Lollato


**24. run MCMC.qmd:** Runs model using a metropolis-hastings sampling algorithm and divide and conquer approach across 12 splits of the data.

**28.visualization.qmd:** This code reproduced all figures of the manuscript.

**functions.R;** Contains functions to adjust plotting aesthetics. And function compare yields and variance among treatments at the expected value of AOPD and minRPD.

**MCMC REFINED.R** Contains hard-coded metropolis-hastings sampling algorithm. We include a set of functions to retrieve model derived quantities and posterior predictive distributions. 
