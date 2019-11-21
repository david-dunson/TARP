# TARP: 
### Targeted Random Projection for Prediction from High-Dimensional Features 

https://www.tandfonline.com/doi/abs/10.1080/01621459.2019.1677240

*Minerva Mukhopadhyay and David B. Dunson*

**Abstract** We consider the problem of computationally-efficient prediction with high dimensional and highly correlated predictors when accurate variable selection is effectively impossible. Direct application of penalization or Bayesian methods implemented with Markov chain Monte Carlo can be computationally daunting and unstable. A common solution is first stage dimension reduction through screening or projecting the design matrix to a lower dimensional hyper-plane. Screening is highly sensitive to threshold choice, while projections often have poor performance in very high-dimensions. We propose TArgeted Random Projection (TARP) to combine positive aspects of both strategies. TARP uses screening to order the inclusion probabilities of the features in the projection matrix used for dimension reduction, leading to data-informed sparsity. We provide theoretical support for a Bayesian predictive algorithm based on TARP, including statistical and computational complexity guarantees. Examples for simulated and real data applications illustrate gains relative to a variety of competitors.

**Instructions**

`TARP.R` computes RIS-RP and RIS-PCR for normal-linear model

## Usage

input: `x` : standardized taining design matrix in nXp format

input: `y` : training response, a n-vector

input: `x_test` : standardized test design matrix in nXp format

## Examples
	RISRP_res=RIS_PCR(x,y,x_test,alpha=0.95)
	y_hat=RIS_PCR(x,y,x_test,alpha=0.95)
	y_hat=RISRP_res[,1]     # predicted value
	lower_PI=RISRP_res[,2]  # lower 100(1-alpha)% confidence limit
	upper_PI=RISRP_res[,3]  # upper 100(1-alpha)% confidence limit

	plot(y_test,type="l",lwd=3,ylim=c(-10,10),xlab="Predicted and observed values")
	par(new=TRUE)
	plot(y_hat,type="l",lwd=3,ylim=c(-10,10),col="blue",xlab="")
	par(new=TRUE)
	plot(lower_PI,type="l",lwd=3,ylim=c(-10,10),col="red",lty=2,xlab="")
	par(new=TRUE)
	plot(upper_PI,type="l",lwd=3,ylim=c(-10,10),col="red",lty=2,xlab="")




