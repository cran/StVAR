\name{StVAR-package}
\alias{StVAR-package}

\alias{StVAR}

\docType{package}

\title{
Student's t Vector Autoregression (StVAR)
}

\description{
Maximum likelihood estimation of StVAR model is the purpose of this function. It can be used to estimate the linear autoregressive function (conditional mean) and the quadratic 
autosckedastic function (conditional variance).
}


\usage{
StVAR(Data, lag = 1, v = 1, maxiter = 1000, meth = "BFGS", hes = "FALSE", init = "na", tre = "C")
}



\arguments{
  \item{Data}{data matrix with at least two columns. Cannot be empty.}

  \item{lag}{a positive integer (default value is 1).}

  \item{v}{scalar (default value is 1) greater than or equal to 1. Degrees of freedom parameter.}

  \item{maxiter}{maximum number of iteration. Must be an integer bigger than 10.}

  \item{meth}{
one of the optimization method from \code{optim} function (default value is BFGS). See details of \code{optim} function.}

  \item{hes}{logical (default value is FALSE). If TRUE produces estimated hessian matrix and the standard errors of estimates.}

  \item{init}{if na (default), initial values for optimization are generated from a uniform distribution. A vector of initial values can also be used
(not recommended). The length of the init vector must be equal to the number of parameters of the joint distribution.}

  \item{tre}{C for constant (default), L for linear trend, Q for quadratic trend.}
}


\details{For the functional form of the autoregressive function and the autoskedastic function, see Spanos (1994). 
\code{StVAR} also produces some mis-specification test results. For detail, see Poudyal (2012). 
The trend (t=1,2,3,...) is divided by 1000 for estimation.}


\value{
\item{beta}{coefficients of the autoregressive function including the coefficients of trends in matrix form.}
\item{coef}{coefficients of the autoregressive function, standard errors and p-values. If some of the standard errors are 
NA's, the \code{StVAR()} function has to be run again.}	
\item{like }{maximum log likelihood value.}
\item{sigma}{contemporary variance covariance matrix.}
\item{cvar}{\code{(v/(v+lag*l-2))*sigma*cvar} is the fitted value of the autoskedastic function where \code{l} is the rank of \code{Data}}
\item{res}{nonstandardized residuals}
\item{fitted }{fitted values of the autoregressive function.}
\item{ms}{mis-specification test results.}	
\item{dist}{distributional test results. }
\item{init}{estimates of the joint distribution parameters. It can be used as new initial value \code{init} in \code{StVAR()} to improve optimization further.}
\item{S}{variance covariance matrix of the joint distribution.}
}


\author{
Niraj Poudyal \email{nirajp6@vt.edu}
}



\references{
Poudyal, N. (2012), Confronting Theory with Data: the Case of DSGE Modeling. Doctoral dissertation, Virginia Tech.

Spanos, A. (1994), On Modeling Heteroskedasticity: the Student's t and Elliptical
Linear Regression Models. Econometric Theory, 10: 286-315.
}


\examples{
##1. Stationary data
set.seed(1234)
X<-cbind(rt(70,3),rt(70,3))
init<-c(0.045,0.289,0.065,0.687,0.030,-0.361,0.186,0.326,-0.770,-0.229,-0.105,0.136,0.0731,0.281)
model<-StVAR(X,v=3,lag=1,tre="C",hes="FALSE",init=init)
plot(X[,1]) ;lines(model$fitted[,1]) # Actual value and fitted value
plot(model$cvar,type='l') #Conditional variance fit

##2. Non-stationary (linear trend) data
set.seed(3)
X<-cbind(0.1*seq(1,80,1)+rt(80,3),-0.2*seq(1,80,1)+rt(80,3))
init<-c(-0.472,0.093,0.254,0.134,0.181,-0.479,0.203,-0.685,-0.408,0.100,0.276,-0.599,-0.247,-0.358,105.514,-200.332)
model<-StVAR(X,v=3,lag=1,tre="L",hes="TRUE",init=init)

## The model captures the trend in data. p-values in mis-specification tests
##are expected to be bigger than 0.05. 
model$ms
model$dist

plot(X[,1]) ;lines(model$fitted[,1]) # Actual value and fitted value
plot(model$cvar,type='l') #Conditional variance fit

##3. Estimating a misspecified (stationary) model when 
## data has linear trends
set.seed(1234)
X<-cbind(0.1*seq(1,80,1)+rt(80,3),-0.2*seq(1,80,1)+rt(80,3))
init<-c(-0.724,0.600,-1.526,-0.551,1.162,-2.148,-0.377,1.294,-1.555,-0.321,0.340,-1.748,4.068,-7.689)
ms.model<-StVAR(X,v=3,lag=1,tre="C",hes="TRUE",init=init)

## Notice the low p-values for t-invariance.
ms.model$ms
ms.model$dist

plot(X[,1]) ;lines(ms.model$fitted[,1]) # Actual value and fitted value
plot(ms.model$cvar,type='l') #Conditional variance fit

}
