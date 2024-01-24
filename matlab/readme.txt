The optdesign package is Copyright (C) Juha Karvanen 2005

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Version: 1.1
Version date: May 6, 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Description:
===========

This package provides the Matlab (5.x or later) functions needed for 
experimental designs described in

	J. Karvanen, J. J. Vartiainen, A. Timofeev, J. Pekola, Experimental 
	designs for binary data in switching measurements on superconducting 
	Josephson junctions. Journal of the Royal Statistical Society: Series C 
	(Applied Statistics) 56 (2), 167–181, 2007.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Background:
==========

1) What is experimental design? optimal design? sequential design?

Experimental design is a plan how a test, a measurement or an experiment is
performed, especially how the values of controllable factors (covariates) are 
selected. Experimental design is the general term; optimal design and 
sequential design are experimental designs with some specific properties.

Optimal design tells which covariate values are the best according 
to some optimality criterion. For instance, often used D-optimality criterion
minimizes the generalized variance of the estimated parameters. 

Sequential design refers to situation where the experiment is performed in
several stages. The parameter estimates are updated after each stage and 
the new estimates are utilized in the design for the next stage. 

More information can be found, for example, from the references of the paper 
mentioned above. 

2) What are switching measurements? Josephson junctions?

The Josephson junctions (JJ) are important non-linear components of superconducting 
electronics. A junction involves two superconducting electrodes separated by an 
insulator gap. In practice, miniature JJ circuits are manufactured by lithographic 
means and operated at cryogenic temperatures, i.e., below 4.2 K. 

An experiment called switching measurement is a common way to probe the properties 
of a JJ circuit sample. In the experiment, sequences of current pulses are applied 
to the sample, while the voltage over the structure is monitored. The physical state 
of the JJ is described by the phase difference over the junction. A sufficiently 
high current pulse causes the  phase difference  to escape by tunneling from the 
local minimum of the potential to the direction of the steepest descent. This is 
called switching. In the experiment the observed voltage pulse indicates that 
switching has occurred.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Installation: Just put all files to a directory along Matlab's search path.
============

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Quick use: Type optdesign_demo for a demonstration.
=========


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Use: 
=== 

A) optdesign function
~~~~~~~~~~~~~~~~~~~~~~~~

The program calculates the optimal design for swithching measurements in 
superconducting Josephson junctions. More generally, the program calculates 
the optimal covariate points for binary response with complementary
log-log link. The optimality is defined as D-optimality and parameters are estimated 
by maximum likelihood estimation. The data is assumed to follow
Gompertz/Gumbel/extreme value distribution

               P(Y=1) = 1-exp(-exp(a*x+b)),

where a and b are the model parameters, x is the covariate variable and Y is the binary 
response. In the case of one Josephson junction this model arises from the laws of
physics. All data collected so far is given as an input to the program which then returns 
the optimal covariate points to be measured in next stage.

Sequential design with optdesign.m:

	1) Obtain the limits of the initial interval: xmin and xmax
	2) Measure at points xmin and xmax. Let x contain the points
	   of measurement and y contain the measured responses. 
	   (Sometimes it can be directly assumed that the response for 
	    xmin is 0 and the response for xmax is 1.)
	3) Call opdesign.m
	   [optx,betaest,betacov,IsMLE]=optdesign(x,y)
	4) If the required accuracy of betaest has been achieved, stop the 
	   procedure; otherwise proceed to the next step.
	5) Measure at the points given by optx. Append the points of measurement 
	   to x. Append the measured responses to y.
	6) Go to step 3.

The function is intended to be used for sequential designs of switching 
measurements but it may used for other suitable purposes as well. 


Uses:          gompertz_logl.m
               gompertz_cov.m
               gompertz_cdf (through gompertz_logl.m)

Input:
xx             Points (pulse heights) measured so far, 
                   size(xx)=[n 1], where n is the total number of points
                   measured so far.
yy             Measured binary response for points xx, 0 or 1, size(yy)=[n 1].
beta0          (optional) initial values of a and b for maximum likelihood
                   estimation, beta0=[a0 b0], lenght(beta0)=2.
transform		(Not in use in this version)	

Output:
optx			Points (pulse heights) to be measured at the next stage, length(optx)=2
                   If IsMLE==1 measure at both optx(1) and optx(2)
                   If IsMLE==0 measure only at optx(1).
betaest		Estimated values of parameters a and b
                   If IsMLE==0, an ad-hoc estimate is returned.
betacov        Estimated covariance matrix for parameters a and b,
                   size(betacov)=[2,2]
                   If IsMLE==0, betacov cannot be calculated and a matrix of NaN:s is returned.
IsMLE          Are the estimates maximum likelihood estimates? 0 or 1.





B) Use of the other functions and the list of files in the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contents.m
  A brief description of all the files in the package.


gompertz_cdf.m
  Calculates the cumulative distribution function of the Gompertz/Gumbel 
  distribution
  
  Syntax: p=gompertz_cdf(beta_,x)

  Input:
	beta_(1)   scale parameter of the Gompertz distribution = a
	beta_(2)   location parameter of the Gompertz distribution = b
	x          vector of points where the value of the cumulative
	           distribution function is calculated.
  Output:
	p          vector of the values of the cumulative distribution
	           function in points x.


gompertz_cov.m
  Calculates the covariance matrix for parameters a and b in complementary log-log regression
  
  Syntax: lcov=gompertz_cov(beta_,x,y)

  Input:
	beta_(1)   scale parameter of the Gompertz distribution = a
	beta_(2)   location parameter of the Gompertz distribution = b
	x          vector of the covariate points
	y          vector of the response values

	Output:
	lcov       The (asymptotic) covariance matrix for beta_ 



gompertz_logl.m
  Calculates the log-likehood for complementary log-log regression

  Syntax: logl=gompertz_logl(beta_,x,y)

  Input:
	beta_(1)   scale parameter of the Gompertz distribution = a
	beta_(2)   location parameter of the Gompertz distribution = b
	x          vector of the covariate points
	y          vector of the response values

  Output:
	logl       The value of the log-likelihood multiplied by -1



gpl.txt
  The GNU General Public License


optdesign_demo.m
  Demonstration of the use of the optdesign package

  Syntax: optdesign_demo	


readme.txt 
  This file.


thetaest.m
  Calculates estimates of theta and lambda, in complementary log-log regression

  Syntax: [thetaest,thetacov,theta_cf,lambda_cf]=thetaest(betaest,betacov) 

  Input:
	betaest		estimates of the parameters a and b
	betacov		covariance matrix for the parameters a and b

  Output:
	thetaest	estimates of parameters theta (the middle point) and lambda (the width of the curve)
	thetacov 	covariance matrix for the parameters theta and lambda
	theta_cf	95% confidence interval for the parameter theta
	lambda_cf	95% confidence interval for the parameter lambda




