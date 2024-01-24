function p=gompertz_cdf(beta_,x)
%Calculates the cumulative distribution function of the Gompertz/Gumbel distribution
%p=gompertz_cdf(beta_,x)
%
%Input:
%beta_(1)   scale parameter of the Gompertz distribution = a
%beta_(2)   location parameter of the Gompertz distribution = b
%x          vector of points where the value of the cumulative
%           distribution function is calculated.
%Output:
%p          vector of the values of the cumulative distribution
%           function in points x.
%
%This function is a part of the optdesign package. See readme.txt for more info.
%Juha Karvanen 2005-03-28
a=beta_(1);
b=beta_(2);
p=1-exp(-exp(a*x+b));