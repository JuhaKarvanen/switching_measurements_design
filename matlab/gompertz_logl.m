function logl=gompertz_logl(beta_,x,y)
%Calculates the log-likehood for complementary log-log regression
%logl=gompertz_logl(beta_,x,y)
%
%Input:
%beta_(1)   scale parameter of the Gompertz distribution = a
%beta_(2)   location parameter of the Gompertz distribution = b
%x          vector of the covariate points
%y          vector of the response values
%
%Output:
%logl       The value of the log-likelihood multiplied by -1
%
%This function is a part of the optdesign package. See readme.txt for more info. 
%Juha Karvanen 2005-04-09
Fx=max(0.00000001,min(0.99999999,gompertz_cdf(beta_,x)));
logl=sum(y.*log(Fx)+(1-y).*log(1-Fx));
logl=-logl;