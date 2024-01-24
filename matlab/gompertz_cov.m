function lcov=gompertz_cov(beta_,x,y)
%Calculates the covariance matrix for parameters a and b in complementary log-log regression
%lcov=gompertz_cov(beta_,x,y)
%
%Input:
%beta_(1)   scale parameter of the Gompertz distribution = a
%beta_(2)   location parameter of the Gompertz distribution = b
%x          vector of the covariate points
%y          vector of the response values
%
%Output:
%lcov       The (asymptotic) covariance matrix for beta_ 
%
%This function is a part of the optdesign package. See readme.txt for more info. 
%Juha Karvanen 2005-05-02
a=beta_(1);
b=beta_(2);
z=a*x+b;
expz=exp(z);
%In order to avoid numerical problems, very small values are truncated to
%zero.
zind=((z>-10) & (z<5)); 
expz_=expz(zind);
d2=zeros(size(z));
d2(zind)=(-expz_+exp(z(zind)+expz_)-exp(2*z(zind)+expz_))./(-1+exp(expz_)).^2;
%The covariance matrix is calculated as the inverse of the observed Fisher
%information.
J(1,1)=-sum(y.*(x.^2).*d2-(1-y).*(x.^2).*expz);
J(1,2)=-sum(y.*x.*d2-(1-y).*x.*expz);
J(2,1)=J(1,2);
J(2,2)=-sum(y.*d2-(1-y).*expz);
lcov=inv(J);
