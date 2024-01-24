function [thetaest,thetacov,theta_cf,lambda_cf]=thetaest(betaest,betacov)
%Calculates estimates of theta and lambda, in complementary log-log regression
%[thetaest,thetacov,theta_cf,lambda_cf]=thetaest(betaest,betacov) 
%
%Here Theta stand for the 50% point of the response curve, while 
%lambda is the difference of the 10% and 90% points, thus characterizing
%the width of the response curve.
%
%Input:
%	betaest		estimates of the parameters a and b
%	betacov		covariance matrix for the parameters a and b
%
%Output:
%	thetaest	estimates of parameters theta (the middle point) and lambda (the width of the curve)
%	thetacov 	covariance matrix for the parameters theta and lambda
%	theta_cf	95% confidence interval for the parameter theta
%	lambda_cf	95% confidence interval for the parameter lambda
%
%This function is a part of the optdesign package. See readme.txt for more info. 
%Juha Karvanen 2005-05-09
q1=log(-log(0.5));
q2=log(-log(0.1))-log(-log(0.9));
a=betaest(1);
b=betaest(2);
theta=(q1-b)/a;
lambda=q2/a;
thetaest=[theta lambda];
A=[(b-q1)/(a^2) -1/a;-q2/(a^2) 0];
thetacov=A*betacov*A';
stds=sqrt(diag(thetacov));
%Normal approximation is used for the confidence intervals. 
theta_cf=[theta-1.96*stds(1) theta+1.96*stds(1)];
lambda_cf=[lambda-1.96*stds(2) lambda+1.96*stds(2)];

thetaest=thetaest';
thetacov=thetacov';
theta_cf=theta_cf';
lambda_cf=lambda_cf';