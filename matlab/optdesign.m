function [optx,betaest,betacov,IsMLE]=optdesign(xx,yy,beta0,transform)
%[optx,betaest,betacov]=optdesign(xx,yy)  - Calculates the optimal design for a switching measurement
%Version 1.1 2006-05-06
%Version 1.0 2005-10-01
%Copyright Juha Karvanen 2005
%The program is distributed under the terms of the GNU General Public License
%For details of the copyright see the files readme.txt and gpl.txt provided in the same package.
%
%The program calculates the optimal design for swithching measurements in 
%superconducting Josephson junctions. More generally, the program calculates 
%the optimal covariate points for binary response with complementary
%log-log link. The optimality is defined as D-optimality and parameters are estimated 
%by maximum likelihood estimation. The data is assumed to follow
%Gompertz/Gumbel/extreme value distribution
%
%               P(Y=1) = 1-exp(-exp(a*x+b)),
%
%where a and b are the model parameters, x is the covariate variable and Y is the binary 
%response. In the case of one Josephson junction this model arises from the laws of
%physics. All data collected so far is given as an input to the program which then returns 
%the optimal covariate points to be measured in next stage. 
%
%Sequential design with optdesign.m: (demonstrated in optdesign_demo.m)
%
%	1) Obtain the limits of the initial interval: xmin and xmax
%	2) Measure at points xmin and xmax. Let x contain the points
%	   of measurement and y contain the measured responses. 
%	   (Sometimes it can be directly assumed that the response for 
%	    xmin is 0 and the response for xmax is 1.)
%	3) Call opdesign.m
%	   [optx,betaest,betacov,IsMLE]=optdesign(x,y)
%	4) If the required accuracy of betaest has been achieved, stop the 
%	   procedure; otherwise proceed to the next step.
%	5) Measure at the points given by optx. Append the points of measurement 
%	   to x. Append the measured responses to y.
%	6) Go to step 3.
%
%The function is intended to be used for sequential designs of switching 
%measurements but it may used for other suitable purposes as well. 
%
%
%Reference:     Juha Karvanen, Juha J. Vartiainen, Andrey Timofeev and Jukka
%               Pekola, Experimental Designs for Binary Data in Switching
%               Measurements on Superconducting Josephson Junctions, 2005,
%               http://ltl.tkk.fi/PICO/optdesign/
%
%Uses:          gompertz_logl.m
%               gompertz_cov.m
%               gompertz_cdf (through gompertz_logl.m)
%
%Input:
%xx             Points (pulse heights) measured so far, 
%                   size(xx)=[n 1], where n is the total number of points
%                   measured so far.
%yy             Measured binary response for points xx, 0 or 1, size(yy)=[n 1].
%beta0          (optional) initial values of a and b for maximum likelihood
%                   estimation, beta0=[a0 b0], lenght(beta0)=2.
%transform		(Not in use in this version)	
%
%Output:
%optx			Points (pulse heights) to be measured at the next stage, length(optx)=2
%                   If IsMLE==1 measure at both optx(1) and optx(2)
%                   If IsMLE==0 measure only at optx(1).
%betaest		Estimated values of parameters a and b
%                   If IsMLE==0, an ad-hoc estimate is returned.
%betacov        Estimated covariance matrix for parameters a and b,
%                   size(betacov)=[2,2]
%                   If IsMLE==0, betacov cannot be calculated and a matrix of NaN:s is returned.
%IsMLE          Are the estimates maximum likelihood estimates? 0 or 1.
%
%Enter optdesign_demo() to run a demonstration of use.


% Checking input %
if nargin<2
    error('Input variables xx and yy are required.')
end
[nx px]=size(xx);
[ny py]=size(yy);
if px~=1 | py~=1 | nx~=ny
    error('Input variables xx and yy should be column vectors of same size: size(xx)==size(yy)==[n 1]')
end
if sum(yy~=0 & yy~=1)>0 
    error('Input variable yy may contain only values 0 and 1.')
end
if not(isreal(xx))
    error('Input variable xx must be real valued.')
end
if (var(xx)==0)
    error('Input variable xx must have at least two distinct values.')
end    
if nargin>=3 %Checking for beta0
    if not(isreal(beta0)) | length(beta0)~=2 
        error('Input variable beta0 must be a real valued vector of two values [a0>0 b0].')
    end
    if beta0(1)<=0 
        error('Input variable beta0 must be a real valued vector of two values [a0>0 b0].')
    end
end

% MATLAB version
verz=version;
% Initialization
betaest=[NaN; NaN];
betacov=[NaN NaN; NaN NaN];
IsMLE=0;
epsilon=1;

[minx minxind]=min(xx);
[maxx maxxind]=max(xx);
if (mean(yy(xx==minx))>=mean(yy(xx==maxx))) 
    warning('Assumption a>0 violated. Check the data.')
end
%Checking the existence of maximum likelihood estimate
[zeromax maxind]=max(xx(yy==0,1));
[onemin minind]=min(xx(yy==1,1));
if zeromax<onemin,
    x1=(zeromax+onemin)/2;
    optx=[x1 x1];
    %When MLE does not exist, an ad-hoc estimate is returned.
    a=9.6/(maxx-minx); %9.6=-log(-log(0.9995))+log(-log(0.0005))
    b=log(-log(0.5))-a*(zeromax+onemin)/2;
	betaest=[a b];
elseif onemin==zeromax, 
    %If onemin same as zeromax, one of the points needed is found and the other is searched in the neighbourhood 
    epsilonsign=sign(0.5-mean(yy));
    if (epsilonsign==0) 
        epsilonsign=2*(rand(1)-0.5)-1; %randomly 1 or -1
    end
    if (epsilonsign==1) 
        epsilon=(min(xx(xx>zeromax | xx==maxx))-zeromax)/2;
    end
    if (epsilonsign==-1) 
        epsilon=(max(xx(xx<zeromax | xx==minx))-zeromax)/2;
    end
    x1=zeromax+epsilon;
    optx=[x1 x1];
   
    %When MLE does not exist, an ad-hoc estimate is returned.
    a=9.6/(maxx-minx); %9.6=-log(-log(0.9995))+log(-log(0.0005))
    b=log(-log(0.5))-a*(zeromax+onemin)/2;
	betaest=[a b];
else
    IsMLE=1;
    if nargin==2 %No beta0 is given as input, an ad-hoc beta0 is calculated.
        a=9.6/(maxx-minx); %9.6=-log(-log(0.9995))+log(-log(0.0005))
        b=log(-log(0.5))-a*(zeromax+onemin)/2;
		beta0=[a b];
	end
	x=xx;
	y=yy;
    %The optimal points (canonical form)
	z2=0.979633;
    z1=-1.33774;
      
    %The name of optimization procedure is different in the earlier
    %MATLAB versions. We use here the Nelder-Mead optimization that slower
    %than the gradient based methods but avoids the risk of the numerical
    %instability of gradient. In our application (switching measurements) 
    %the stability is critical.
    if str2num(verz(1))<=5, 
        fminopt=FOPTIONS;
        fminopt(14)=4000;
        betaest=fmins('gompertz_logl',beta0,fminopt,[],x,y);
    else
        fminopt=optimset('MaxIter',4000);
        betaest=fminsearch('gompertz_logl',beta0,fminopt,x,y);
    end
 	a=betaest(1);
	b=betaest(2);
	x1=(z1-b)/a;
	x2=(z2-b)/a;
	optx=[x1 x2];
    betacov=gompertz_cov(betaest,x,y);
end


