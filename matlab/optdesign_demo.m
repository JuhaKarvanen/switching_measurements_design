%optdesign_demo.m - Demonstation of the optdesign package
%This function is a part of the optdesign package. See readme.txt for more info. 
%Juha Karvanen 2005-09-28
fprintf(['\nThis is a demonstration of the optdesign package.\n']);

clear all
a_true=0.24;
b_true=-61;
z1=0.979633;
z2=-1.33774;
nstage=40;
xplot=230:0.1:270;

optx=zeros(2,nstage);
betas=zeros(2,nstage);
betacov=zeros(2,2,nstage);
thetas=zeros(2,nstage);
thetacov=zeros(2,2,nstage);
theta_cf=zeros(2,nstage);
lambda_cf=zeros(2,nstage);
ismles=-ones(nstage);
totaln=zeros(nstage);
optx_c=zeros(2,1);
beta_c=zeros(2,1);
minlimit=200;
maxlimit=300;

fprintf('Consider a measurement where the binary response Y is a random variable \n')
fprintf('that follows the complementary log-log regression model\n\n')
fprintf('P(Y=1) = 1-exp(-exp(a*x+b)),\n\n')
fprintf('where a and b are the unknown model parameters, x is the covariate variable.\n\n')
fprintf('In the measurement, we choose the point of measurement x and measure the response Y.\n')
fprintf('Our objective is the estimate the parameters a and b in optimal manner (D-optimality).\n\n')
fprintf('At the beginning we know only that the response curve rises from 0.01 to 0.99 \n')
fprintf('somewhere in the interval [%g %g].\n',minlimit,maxlimit)
fprintf('In the following we apply sequential design %d of stages to estimate \n' ,nstage)
fprintf('the parameters of the response curve\n')
fprintf('\nPress any key to continue...\n');
pause;

x=[minlimit; maxlimit];
y=rand(2,1)<gompertz_cdf([a_true b_true],x);
n=2;
ismle=0;
ismle_old=0;
xpoint=x;
ypoint=y;
figure

for i=1:nstage,
      if n>2 
          [optx_c beta_c betacov_c ismle]=optdesign(x,y,beta_c);
      else 
          [optx_c beta_c betacov_c ismle]=optdesign(x,y);
      end
      optx(:,i)=optx_c';
      betas(:,i)=beta_c';
      betacov(:,:,i)=betacov_c;
      ismles(:,i)=ismle;
      if ismle==0, 
         if n==2, n_inc=23;
         else n_inc=25;
         end
         optx_c=min(maxlimit,max(minlimit,optx_c));
         x_inc=repmat(optx_c(1),n_inc,1);
         y_inc=rand(n_inc,1)<gompertz_cdf([a_true b_true],x_inc);
         n_inc=50;
         xpoint=[xpoint;optx_c(1)];
         ypoint=[ypoint;mean(y_inc)];
      else 
         [thetas(:,i),thetacov(:,:,i),theta_cf(:,i),lambda_cf(:,i)]=thetaest(beta_c,betacov_c);
         n_inc=round(1.1*n_inc)+mod(round(1.1*n_inc),2);
         optx_c=min(maxlimit,max(minlimit,optx_c));
         x_inc=repmat(optx_c',n_inc/2,1);
         y_inc=rand(n_inc,1)<gompertz_cdf([a_true b_true],x_inc);
         xpoint=[xpoint;optx_c(1)];
         ypoint=[ypoint;mean(y_inc(1:2:(n_inc-1)))];
         xpoint=[xpoint;optx_c(2)];
         ypoint=[ypoint;mean(y_inc(2:2:(n_inc)))];
      end
      
        fprintf('\nStage: %d, Total n: %d\n',i,n)
        fprintf('Point(s) of measurement: %g %g\n',optx_c(1),optx_c(2))
        if ismle==1 & ismle_old==0
            fprintf('Maximum likelihood estimates found!\n')
        end
        if ismle==1
            fprintf('Parameter  Estimate        std\n')
            fprintf('a           %1.5f       %1.5f\n',betas(1,i),sqrt(betacov(1,1,i)))
            fprintf('b          %1.4f       %1.5f\n',betas(2,i),sqrt(betacov(2,2,i)))
            fprintf('theta       %1.3f       %1.5f\n',thetas(1,i),sqrt(thetacov(1,1,i)))
            fprintf('lambda      %1.4f       %1.5f\n',thetas(2,i),sqrt(thetacov(2,2,i)))
        end
            
        lwidth=1.5;
        msize=8;
        fsize=14;
        xpointlast=length(xpoint);
        plottitle=['Stage ',num2str(i),': Response curve'];
        clf
        plot(xplot,gompertz_cdf([a_true b_true],xplot),'k',...
            xplot,gompertz_cdf(beta_c,xplot),'b',...
            xpoint,ypoint,'kx','MarkerSize',msize)
        hold on
        plot(xpoint(xpointlast-1:xpointlast),ypoint(xpointlast-1:xpointlast),'rx','MarkerSize',msize);
        title(plottitle,'FontSize',fsize)
        legend('True ','Estimated','Measured point',2)
        legend('boxoff');
        set(gca, 'Ylim',[0,1]);
        set(gca, 'Xlim',[220,280]);
        xlabel('x','FontSize',fsize);
        ylabel('P','FontSize',fsize,'Rotation',1);
        set(gca,'LineWidth',lwidth);
        set(gca,'FontSize',fsize);
        h=findobj('Type','line');
        set(h,'LineWidth',lwidth);
        
        pause(min(2.5*(i<15),i));
     
        ismle_old=ismle;
        n=n+n_inc;
   		x=[x;x_inc];
      	y=[y;y_inc];
     
end

 
