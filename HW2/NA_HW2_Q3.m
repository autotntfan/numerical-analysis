clear
close all

%% (a) random number generator using Box–Muller transform
N=1e6;
u1=rand(1,N);
u2=rand(1,N);
r=sqrt(-2.*log(u1));
noise1=r.*cos(2*pi*u2);
noise2=r.*sin(2*pi*u2);
load noise1.mat
load noise2.mat
histogram(noise1,'BinWidth',0.1,'Normalization','pdf')
ylim([0 0.45])
ylabel('probability')
title('myself generator')
hold on
% verify RNG with pdf of N(0,1)
x=-5:0.1:5;
y=exp(-x.^2./2)./sqrt(2*pi);
plot(x,y,'r');
legend('generator','N(0,1)')
%% (b) 
r=200;theta=pi/4;u_theta=zeros(1,5);sigma_theta=zeros(1,5);
signal_re=r*cos(theta).*ones(1,N);
signal_im=r*sin(theta).*ones(1,N);
new_theta=atan((signal_im+noise2)./(signal_re+noise1));
figure
histogram(new_theta)
u_theta(1)=mean(new_theta);
sigma_theta(1)=std(new_theta);
msg=sprintf('distribution of \\theta follows ~N(%.4f,%.4f)',[u_theta(1),sigma_theta(1)]);
title(msg)
ylabel('count')

%% (c)
figure
r=[80 40 20 10];theta=pi/4;
for ii = 1:length(r)
    signal_re=r(ii)*cos(theta).*ones(1,N);
    signal_im=r(ii)*sin(theta).*ones(1,N);
    new_theta=atan((signal_im+noise2)./(signal_re+noise1));
    % in the report, I drew each figures in individual windows to clearly identify them. 
    subplot(2,2,ii)
    histogram(new_theta)
    u_theta(1+ii)=mean(new_theta);
    sigma_theta(1+ii)=std(new_theta);
    msg=sprintf('distribution of \\theta follows ~N(%.4f,%.4f) as r = %d',[u_theta(1+ii),sigma_theta(1+ii),r(ii)]);
    title(msg)
    ylabel('count')
    set(gca,'XLim',[0 1.3],'YLim',[0 18000])
end

%% (d)
figure
r=[200 r];
plot(r,sigma_theta,'o')
xlabel('r')
ylabel('\sigma_{\theta}','Fontsize',16)
hold on
plot(r,1./r,'r')
legend('$\sigma_{\theta}$','$\sigma_{\theta}=\frac{1}{r}$','interpreter','latex','Fontsize',14)

%% (e)
t=[10 50 90]';
S0=100;T2=32;w=pi/200;theta0=pi/4;
S=S0*exp(-t./T2+1i*(w*t+theta0));
S_re=real(S)*ones(1,N)+noise1;  % real signal component with noise
S_im=imag(S)*ones(1,N)+noise2;  % imaginary component with noise
% atan(X) returns values in the interval [-π/2, π/2],so we need to change
% it to [0, 2π]
S_theta=atan(S_im./S_re);
S_theta=S_theta+pi*(S_re<0 & S_im>0); % II quadrant:[-π/2,0]->[π/2,π]
S_theta=S_theta-pi*(S_re<0 & S_im<0); % III quadrant:[0,π/2]->[-π,-π/2]
S_theta=2*pi*(S_theta<=0)+S_theta;    % [-π π]->[0 2π]
% regression: least square solution:Ax=b=>A'Ax=A'b=>x=inv(A'A)A'b
A=[ones(3,1) t];         % A=[1 t1;1 t2;1 t3]=[1 10;1 50;1 90]
X_hat=(A'*A)\A'*S_theta; % least square solution X=[a0 a1]=> y=a0+a1x=XA
w=X_hat(2,:);            % a1=w(i.e. omega)
figure
H=histogram(w,'Binwidth',1e-4);
set(gca,'YLim',[0 3.5*1e4]);
% find the most probable omega to evaluate our model
omega=findomega(H);
% find each mean of the N group
u_omega=mean(w);
sigma_omega=std(w);
msg=sprintf('distribution of  \\omega  follows ~N(%.4f,%.4f)',[u_omega,sigma_omega]);
title(msg)
% plot the true omega versus our finding omega
figure
plot(t,pi/4+pi*t./200,'r')
hold on
plot(t,pi/4+omega*t,'b--')
legend('true \omega','our \omega')
% Do NOT run the following code! Matlab will break down!
% if you want to check it,please reduce N at line five,such as N=1e4,then rerun the kernel.
%%% figure
%%% plot(t(1),S_theta(1,:),'o',t(2),S_theta(2,:),'o',t(3),S_theta(3,:),'o')
%% (d)
W=[0.6 0 0;0 0.3 0;0 0 0.1];
A_prime=W*A;
X_hat=(A_prime'*A_prime)\A_prime'*(W*S_theta);
E1=W*S_theta-A_prime*X_hat;
w=X_hat(2,:);
figure
H=histogram(w,'Binwidth',1e-4);
u_omega=mean(w);
sigma_omega=std(w);
msg=sprintf('distribution of  \\omega  follows ~N(%.4f,%.4f)',[u_omega,sigma_omega]);
title(msg)
figure
plot(t,pi/4+pi*t./200,'r')
hold on
plot(t,pi/4+omega*t,'b--')
hold on
omega=findomega(H);
plot(t,pi/4+omega*t,'g--')
legend('true \omega','our \omega','wighted \omega')
function omega=findomega(fig)
    [~,loc]=max(fig.Values);
    omega=fig.BinEdges(loc)+5*1e-5;
end
