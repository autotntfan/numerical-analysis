clear
close all

%% (a) random number generator using Box–Muller transform
N=1e6;
u1=rand(1,N);
u2=rand(1,N);
r=sqrt(-2.*log(u1));
noise1=r.*cos(2*pi*u2);
noise2=r.*sin(2*pi*u2);
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
    ylim=([0 18000]);
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
S_re=real(S)*ones(1,N);
S_im=imag(S)*ones(1,N);
S_theta=atan((S_im+noise2)./(S_re+noise1));
% regression
A=[ones(3,1) t]; 
X_hat=(A'*A)\A'*S_theta;
figure
histogram(X_hat(2,:))
u_omega=mean(X_hat(2,:));
sigma_omega=std(X_hat(2,:));
msg=sprintf('distribution of \\omega follows ~N(%.4f,%.4f)',[u_omega,sigma_omega]);
title(msg)
%% (d)
sigma_theta=std(S_theta,0,2);
mag=abs(S);
figure
plot(sigma_theta,mag,'ro')
A=[ones(3,1) sigma_theta];
X_hat=(A'*A)\A'*mag;
hold on
plot(sigma_theta,A*X_hat,'b')
W=[1 0 0;0 1 0;0 0 1];
A_prime=W*A;
X_hat=(A_prime'*A_prime)\A_prime'*(W*mag);
hold on
plot(sigma_theta,A_prime*X_hat,'go')
xlabel('\sigma_{theta}')
ylabel('|S(t)|')