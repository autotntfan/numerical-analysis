% HW2 Q2
clear
close all

% meticulousrun()
load epsilon.mat
figure
plot(epsilon)
ylabel('absolute error(%)');
xlabel('N(x1000)');
title('meticulous result in multiples of 1000');
load epsilon1.mat
% algorithm2()
load epsilon1.mat
figure

plot(epsilon1);
xlabel('N')
ylabel('absolute error(%)')
title('variant repeated time')

N=linspace(1,10000,100);
num=1000;
epsilon2=zeros(1,length(N));
for jj = 1:length(N)
    n=N(jj);
    pi_prob=zeros(1,num);
    for ii=1:num
        uniform_x=2*rand(1,n)-1;
        uniform_y=2*rand(1,n)-1;
        uniform_r=sqrt(uniform_x.^2+uniform_y.^2);
        pi_prob(ii)=4*sum(uniform_r<=1)/n;
    end
    avg_pi=mean(pi_prob);
    epsilon2(jj)=100*abs(avg_pi-pi)/pi;
end
figure
plot(N,epsilon2);
xlabel('N')
ylabel('absolute error(%)')
title(['repeat ' num2str(num) ' times'])      

function algorithm2()
    N=linspace(10,10000,1000);
    proudctnum=10000000;
    epsilon1=zeros(1,length(N));
    for jj = 1:length(N)
        n=N(jj);
        num=round(proudctnum/n);
        pi_prob=zeros(1,num);
        for ii=1:num
            uniform_x=2*rand(1,n)-1;
            uniform_y=2*rand(1,n)-1;
            uniform_r=sqrt(uniform_x.^2+uniform_y.^2);
            pi_prob(ii)=4*sum(uniform_r<=1)/n;
        end
        avg_pi=mean(pi_prob);
        epsilon1(jj)=100*abs(avg_pi-pi)/pi;
    end
    save('epsilon1.mat','epsilon1');
end

function meticulousrun()
    format long;
    N=[1000,5000,10000,1e5,1e6,1e7];
    % To "parfor" work successfully
    n_space=1000:1000:max(N);
    x=linspace(-1,1,max(N));
    y_plus=sqrt(1-x.^2);
    y_minus=-sqrt(1-x.^2);
    epsilon=ones(1,max(N)/1000);
    parfor ii=1:max(N)/1000
        n=n_space(ii);
        uniform_x=2*rand(1,n)-1;
        uniform_y=2*rand(1,n)-1;
        [logic,loc]=ismember(n,N);
        uniform_r=sqrt(uniform_x.^2+uniform_y.^2);
        pi_prob=4*sum(uniform_r<=1)/n;
        if logic==true
            figure(loc)
            plot(uniform_x,uniform_y,'b.');
            hold on
            plot(x,y_plus,'r');
            hold on
            plot(x,y_minus,'r');
            title(['N = ' num2str(n) '       \pi \approx ' num2str(pi_prob,16)]);
            xlabel('x');
            ylabel('y');
            saveas(gcf,['scatter' num2str(n) '.jpg']);
        end
        epsilon(ii)=100*abs(pi_prob-pi)/pi;
    end
    save('epsilon.mat','epsilon');
end




    






