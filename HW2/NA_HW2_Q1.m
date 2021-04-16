% HW2 Q1 
clear;
close all;
% (a) Jacobian method
% initial condition
u=0;
v=0;
count=1;
epsilon_u=1;
epsilon_v=1;
% stopping criterion
epsilon_t=0.001;
% Exact solution,let A=[47 -25;-25 70] b=[0 100],and determine x so that Ax=b
A=[47 -25;-25 70];b=[0;100];
x=A\b;
% recursive function
while epsilon_u>=epsilon_t & epsilon_v>=epsilon_t
    u(count+1)=25*v(count)/47;
    v(count+1)=(100+25*u(count))/70;
    epsilon_u(count+1)=abs(u(count+1)-x(1))/x(1);
    epsilon_v(count+1)=abs(v(count+1)-x(2))/x(2);
    count=count+1;
end
figure('Name','Jacobi')
subplot(2,1,1)
plot(linspace(0,count,count),u,'b');
hold on;
plot(linspace(0,count,count),v,'r');
legend('$i_{1} $','$i_{2} $','Interpreter','latex','FontSize',15)
xlabel('iteration')
ylabel('solution')
title('solution v.s. iteration');

subplot(2,1,2)
plot(linspace(0,count,count),epsilon_v*100,'b');
hold on;
plot(linspace(0,count,count),epsilon_u*100,'r');
legend('$i_{1} $','$i_{2} $','Interpreter','latex','FontSize',15)
xlabel('iteration')
ylabel('absolute error(%)')
title('absolute error v.s. iteration');

% (b) Gauss-Seidel method
% initial condition
u=0;
v=0;
count=1;
epsilon_u=1;
epsilon_v=1;
% stopping criterion
epsilon_t=0.001;
while epsilon_u>=epsilon_t & epsilon_v>=epsilon_t
    u(count+1)=25*v(count)/47;
    % u(count)->u(count+1)
    v(count+1)=(100+25*u(count+1))/70;
    epsilon_u(count+1)=abs(u(count+1)-x(1))/x(1);
    epsilon_v(count+1)=abs(v(count+1)-x(2))/x(2);
    count=count+1;
end
figure('Name','Gauss-Seidal')
subplot(2,1,1)
plot(linspace(0,count,count),u,'b');
hold on;
plot(linspace(0,count,count),v,'r');
legend('$i_{1} $','$i_{2} $','Interpreter','latex','FontSize',15)
xlabel('iteration')
ylabel('solution')
title('solution v.s. iteration');

subplot(2,1,2)
plot(linspace(0,count,count),epsilon_v*100,'b');
hold on;
plot(linspace(0,count,count),epsilon_u*100,'r');
legend('$i_{1} $','$i_{2} $','Interpreter','latex','FontSize',15)
xlabel('iteration')
ylabel('absolute error(%)')
title('absolute error v.s. iteration');

% (d) redo Jacobi method and G-S method with relaxation
colors=['r','g','b'];
lambdas=[0.3 0.7 1.2];
for ii =1:length(lambdas)
    % initial condition
    u=0;
    v=0;
    count=1;
    epsilon_u=1;
    epsilon_v=1;
    % stopping criterion
    epsilon_t=0.001;
    while epsilon_u>=epsilon_t & epsilon_v>=epsilon_t
        u(count+1)=lambdas(ii)*(25*v(count)/47)+(1-lambdas(ii))*u(count);
        v(count+1)=lambdas(ii)*(100+25*u(count))/70+(1-lambdas(ii))*v(count);
        epsilon_u(count+1)=abs(u(count+1)-x(1))/x(1);
        epsilon_v(count+1)=abs(v(count+1)-x(2))/x(2);
        count=count+1;
    end
    figure(4)
    hold on;
    plot(linspace(0,count,count),epsilon_v*100,colors(ii));
    xlabel('iteration')
    ylabel('absolute error(%)')
    title('Jacobi absolute error v.s. iteration ($i_{1} $)','Interpreter','latex','FontSize',15);
    figure(5)
    hold on;
    plot(linspace(0,count,count),epsilon_u*100,colors(ii));
    xlabel('iteration')
    ylabel('absolute error(%)')
    title('Jacobi absolute error v.s. iteration  ($i_{2} $)','Interpreter','latex','FontSize',15);
    % initial condition
    u=0;
    v=0;
    count=1;
    epsilon_u=1;
    epsilon_v=1;
    % stopping criterion
    epsilon_t=0.001;
    while epsilon_u>=epsilon_t & epsilon_v>=epsilon_t
        u(count+1)=lambdas(ii)*(25*v(count)/47)+(1-lambdas(ii))*u(count);
        % G-S:u(count)->u(count+1)
        v(count+1)=lambdas(ii)*(100+25*u(count+1))/70+(1-lambdas(ii))*v(count);
        epsilon_u(count+1)=abs(u(count+1)-x(1))/x(1);
        epsilon_v(count+1)=abs(v(count+1)-x(2))/x(2);
        count=count+1;
    end
    figure(6)
    hold on;
    plot(linspace(0,count,count),epsilon_v*100,colors(ii));
    xlabel('iteration')
    ylabel('absolute error(%)')
    title('Gauss-Seidel absolute error v.s. iteration ($i_{1} $)','Interpreter','latex','FontSize',15);
    figure(7)
    hold on;
    plot(linspace(0,count,count),epsilon_u*100,colors(ii));
    xlabel('iteration')
    ylabel('absolute error(%)')
    title('Gauss-Seidel absolute error v.s. iteration  ($i_{2} $)','Interpreter','latex','FontSize',15);
end
figure(4)
legend('\lambda=0.3 ','\lambda=0.7 ','\lambda=1.2 ');
figure(5)
legend('\lambda=0.3 ','\lambda=0.7 ','\lambda=1.2 ');
figure(6)
legend('\lambda=0.3 ','\lambda=0.7 ','\lambda=1.2 ');
figure(7)
legend('\lambda=0.3 ','\lambda=0.7 ','\lambda=1.2 ');


