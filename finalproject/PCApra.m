clear
close all

% 2D(m) dataset with 20 samples(n)
A = [2.89  0.32  5.8  -6.52  3.94 -4.21  0.45  2.14  1.3  -4.98 -2.4  -3.1 ...
   0.69 -1.59 -3.64 -0.24  6.81  4.63 -2.24 -0.06;1.52  0.91  1.52 -0.88 ...
   -0.03 -1.26 -0.25  0.96 -0.89 -0.45 -0.88 -1.12 -0.86  0.13 -1.53  0.51 ...
   2.66  1.28 -0.14 -1.19];
% centre 
A = A - mean(A(:));
% covariance matrix = A*A'/(# sampls-1) whose size is m
C = (A*A') / (size(A,2)-1);
% pricipal directions and components
[V,D] = eig(C);
% sort
[~,ind] = sort(diag(D),'descend');
V = V(:,ind);                       % negative or positive direction
% dataset
plot(A(1,:),A(2,:),'o')
x = -8:0.01:8;
% principal direction
y = V(2,1) .*x ./ V(1,1);
hold on;
plot(x,y,'r')
% project each point onto pricipal direction
P =  V(:,1)' * A .* V(:,1); 
for ii = 1:20
    hold on
    plot([A(1,ii) P(1,ii)],[A(2,ii) P(2,ii)],'b')
end
hold on
q = quiver(0,0,V(1,1),V(2,1),'black','linewidth',1.5);
legend(q,'$\vec{v}$ ','Interpreter','latex','Fontsize',14);
grid on
title('projection')
figure
subplot(2,1,1)
plot(A(1,:),A(2,:),'o')
title('original');
% transformation
X_T = V'*A;
subplot(2,1,2)
plot(X_T(1,:),X_T(2,:),'go')
ylim([-2 3]);
title('uncorrelated');
% SVD
[U1,S1,V1] = svd(A');
pc = (U1*S1)';