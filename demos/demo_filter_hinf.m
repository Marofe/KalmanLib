close all
clear all
clc
addpath('../../kalmanlib')
%% Hinf Filter
n=2;
p=1;
%% system model
A=randn(n);
A=A*1/(1.01*max(abs(eig(A))));
abs(eig(A))
S=eye(p);
L=eye(p);
H=randn(p,n);
sigma=0.01*eye(n);
nu=0.1*eye(p);
Q=sigma*sigma';
R=nu*nu';
%% simulation parameters
N=100;
M=10;
%% simulation
for j=1:M
    %% parameters
    x(:,1)=ones(n,1);
    hx(:,1)=x(:,1);
    hx2(:,1)=x(:,1);
    P=10*eye(n);
    P2=eye(n);
    for k=1:N
        %% real system
        w(:,k)=sigma*ones(n,1)*cos(10*k)*(0.9)^k*randn;
        x(:,k+1)=A*x(:,k)+w(:,k);
        y(:,k+1)=H*x(:,k)+nu*randn(p,1);
        %% Kalman Filter
        [hx(:,k+1) P]=kalman_filter(A,H,Q,R,hx(:,k),P,y(:,k+1),0);
        mse_kf(j,k+1)=norm(x(:,k+1)-hx(:,k+1))^2;
        %% Hinf Filter
        lambda=5;
        [hx2(:,k+1) P2]=hinf_filter(A,H,L,S,lambda,Q,R,hx2(:,k),P2,y(:,k+1),0);
        mse_hinf(j,k+1)=norm(x(:,k+1)-hx2(:,k+1))^2;
    end
end
%% performance
mse_kf = mean(mse_kf);
mse_hinf = mean(mse_hinf);
gain = (1-sum(mse_hinf)/sum(mse_kf))*100;
disp([num2str(gain) '%'])
%% plot
p1=plot(10*log10(mse_kf));
hold on
p2=plot(10*log10(mse_hinf));
legend([p1 p2],'KF','Hinf')
grid on
xlabel('time (k)')
ylabel('MSE(k)_{dB}')
figure
subplot(2,1,1)
plot(w(1,:))
title('Disturb')
subplot(2,1,2)
autocorr(w(1,:))