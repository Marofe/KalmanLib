close all
clear all
clc
addpath('kalmanlib')
%% Demo Hinf Filter
n=2; %number of states
p=1; %number of outputs

%% System Matrices
A0=[0.9802 0.0196;0 0.9802];
A1=[0 0.099;0 0];
H=[1 -1];
G=eye(n);
L=eye(n);
Ef=[0 5];
Eg=[0 0];
M=[0.0198;0];

%% Noise Covariances
sigma=sqrtm([1.9608 0.0195;0.0195 1.9605]);
nu=eye(p);
Q = sigma*sigma';
R = nu*nu';
%% Hinf parameters
gamma=100;
%%
N=1000;
for j=1:100
    delta(j)=(1-2*rand);
    A=A0+delta(j)*A1;
 clear x hx 
x(:,1)=1*ones(n,1);
hx(:,1)=zeros(n,1);
hx2=hx;
P=eye(n);
P2=P;
P0=P;
hx0=hx;
trP(1)=trace(P);
trP2=trP;
trP0=trP;
e_kf(j,1)=norm(x(:,1)-hx(:,1))^2;
e_hinf(j,1)=norm(x(:,1)-hx2(:,1))^2;
e_opt(j,1)=norm(x(:,1)-hx0(:,1))^2;
for k=1:N-1
    %% real system
    x(:,k+1)=A*x(:,k)+sigma*randn(n,1);
    y(:,k+1)=H*x(:,k+1)+nu*randn(p,1);
    %% Kalman Filter
    [hx0(:,k+1),P0]=kalman_filter(A0,H,Q,R,hx0(:,k),P0,y(:,k+1));
    %% Kalman Filter with Uncertainty
    [hx(:,k+1),P]=kalman_filter(A,H,Q,R,hx(:,k),P,y(:,k+1));
    %% Hinf Filter
    [hx2(:,k+1),P2]=hinffilter(A0,H,G,Q,R,hx2(:,k),P2,y(:,k+1),L,gamma);
    %% erro
    e_kf(j,k+1)=norm(x(:,k+1)-hx(:,k+1))^2;
    e_hinf(j,k+1)=norm(x(:,k+1)-hx2(:,k+1))^2;
    e_opt(j,k+1)=norm(x(:,k+1)-hx0(:,k+1))^2;
    trP(k+1)=trace(P);
    trP2(k+1)=trace(P2);
    trP0(k+1)=trace(P0);
end
end
%% show 
close all
plot(trP,'linewidth',2)
hold on
plot(trP0,'--k','linewidth',2)
plot(trP2,'linewidth',2)
title('Trace')
figure
semilogx(10*log10(mean(e_kf)),'linewidth',2)
hold on
semilogx(10*log10(mean(e_hinf)),'linewidth',2)
semilogx(10*log10(mean(e_opt)),'linewidth',2)
legend('KF','Hinf','Optimal')
grid on
figure
plot(hx')
hold on
plot(hx2','--')
figure
histogram(delta)