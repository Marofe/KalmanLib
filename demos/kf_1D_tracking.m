close all
clear all
clc
addpath('../kalmanlib')
addLibrary('../kalmanlib')
rng(101) %fix the seed of the random generator for reproducibility
%% System
% 1D tracking example 
% x=[p v] -> state vector
% y=p -> measurement vector
n=2; %dimension of state vector
m=1; %dimension of measurement vector
%% linear dynamic model
dt=1; %sample time
A=[1 dt;0 1];
G=[0;1]; 
H=[1 0];
%% Noise Properties
Q=.1^2; %process noise cov (velocity fluctuation)
stdQ=sqrt(Q); %process noise std
R=2^2; %measurement noise cov (radar noise)
stdR=sqrt(R); %measurement noise std
%% Time-horizon
N=50;
time=0:dt:(N-1)*dt;
%% Allocate memory
x=zeros(n,N);
y=zeros(m,N);
hx=zeros(n,N);
P=zeros(n,n,N);
%% Initial condition
x(:,1)=[90 3]'; %true
hx(:,1)=[100 5]'; %initial guess
P(:,:,1)=diag([5^2,1^2]);
sys=struct('A',A,'H',H,'Q',G*Q*G','R',R);
%% Simulation
    for k=1:N-1
        %% true system
        x(:,k+1)=A*x(:,k)+G*stdQ*randn;
        y(:,k+1)=H*x(:,k+1)+stdR*randn(m,1);
        %% Kalman Filter
        [hx(:,k+1) P(:,:,k+1)]=kalmanFilter(hx(:,k),P(:,:,k),y(:,k+1),sys);
    end
%% plot
figure
subplot(2,1,1)
p1=plot(time,x(1,:),'linewidth',1.5);
hold on
p2=plot(time,hx(1,:),'linewidth',1.5);
stdPos=squeeze(P(1,1,:))';
ciY=[hx(1,:)+3*stdPos;hx(1,:)-3*stdPos];
plot(time,ciY','color',colors('red'))
p3=plotConfidence(time,ciY,colors('red'));
legend([p1,p2,p3],'Ground-truth','Kalman Filter','3\sigma')
grid on
title('Position')
xlabel('time')
subplot(2,1,2)
p1=plot(time,x(2,:),'linewidth',1.5);
hold on
p2=plot(time,hx(2,:),'linewidth',1.5);
stdVel=squeeze(P(2,2,:))';
ciY=[hx(2,:)+3*stdVel;hx(2,:)-3*stdVel];
plot(time,ciY','color',colors('red'))
p3=plotConfidence(time,ciY,colors('red'));
legend([p1,p2,p3],'Ground-truth','Kalman Filter','3\sigma')
grid on
title('Velocity')
xlabel('time')