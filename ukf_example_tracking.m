close all
clear all
clc
addLibrary('kalmanlib')
rng(101)
%% KALMAN FILTER 

Ts = 0.1; %sampling time

%x=[px py vx vy]';

%state matrix for 2D motion
A= [1 0 Ts 0;...
    0 1 0 Ts;...
    0 0 1 0;...
    0 0 0 1]
%input matrix for uniformly accelerated motion
B=[0.5*Ts*Ts 0;...
    0 0.5*Ts;...
    Ts 0;...
    0 Ts]

%radar in the origin (0,0)
xr=0;
yr=0;
%% dynamic model
f=@(x,u)A*x+B*u;

%% measurement model (y=h(x)+eps)
h=@(x) [
    sqrt((xr-x(1))^2+(yr-x(2))^2);...
    atan2(x(2)-yr,x(1)-xr);
    ];

%% define state covariance matrix (model accuracy)
Q=eye(4)*10e-4;

%define sensor covariance matrix (sensor quality)
R=0.01*eye(2);

N=100;
%% memory allocation
x=zeros(4,N);
hx=zeros(4,N);
u=zeros(2,N);
y=zeros(2,N);
P=zeros(4,4,N);
mse=zeros(1,N);
%%
%define initial state covariance matrix
P(:,:,1)=blkdiag(1e4,1e4,1,1);
%% UT parameters
UT.alpha=1e-3;
UT.beta=2; %Gaussian
UT.kappa=0; %kurtosis
%% simulate system
%define initial conditions
x(:,1)=[1;1;0;0]; %true state
hx(:,1)=zeros(4,1); %estimated state
mse(1)=norm(x(:,1)-hx(:,1))^2; %initial mse
trP(1)=trace(P(:,:,1)); %initial covariance trace
tic
for k=1:N
    %generate input signal
    u(:,k) = [1 1]';
    
    %generate true system signal
    x(:,k+1) = A*x(:,k)+B*u(:,k)+sqrt(Q)*randn(4,1);
    y(:,k+1) = h(x(:,k+1))+sqrt(R)*randn(2,1);
    
    %% Filtering (UKF)
    %Step 1 - Prediction
    [hx(:,k+1),P(:,:,k+1)] = UKF_prediction(f,hx(:,k),u(:,k),P(:,:,k),Q,UT);
    %Step 2 - Update    
    [hx(:,k+1),P(:,:,k+1)] = UKF_update(h,hx(:,k+1),P(:,:,k+1),R,y(:,k+1),UT);
    %% error
    mse(k+1)=norm(x(:,k+1)-hx(:,k+1))^2;
    trP(k+1)=trace(P(:,:,k+1));
end
time_elapsed=toc
%% show the results
figure
subplot(1,2,1)
plot(x(1,:),'g','linewidth',2)
hold on
plot(hx(1,:),'r','linewidth',2)
title('Distance in X')
legend('Real','UKF')
xlabel('Time (k)')
ylabel('x_k')
grid on
subplot(1,2,2)
plot(x(2,:),'g','linewidth',2)
hold on
plot(hx(2,:),'r','linewidth',2)
title('Distance in Y')
legend('Real','UKF')
xlabel('Time (k)')
ylabel('y_k')
grid on
figure
subplot(1,2,1)
plot(x(3,:),'g','linewidth',2)
hold on
plot(hx(3,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in X')
legend('Real','UKF')
xlabel('Time (k)')
ylabel('v^x_k')
grid on
subplot(1,2,2)
plot(x(4,:),'g','linewidth',2)
hold on
plot(hx(4,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in Y')
legend('Real','UKF')
xlabel('Time (k)')
ylabel('v_k^y')
grid on
figure
plot(x(1,:),x(2,:),'g','linewidth',2)
hold on
plot(hx(1,:),hx(2,:),'r-','linewidth',2)
grid on
legend('Real','UKF')
title('Trajectory')
xlabel('x_k')
ylabel('y_k')
%%
figure
plot(mse)
% Calculate the mean squared error over time
meanMSE = mean(mse);
disp(['Mean Squared Error: ', num2str(meanMSE)]);