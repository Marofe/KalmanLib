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

%radar position
xr=0;
yr=0;
%% dynamic model
f=@(x,u)A*x+B*u;
Jf=@(x)A; %jacobian matrix of the dynamical model
%% measurement model (y=h(x)+eps)
h=@(x) [
    sqrt((xr-x(1))^2+(yr-x(2))^2);...
    atan2(x(2)-yr,x(1)-xr);
    ];
Jh=@(x)[(x(1)-xr)/sqrt((xr-x(1))^2+(yr-x(2))^2) (x(2)-yr)/sqrt((xr-x(1))^2+(yr-x(2))^2) 0 0;...
    -(x(2)-yr)/((xr-x(1))^2+(yr-x(2))^2) (x(1)-xr)/((xr-x(1))^2+(yr-x(2))^2) 0 0];
N=500;
t=1:N;
%% define state covariance matrix (model accuracy)
Q=blkdiag(zeros(2),eye(2)*10e-5);
%Q=eye(4)*10e-4;
%define sensor covariance matrix (sensor quality)
rr1=0.01+0.05*(square(pi/100*t)+1)/2;
rr2=0.005+0.02*(square(pi/100*t)+1)/2;
R0=diag([0.01 0.005]);
figure
plot(rr1)
hold on
plot(rr2)
%% memory allocation
x=zeros(4,N);
hx=zeros(4,N);
u=zeros(2,N);
y=zeros(2,N);
z=zeros(2,N);
hatR=zeros(2,2,N);
P=zeros(4,4,N);
mse=zeros(1,N);
trP=zeros(1,N);
%%
%define initial state covariance matrix
P(:,:,1)=blkdiag(1e4,1e4,1,1);

%% simulate system
%define initial conditions
x(:,1)=[1;0;0;0]; %true state
hx(:,1)=zeros(4,1); %estimated state
mse(1)=norm(x(:,1)-hx(:,1))^2; %initial mse
trP(1)=trace(P(:,:,1)); %initial covariance trace
tic
Nw=25;
for k=1:N
    %generate input signal
    u(:,k) = [1 0]';
    
    %generate true system signal
    x(:,k+1) = A*x(:,k)+B*u(:,k)+sqrt(Q)*randn(4,1);
    R=diag([rr1(k),rr2(k)]);
    y(:,k+1) = h(x(:,k+1))+sqrt(R)*randn(2,1);
    
    %% Filtering (EKF)
    %Step 1 - Prediction
    [hx(:,k+1),P(:,:,k+1)] = EKF_prediction(f,hx(:,k),u(:,k),P(:,:,k),Jf,Q);
    %% Adapt
    z(:,k+1)=y(:,k+1)-h(hx(:,k+1)); %innovation
    if k>Nw
        S=1/Nw*z(:,k+2-Nw:k+1)*z(:,k+2-Nw:k+1)';
        H=Jh(hx(:,k+1));
        hatR(:,:,k+1)=S-H*P(:,:,k+1)*H';
        if min(eig(hatR(:,:,k+1)))<=0
            hatR(:,:,k+1)=R0;
            warning('hatR inconsistent!')
        end
    else
        hatR(:,:,k+1)=R0;
    end
      
    %Step 2 - Update    
    [hx(:,k+1),P(:,:,k+1)] = EKF_update(h,hx(:,k+1),P(:,:,k+1),Jh,hatR(:,:,k+1),y(:,k+1));
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
legend('Real','EKF')
xlabel('Time (k)')
ylabel('x_k')
grid on
subplot(1,2,2)
plot(x(2,:),'g','linewidth',2)
hold on
plot(hx(2,:),'r','linewidth',2)
title('Distance in Y')
legend('Real','EKF')
xlabel('Time (k)')
ylabel('y_k')
grid on
figure
subplot(1,2,1)
plot(x(3,:),'g','linewidth',2)
hold on
plot(hx(3,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in X')
legend('Real','EKF')
xlabel('Time (k)')
ylabel('v^x_k')
grid on
subplot(1,2,2)
plot(x(4,:),'g','linewidth',2)
hold on
plot(hx(4,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in Y')
legend('Real','EKF')
xlabel('Time (k)')
ylabel('v_k^y')
grid on
figure
plot(x(1,:),x(2,:),'g','linewidth',2)
hold on
plot(hx(1,:),hx(2,:),'r-','linewidth',2)

grid on
legend('Real','EKF')
title('Trajectory')
xlabel('x_k')
ylabel('y_k')

%%
figure
plot(mse)
% Calculate the mean squared error over time
meanMSE = mean(mse);
disp(['Mean Squared Error: ', num2str(meanMSE)]);
%%
figure
subplot(2,1,1)
plot(rr1,displayName='true',linewidth=2)
hold on
plot(squeeze(hatR(1,1,:)),'--',DisplayName='hatR(1,1)',linewidth=2)

legend
subplot(2,1,2)
plot(rr2,displayName='true',linewidth=2)
hold on
plot(squeeze(hatR(2,2,:)),'--',DisplayName='hatR(2,2)',linewidth=2)

legend