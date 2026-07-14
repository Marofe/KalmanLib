close all
clear all
clc
addLibrary('kalmanlib')
rng(101)
%% KALMAN FILTER 

Ts = 0.1; %tempo de amostragem

%x=[px py vx vy]';

%matriz de estados para movimento 2D
A= [1 0 Ts 0;...
    0 1 0 Ts;...
    0 0 1 0;...
    0 0 0 1]
%matriz de entrada para MRUV
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
%% define matriz de covariance dos estados (precisao do modelo)
Q=eye(4)*10e-4;

%define matriz de covariance do sensor (qualidade do sensor)
R=0.01*eye(2);

N=100; %time window
Np=200; %number of particles

%% allocação memoria
x=zeros(4,N);
hx=zeros(4,N);
u=zeros(2,N);
y=zeros(2,N);
P=zeros(4,4,N);
mse=zeros(1,N);
Xi=zeros(4,Np);
%%
%define matriz de covariancia de estado inicial
P(:,:,1)=blkdiag(1e4,1e4,1,1);

%% simula sistema
%define condições iniciais
x(:,1)=[1;1;0;0]; %estado real
hx(:,1)=zeros(4,1); %estado estimado
mse(1)=norm(x(:,1)-hx(:,1))^2; %mse inicial
trP(1)=trace(P(:,:,1)); %trace covariancia inicial
tic
for k=1:N
    %gera sinal de entrada
    u(:,k) = [1 1]';
    
    %gera sinal do sistema real
    x(:,k+1) = A*x(:,k)+B*u(:,k)+sqrt(Q)*randn(4,1);
    y(:,k+1) = h(x(:,k+1))+sqrt(R)*randn(2,1);
    
    %% Filtragem (UKF)
    %Step 1 - Prediction
    [Xi,hx(:,k+1),P(:,:,k+1)] = PF_prediction(f,Xi,u(:,k),Q);
    %Step 2 - Update    
    [Xi,hx(:,k+1),P(:,:,k+1)] = PF_update(h,Xi,hx(:,k+1),R,y(:,k+1));
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
title('Distancia em X')
legend('Real','PF')
xlabel('Tempo (k)')
ylabel('x_k')
grid on
subplot(1,2,2)
plot(x(2,:),'g','linewidth',2)
hold on
plot(hx(2,:),'r','linewidth',2)
title('Distancia em Y')
legend('Real','PF')
xlabel('Tempo (k)')
ylabel('y_k')
grid on
figure
subplot(1,2,1)
plot(x(3,:),'g','linewidth',2)
hold on
plot(hx(3,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocidade em X')
legend('Real','PF')
xlabel('Tempo (k)')
ylabel('v^x_k')
grid on
subplot(1,2,2)
plot(x(4,:),'g','linewidth',2)
hold on
plot(hx(4,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocidade em Y')
legend('Real','PF')
xlabel('Tempo (k)')
ylabel('v_k^y')
grid on
figure
plot(x(1,:),x(2,:),'g','linewidth',2)
hold on
plot(hx(1,:),hx(2,:),'r-','linewidth',2)
grid on
legend('Real','PF')
title('Trajetoria')
xlabel('x_k')
ylabel('y_k')
%%
figure
plot(mse)
% Calculate the mean squared error over time
meanMSE = mean(mse);
disp(['Mean Squared Error: ', num2str(meanMSE)]);