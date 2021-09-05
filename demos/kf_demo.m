close all
clear all
clc

%% real system
% Car tracking example 4.3 from book Särkkä, S. (2013). Bayesian Filtering and Smoothing.
% x=[x y dx dy] -> state vector
% y=[x;y] -> measurement vector
n=4; %number of states
p=2; %number of measurements
%% linear dynamic model
dt=0.1;
A=[1 0 dt 0;0 1 0 dt;0 0 1 0;0 0 0 1];
H=[1 0 0 0;0 1 0 0];
%% statistic properties
q1=1;
q2=1;
nu=5*eye(p);
% covariances
Q = [q1*dt^3 0 q1*dt^2 0;0 q2*dt^3/3 0 q2*dt^2/2;q1*dt^2/2 0 q1*dt 0;0 q2*dt^2/2 0 q2*dt];
sigma=sqrtm(Q);
R = nu*nu';
%% steps
N=500;
%% initial condition
x(:,1)=[0 0 1 1]';
hx(:,1)=x(:,1);
y(:,1)=H*x(:,1)+nu*randn(p,1);
P{1}=1e-3*eye(n);
%% Monte Carlo Simulation
for i=1:1
    % Forward
    for k=1:N-1
        %% real system
        x(:,k+1)=A*x(:,k)+sigma*randn(n,1);
        y(:,k+1)=H*x(:,k+1)+nu*randn(p,1);
        %% Kalman Filter
        [hx(:,k+1) P{k+1}]=kalman_filter(A,H,Q,R,hx(:,k),P{k},y(:,k+1));
    end
    % Backward
    hxs(:,N)=hx(:,N);
    Ps{N}=P{N};
    for k=N-1:-1:1
       %% Kalman Smooth
        [hxs(:,k) Ps{k}]=kalman_smooth(A,H,Q,R,hx(:,k),P{k},hxs(:,k+1),Ps{k+1}); 
    end
    %fprintf("%d %.2f %.2f\n",i,mean(mse_ukf(i,:)),mean(mse_ekf(i,:)))
end

%% plot
figure
plot(x(1,:),x(2,:),'linewidth',2)
hold on
plot(hx(1,:),hx(2,:),'linewidth',2)
plot(hxs(1,:),hxs(2,:),'linewidth',2)
legend('real','Filter','Smoothing')
grid on