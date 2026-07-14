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
%% measurement model 
H=[eye(2) zeros(2)];

N=500;
t=1:N;
%% define state covariance matrix (model accuracy)
Q=blkdiag(zeros(2),eye(2)*10e-5);

%define sensor covariance matrix (sensor quality)
rr1=0.1+0.5*(square(pi/100*t)+1)/2;
rr2=0.05+0.2*(square(pi/50*t)+1)/2;
R0=diag([0.1 0.05]);
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
x(:,1)=[1;1;1;1]; %true state
hx(:,1)=zeros(4,1); %estimated state
mse(1)=norm(x(:,1)-hx(:,1))^2; %initial mse
trP(1)=trace(P(:,:,1)); %initial covariance trace
tic
Nw=11; 
for k=1:N
    %generate input signal (acceleration)
    u(:,k) = [1 1]'; 
    
    %generate true system signal
    x(:,k+1) = A*x(:,k)+B*u(:,k)+sqrt(Q)*randn(4,1);
    R=diag([rr1(k),rr2(k)]);
    y(:,k+1) = H*x(:,k+1)+sqrt(R)*randn(2,1);
    
    %% Filtering (KF)
    %Step 1 - Prediction
    [hx(:,k+1),P(:,:,k+1)] = KF_prediction(hx(:,k),u(:,k),P(:,:,k),A,B,Q);
    %% Adapt
    z(:,k+1)=y(:,k+1)-H*hx(:,k+1); %innovation
    if k>Nw
        S=1/Nw*z(:,k+2-Nw:k+1)*z(:,k+2-Nw:k+1)';
        hatR(:,:,k+1)=S-H*P(:,:,k+1)*H';
        if min(eig(hatR(:,:,k+1)))<=0
            hatR(:,:,k+1)=R0;
            warning('hatR inconsistent!')
        end
    else
        hatR(:,:,k+1)=R0;
    end
    %Step 2 - Update    
    [hx(:,k+1),P(:,:,k+1)] = KF_update(hx(:,k+1),P(:,:,k+1),H,hatR(:,:,k+1),y(:,k+1));

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
plot(y(1,:),'o-',Color=colors('gray'))
title('Distance in X')
legend('Real','KF')
xlabel('Time (k)')
ylabel('x_k')
grid on
subplot(1,2,2)
plot(x(2,:),'g','linewidth',2)
hold on
plot(hx(2,:),'r','linewidth',2)
plot(y(2,:),'o-',Color=colors('gray'))
title('Distance in Y')
legend('Real','KF')
xlabel('Time (k)')
ylabel('y_k')
grid on
figure
subplot(1,2,1)
plot(x(3,:),'g','linewidth',2)
hold on
plot(hx(3,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in X')
legend('Real','KF')
xlabel('Time (k)')
ylabel('v^x_k')
grid on
subplot(1,2,2)
plot(x(4,:),'g','linewidth',2)
hold on
plot(hx(4,:),'color',[0, 0.4470, 0.7410],'linewidth',2)
title('Velocity in Y')
legend('Real','KF')
xlabel('Time (k)')
ylabel('v_k^y')
grid on
%%
figure
plot(x(1,:),x(2,:),'g','linewidth',2)
hold on
plot(hx(1,:),hx(2,:),'ro-','linewidth',2)
plot(y(1,:),y(2,:),'ok',linewidth=1)
%% velocity vectors (quiver) - decimated for readability
step = 5; %plot a velocity arrow every 'step' samples
idx = 1:step:size(x,2);
scale = 0.3; %arrow scale factor
q1=quiver(x(1,idx), x(2,idx), scale*x(3,idx), scale*x(4,idx), 0, ...
    'Color',[0, 0.4470, 0.7410],'LineWidth',1.5,'MaxHeadSize',1.0);
q2=quiver(hx(1,idx), hx(2,idx), scale*hx(3,idx), scale*hx(4,idx), 0, ...
    'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5,'MaxHeadSize',1.0);
grid on
legend('Real trajectory','KF trajectory','Real velocity','KF velocity')
title('Trajectory with Velocity Vectors')
xlabel('x_k')
ylabel('y_k')
%%
figure
plot(mse)
hold on
plot(trP,'--')
% Calculate the mean squared error over time
meanMSE = mean(mse);
disp(['Mean Squared Error: ', num2str(meanMSE)]);
%%
figure
subplot(2,1,1)
plot(rr1,displayName='true')
hold on
plot(squeeze(hatR(1,1,:)),'--',DisplayName='hatR(1,1)')

legend
subplot(2,1,2)
plot(rr2,displayName='true')
hold on
plot(squeeze(hatR(2,2,:)),'--',DisplayName='hatR(2,2)')

legend