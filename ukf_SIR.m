%% SIR Model
T=60; %time-horizon (days)
N=1000; %Population
S=zeros(T,1);
I=zeros(T,1);
R=zeros(T,1);
Y=zeros(T,1);
y=zeros(T,1);
gamma=1/7;
beta=zeros(T,1);
%% Model Initial Condition
beta(1)=1; %Transmission rate
S(1)=N; %Suceptible
I(1)=2; %Infected
%%
f=@(x,u)[x(1)-x(4)*x(2)*x(1)/N;...
    (1-gamma)*x(2)+x(4)*x(2)*x(1)/N;...
    x(3)+gamma*x(2);...
    x(4)]; % Note: unscented_kf expects column vector if hx is column
h=@(x)x(4)*x(2)*x(1)/N;

% UT parameters
alpha=1e-3;
kappa=0;
beta_ut=2; % renamed to avoid conflict with SIR beta

hx=zeros(T,4);
sigmaR=2;
Q=zeros(4);
Q(4,4)=0.1;

x_real=zeros(T,4);
hx_est=zeros(T,4);
y_meas=zeros(T,1);

x_real(1,:)=[N 2 0 1];
hx_est(1,:)=[N 1 0 1];
P=eye(4);

for k=1:T-1
    %% Real system
    beta(k+1)=beta(k)+0.1*randn;
    if beta(k+1)<=0
        beta(k+1)=0.1*rand;
    end
    S(k+1)=S(k)-beta(k)*I(k)*S(k)/N;
    I(k+1)=(1-gamma)*I(k)+beta(k)*I(k)*S(k)/N;
    R(k+1)=R(k)+gamma*I(k);
    Y(k+1)=beta(k+1)*I(k+1)*S(k+1)/N;
    x_real(k+1,:)=[S(k+1) I(k+1) R(k+1) beta(k+1)];
    y_meas(k+1)=Y(k+1)+sqrt(sigmaR)*randn;
    
    %% UKF using unscented_kf.m
    % function [hx,P] = unscented_kf(f,h,Q,R,alpha,kappa,beta,x0,P0,y,u)
    [new_hx, P] = unscented_kf(f, h, Q, sigmaR, alpha, kappa, beta_ut, hx_est(k,:)', P, y_meas(k+1), []);
    hx_est(k+1,:) = new_hx';
end

%% Plotting
f1=figure;
p1=plot(S,'linewidth',2);
hold on
p2=plot(I,'linewidth',2);
p3=plot(R,'linewidth',2);
plot(hx_est(:,1:3),'k--','linewidth',2)
xlabel('Days')
ylabel('Individuals')
legend([p1 p2 p3], 'S (Susceptible)','I (Infected)','R (Recovered)')
grid on
%saveas(f1, 'unscented_kf_sir.png');

f2 = figure;
plot(x_real(:,4), 'b', 'linewidth', 2)
hold on
plot(hx_est(:,4), 'r--', 'linewidth', 2)
legend('\beta real', '\beta estimate')
title('Parameter Estimation with unscented\_kf.m')
grid on
%saveas(f2, 'unscented_kf_beta.png');
