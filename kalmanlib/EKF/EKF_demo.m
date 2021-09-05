close all
clear all
clc

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
R0=beta(1)/gamma;
%State x=[S I R beta]
f=@(x)[x(1)-x(4)*x(2)*x(1)/N;...
        (1-gamma)*x(2)+x(4)*x(2)*x(1)/N;...
        x(3)+gamma*x(2);...
        x(4)]';
h=@(x)x(4)*x(2)*x(1)/N;
%% Jacobian
Ax=@(x)[(1-x(4)*x(2)/N) -x(4)*x(1)/N 0 -x(2)*x(1)/N;...
                x(4)*x(2)/N (1-gamma+x(4)*x(1)/N) 0 x(2)*x(1)/N;...
                0 gamma 0 0;...
                0 0 0 1];
Hx=@(x)[x(4)*x(2)/N x(4)*x(1)/N 0 x(2)*x(1)/N]; 
n=4;
m=1;
%% Initial Condition for EKF
M=1000; %Monte Carlo Experiments
hx=zeros(T,4);
sigmaR=2;
Q=zeros(4);
Q(4,4)=0.1;
err=zeros(M,T);
trP=zeros(M,T);
for j=1:M
    x=zeros(T,n);
    hx=zeros(T,n);
    y=zeros(T,m);
   x(1,:)=[N 2 0 1];
   hx(1,:)=[N 1 0 1];
   P=eye(4);
   trP(j,1)=trace(P);
   err(j,1)=norm(x(1,:)-hx(1,:))^2;
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
        x(k+1,:)=[S(k+1) I(k+1) R(k+1) beta(k+1)];
        y(k+1)=Y(k+1)+sqrt(sigmaR)*randn;
    %% EKF
    [hx(k+1,:),P]=EKF(f,h,Ax,Hx,Q,sigmaR,hx(k,:),P,y(k+1,:));
    %% error
     trP(j,k+1)=trace(P);
     err(j,k+1)=norm(x(k+1,:)-hx(k+1,:))^2;
end
    if ~mod(j,round(M/100))
        clc
    fprintf("Monte Carlo simulation... %.2f%%\n",j/M*100)
    end
end
 %%
 close all
 figure
 p1=plot(S,'linewidth',2);
 hold on
 p2=plot(I,'linewidth',2);
 p3=plot(R,'linewidth',2);
 plot(hx(:,1:3),'k--','linewidth',2)
 xlabel('Days')
 ylabel('Individuals')
 legend([p1 p2 p3], 'S (Susceptible)','I (Infected)','R (Recovered)')
 grid on
 hold on
 figure
 plot(beta)
 plot(hx(:,4))
 line([1 T],[1 1],'Color','k')
 legend('$\beta$','$\hat{\beta}$','Interpreter','latex')
 grid on
 %% errors
plot(10*log10(mean(trP)),'Color',[0.8500, 0.3250, 0.0980],'linewidth',3)
hold on
plot(10*log10(mean(err)),'-','Color',[0, 0.4470, 0.7410],'linewidth',2)
xlabel('Time (k)')
ylabel('MSE (dB)')
legend('tr(P)','MSE')
grid on