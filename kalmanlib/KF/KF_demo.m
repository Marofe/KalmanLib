close all
clear all
clc

%% System Model
a=0.5;
b=0.5;
A=[0.9 0.1+0.06*a;0.01+0.06*b 0.9];
B1=[1 0 0;0 1 0];
C2=[1 0];
C1=[1 1];
D11=[0 0 0];
D21=[0 0 1.414];
H=C2;
%% Filter covariances
[m,n]=size(H);
sqrtQ=B1;
sqrtR=D21;
Q = sqrtQ*sqrtQ';
R = sqrtR*sqrtR';
%%
N=100; %time-horizon
M=1000; %Monte Carlo Experiments
err=zeros(M,N);
trP=zeros(N,1);
for j=1:M
    x=zeros(N,n);
    hx=zeros(N,n);
    y=zeros(N,m);
    x(1,:)=10*ones(1,n);
    hx(1,:)=zeros(1,n);
    P=eye(n)*100;
    trP(1)=trace(P);
    err(j,1)=norm(x(1,:)-hx(1,:))^2;
    for k=1:N-1
        %% real system
        x(k+1,:)=(A*x(k,:)'+sqrtQ*randn(3,1))';
        y(k+1,:)=(H*x(k+1,:)'+sqrtR*randn(3,1))';
        [hx(k+1,:),P]=KF(A,H,Q,R,hx(k,:),P,y(k+1,:));
        if j==1 %evaluate only on the first Monte Carlo Experiment
            trP(k+1)=trace(P);
        end
        %% erro
        err(j,k+1)=norm(x(k+1,:)-hx(k+1,:))^2;
    end
    if ~mod(j,round(M/100))
        clc
    fprintf("Monte Carlo simulation... %.2f%%\n",j/M*100)
    end
end
%% show
close all
plot(10*log10(trP),'Color',[0.8500, 0.3250, 0.0980],'linewidth',3)
hold on
plot(10*log10(mean(err)),'-','Color',[0, 0.4470, 0.7410],'linewidth',2)
xlabel('Time (k)')
ylabel('MSE (dB)')
legend('tr(P)','MSE')
grid on