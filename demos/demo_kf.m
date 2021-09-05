close all
clear all
clc

%%
n=2;
p=1;
a=0.5;
b=0.5;
A=[0.9 0.1+0.06*a;0.01+0.06*b 0.9]
B1=[1 0 0;0 1 0];
C2=[1 0];
C1=[1 1];
D11=[0 0 0];
D21=[0 0 1.414];
%%
sigma=B1;
nu=D21;
%%
Q = sigma*sigma';
R = nu*nu';
H=C2;
%%
%%
N=100;
for j=1:100
clear x hx 
x(:,1)=10*ones(n,1);
hx(:,1)=zeros(n,1);
P=eye(n)*100;
trP(1)=trace(P);
e(j,1)=norm(hx(:,1)-x(:,1))^2;
for k=1:N-1
    %% real system
    x(:,k+1)=A*x(:,k)+sigma*randn(3,1);
    y(:,k+1)=H*x(:,k+1)+nu*randn(3,1);
    [hx(:,k+1),P]=kalman_filter(A,H,Q,R,hx(:,k),P,y(:,k+1));
    if j==1
    trP(k+1)=trace(P);
    end
    %% erro
    e(j,k+1)=norm(x(:,k+1)-hx(:,k+1))^2;
end
end
%% show 
close all
plot(trP,'Color',[0.8500, 0.3250, 0.0980],'linewidth',3)
hold on
plot(mean(e),'-','Color',[0, 0.4470, 0.7410],'linewidth',2)
xlabel('Time (k)')
ylabel('Erro Quadrado Médio')
legend('tr(P)','MSE')