close all
clear all
clc

%%
n=2;
p=1;
a=0.5;
b=0.5;
A=[0.9 0.1+0.06*a;0.01+0.06*b 0.9]
iA=inv(A);
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
hxf(:,1)=zeros(n,1);
Pf=eye(n)*100;
trPf(1)=trace(Pf);
ef(j,1)=norm(x(:,1)-hxf(:,1))^2;
hxb(:,N)=zeros(n,1);
Pb=eye(n)*100;
trPb(N)=trace(Pb);
%% forward
for k=1:N-1
    %% real system
    x(:,k+1)=A*x(:,k)+sigma*randn(3,1);
    y(:,k+1)=H*x(:,k+1)+nu*randn(3,1);
    [hxf(:,k+1),Pf]=kalman_filter(A,H,Q,R,hxf(:,k),Pf,y(:,k+1));
    if j==1
    trPf(k+1)=trace(Pf);
    end
    %% erro
    ef(j,k+1)=norm(x(:,k+1)-hxf(:,k+1))^2;
end
%% backward
eb(j,N)=norm(x(:,N)-hxb(:,N))^2;
for k=N-1:-1:1
    [hxb(:,k),Pb]=kalman_filter(iA,H,Q,R,hxb(:,k+1),Pb,y(:,k));
    if j==1
    trPb(k)=trace(Pb);
    end
    %% erro
    eb(j,k)=norm(x(:,k)-hxb(:,k))^2;
end
end
%% show 
close all
plot(trPf,'Color',[0.8500, 0.3250, 0.0980],'linewidth',2)
hold on
plot(mean(ef),'--','Color',[0, 0.4470, 0.7410],'linewidth',2)
plot(trPb,'linewidth',2)
plot(mean(eb),'--','linewidth',2)

xlabel('Time (k)')
ylabel('Mean Square Error')
legend('tr(Pf)','MSE_f','tr(Pb)','MSE_b')
figure
plot(hxf(1,:),'linewidth',2)
hold on
plot(hxb(1,:),'--','linewidth',2)
legend('Forward','Backward')