close all
clear all
clc

%% real system
n=2;p=1;
% statistic properties
global tau
sigma=0.01*eye(n);
nu=0.1*eye(p);
tau=0.01;
% covariances
Q = sigma*sigma';
R = nu*nu';
%% steps
N=500;
%% Monte Carlo Simulation
for i=1:100
    clear x u hx hx2 y
    x=zeros(n,N);
    u=0;
    hx=x; %EKF
    hx2=x; %UKF
    x(:,1)=[pi/4;0]; %initial condition
    hx(:,1)=0.9*x(:,1);
    hx2(:,1)=hx(:,1);
    hx3(:,1)=hx(:,1);
    P=0.1*eye(n);
    P2=P;
    P3=P;
    % Time evaluation
    for k=1:N-1
        % real system
        %clc
        x(:,k+1)=f(x(:,k))+sigma*randn(n,1);
        y(:,k+1)=h(x(:,k+1))+nu*randn(p,1);
        %         x(:,k+1)
        %         y(:,k+1)
        %% EKF
        [hx(:,k+1),P,Ppred,K,hy]=extend_kf(@f,@h,@Ai,@Hi,Q,R,hx(:,k),P,y(:,k+1),0);
        %         P
        %         Ppred
        %         K
        %         hx(:,k+1)
        %         hy
        mse_ekf(i,k+1)=norm(hx(:,k+1)-x(:,k+1))^2;
        %% UKF
        [hx2(:,k+1),P2]=unscented_kf(@f,@h,Q,R,0.001,0,2,hx2(:,k),P2,y(:,k+1),0);
        mse_ukf(i,k+1)=norm(hx2(:,k+1)-x(:,k+1))^2;
        %                 clc
        %                 P
        %                 P2
        %                 hx(:,k+1)
        %                 hx2(:,k+1)
    end
    fprintf("%d %.2f %.2f\n",i,mean(mse_ukf(i,:)),mean(mse_ekf(i,:)))
end

%% plot
figure
plot(rad2deg(x(1,:)),'-k')
hold on
plot(rad2deg(hx(1,:)),'--g')
plot(rad2deg(hx2(1,:)),'--b')
hold off
figure
media1=mean(mse_ekf);
media2=mean(mse_ukf);
desvio1=std(mse_ekf)/sqrt(size(mse_ekf,1));
desvio2=std(mse_ukf)/sqrt(size(mse_ukf,1));
p2=plot(10*log10(media2),'Color',[0, 0.4470, 0.7410],'linewidth',2);
hold on
p1=plot(10*log10(media1),'Color',[0.8500, 0.3250, 0.0980],'linewidth',2);
F=1:500;
lx=10*log10(media2)-10*0.434*2*desvio2./media2;
lx(1)=0;
lx2=10*log10(media2)+10*0.434*2*desvio2./media2;
lx2(1)=0;
fill([F fliplr(F)],[lx fliplr(lx2)],[0, 0.4470, 0.7410], 'FaceAlpha', 0.2,'linestyle','none');
hold on
lx=10*log10(media1)-10*(0.434*2*desvio1./media1);
lx(1)=0;
lx2=10*log10(media1)+10*(0.434*2*desvio1./media1);
lx2(1)=0;
fill([F fliplr(F)],[lx fliplr(lx2)],[0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.2,'linestyle','none');

grid on
xlabel('Tempo (k)')
ylabel('MSE_{dB}')
axis([0 N min(10*log10(media2))-1 max(10*log10(media2))+1])
legend([p1,p2],'EKF','UKF')

%%
function x1=f(x,u)
global tau
%nonlinear dynamic system
x1(1,1)=x(1)+x(2)*tau;
x1(2,1)=x(2)-9.81*sin(x(1))*tau;
end

function y=h(x)
%nonlinear measurement from a radar
y=sin(x(1));
end

function M=Ai(x,u)
global tau
M=[1 tau;-9.81*cos(x(1))*tau 1];
end

function M=Hi(x)
M=[cos(x(1)) 0];
end