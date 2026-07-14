close all
clear all
clc

%% real system
n=3;p=2;m=2;
% statistic properties
global tau
sigma=0.01*eye(n);
nu=diag([0.5 0.5]);
tau=0.1;
% covariances
Q = sigma*sigma';
R = nu*nu';
%% steps
N=500;
%% Monte Carlo Simulation
for i=1:500
    clear x u hx hx2 y
    x=zeros(n,N);
    hx=x; %EKF
    hx2=x; %UKF
    x(:,1)=[50 50 0]'; %initial condition
    hx(:,1)=0.9*x(:,1);
    hx2(:,1)=hx(:,1);
    hx3(:,1)=hx(:,1);
    P=eye(n);
    P2=P;
    P3=P;
    % Time evaluation
    for k=1:N-1
        % real system
        %clc
        if k<N/2
            u(:,k)=[1;-0.2];
        else
            u(:,k)=[1;0.2];
        end
        x(:,k+1)=f(x(:,k),u(:,k))+sigma*randn(n,1);
        y(:,k+1)=h(x(:,k+1))+nu*randn(p,1);
        %% EKF
        [hx(:,k+1),P]=extend_kf(@f,@h,@Ai,@Hi,Q,R,hx(:,k),P,y(:,k+1),u(:,k));
        mse_ekf(i,k+1)=norm(hx(:,k+1)-x(:,k+1))^2;
        %% UKF
        [hx2(:,k+1),P2]=unscented_kf(@f,@h,Q,R,0.001,0,2,hx2(:,k),P2,y(:,k+1),u(:,k));
        mse_ukf(i,k+1)=norm(hx2(:,k+1)-x(:,k+1))^2;
    end
    fprintf("%d %.2f %.2f\n",i,mean(mse_ukf(i,:)),mean(mse_ekf(i,:)))
end

%% plot
figure
plot(x(1,:),x(2,:),'k--','linewidth',1)
hold on
plot(hx(1,:),hx(2,:),'-r','linewidth',1)
plot(hx2(1,:),hx2(2,:),'-g','linewidth',1)


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
axis([0 N min(10*log10(media2))-0.1*abs(min(10*log10(media2))) max(10*log10(media2))+0.1*abs(max(10*log10(media2)))])
legend([p1,p2],'EKF','UKF')

%% Model
function x1=f(x,u)
global tau
%nonlinear dynamic system
x1(1,1)=x(1)+u(1)*cos(x(3))*tau;
x1(2,1)=x(2)+u(1)*sin(x(3))*tau;
x1(3,1)=x(3)+u(2)*tau;
end

function M=Ai(x,u)
global tau
M=[1 0 -u(1)*sin(x(3))*tau;...
    0 1 u(1)*cos(x(3))*tau;...
    0 0 1];
end

function y=h(x)
%nonlinear measurement
y(1,1)=sqrt(x(1)^2+x(2)^2);
y(2,1)=atan(x(2)/x(1));
end

function M=Hi(x)
H11=x(1)/sqrt(x(1)^2+x(2)^2);
H12=x(2)/sqrt(x(1)^2+x(2)^2);
H21=-2*x(2)^2/(x(1)^3+x(2)^2*x(1));
H22=2*x(2)/(x(1)^2+x(2)^2);
M=[H11 H12 0;H21 H22 0];
end