clc
clear
close all

%Mass spring damper
% m = 0.5; %kg
% c = 0.1;
% k = 0.3;
% 
% A = 0;
% omega = 0.5;
% 
% tGrid = 0:0.01:50; %s
% 
% ut = tGrid;
% u = A*sin(omega*ut); %Forcing term
% 
% x0 = [1,0]; %Initial condition
% 
% [t,x] = ode45(@(t,x)massSpringDamper(t,x,ut,u,m,k,c),tGrid,x0);
% 
% figure
% plot(t,x)
% grid on

%Balancing System
M = 1;
m = 0.5; %kg
J = 0.3;
l = 0.5;
c = 0.3;
g = 9.81;
gamma = 0.2;

A = 0.1;
omega = 0.5;

tGrid = 0:0.1:50; %s

ut = tGrid;
%u = A*sin(omega*ut); %Forcing term
u = zeros(1,length(ut));
u(1) = 0.01;

x0 = [0,0.01,0,0]; %Initial condition

[t,x] = ode45(@(t,x)balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma),[0 10],x0);

figure
hold on
subplot(2,2,1);
plot(t,x(:,1:2)) %- inga kimozdulása id? függvényében
xlabel('ido')
ylabel('inga elmozdulása')
title('State Space Models - Ordinary Differential Equations')
subplot(1,2,2);
plot(t,x(:,1),t,rad2deg(x(:,2))) %- inga kimozdulása szögértékben
xlabel('ido')
ylabel('inga elmozdulása szögekben')
title('State Space Models - Ordinary Differential Equations')
hold off
grid on

%Vehicle steering
% v0 = 1;
% a = 0.3;
% b = 0.7;
% 
% A = 0.1;
% omega = 0.5;
% 
% tGrid = 0:0.1:50; %s
% 
% ut = tGrid;
% u = A*sin(omega*ut); %Forcing term
% %u = ones(1,length(ut))*0.1;
% % u(1) = 0.01;
% 
% x0 = [0,0,0]; %Initial condition
% 
% [t,x] = ode45(@(t,x)vehicleSteering(t,x,ut,u,v0,a,b),[0 50],x0);
% 
% figure
% plot(t,x)
% grid on
% 
% figure
% plot(x(:,1),x(:,2))
% grid on

A = [0,1,0,0,0;1,0,1,1,1;0,1,0,1,0;0,1,1,0,0;0,1,0,0,0];
D = diag([1,4,2,2,1]);
gamma = 0.01;

x0 = normrnd(5,0.1,5,1);

xk = x0;
xkk = x0;

t = 1:500;

xLog = zeros(length(t),length(x0));

for i=t
    xLog(i,:) = xk(:)';
    xkk = xk - gamma*(D-A)*xk;
    xk = xkk;   
end

figure
plot(t,xLog)
grid on