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
%u(1) = 0.01;

x0 = [1,deg2rad(10),0,0]; %Initial condition, (kocsi,inga,kocsi_kezdo_sebessége,inga_kezdo_sebessege)
[t,x] = ode45(@(t,x)balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma),[0 10],x0);

figure
hold on
subplot(1,2,1);
plot(t,x(:,1),t,rad2deg(x(:,2))) %- inga kimozdulása szögértékben
legend('kiskocsi','inga');
xlabel('ido')
ylabel('elmozdulás szögben')
title('inga és kiskocsi mozása zárt rendszerben ha az ingának és kiskocsinak is van elmozdulása')

subplot(1,2,2);

x0 = [0,0,1,0]; %Initial condition, (kocsi,inga,kocsi_kezdo_sebessége,inga_kezdo_sebessege)
[t,x] = ode45(@(t,x)balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma),[0 10],x0);

plot(t,x(:,1),t,rad2deg(x(:,2))) %- inga kimozdulása szögértékben
legend('kiskocsi','inga');
xlabel('ido')
ylabel('elmozdulás szögben')
title('inga és kiskocsi mozása zárt rendszerben ha csak a kocsinak van kezdosebessége')
hold off
grid on
