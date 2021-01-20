clc
clear
close all

% %Mass spring damper
m = 0.5; %kg
c = 0.1;
k = 0.3;

A = 0;
omega = 0.5;

tGrid = 0:0.01:7; %s

ut = tGrid;
u = A*sin(omega*ut); %Forcing term

x0 = [0.01,0]; %Initial condition

[t,x] = ode45(@(t,x)massSpringDamper(t,x,ut,u,m,k,c),tGrid,x0);


xRange = -1:0.1:1;
yRange = -1:0.1:1;

[X,Y] = meshgrid(xRange,yRange); % az összes xRange és yRange közti értékre elvégzi a Descartes szorzást

dxdt = Y;
dydt = -(c/m)*Y - (k/m)*X;

figure
plot(t,x)
grid on

figure
plot(x(:,1),x(:,2)) %fázisállapot, ahol x tengelyen az els? állapotváltozó, y tengelyen a második
title('fázisállaopotok')
xlabel("elso állapot változás")
ylabel("második állapot változása")
grid on

figure
quiver(X,Y,dxdt,dydt) %minden egyes pontpárhoz hozzárendel vektort
title('vektor rendelése minden egyes pontpárhoz')
grid on
hold on
plot(x(:,1),x(:,2))
hold off