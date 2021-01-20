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

[X,Y] = meshgrid(xRange,yRange); % az �sszes xRange �s yRange k�zti �rt�kre elv�gzi a Descartes szorz�st

dxdt = Y;
dydt = -(c/m)*Y - (k/m)*X;

figure
plot(t,x)
grid on

figure
plot(x(:,1),x(:,2)) %f�zis�llapot, ahol x tengelyen az els? �llapotv�ltoz�, y tengelyen a m�sodik
title('f�zis�llaopotok')
xlabel("elso �llapot v�ltoz�s")
ylabel("m�sodik �llapot v�ltoz�sa")
grid on

figure
quiver(X,Y,dxdt,dydt) %minden egyes pontp�rhoz hozz�rendel vektort
title('vektor rendel�se minden egyes pontp�rhoz')
grid on
hold on
plot(x(:,1),x(:,2))
hold off