
%Inverted pendulum
c = 0.2;

tEnd = 100;

tGrid = 0:0.01:tEnd; %s

x0 = [deg2rad(1),0]; %Initial condition

[t,x] = ode45(@(t,x)invertedPendulum(t,x,c),[0 tEnd],x0);

xRange = -3*pi:0.1:3*pi;
yRange = -2:0.1:2;

[X,Y] = meshgrid(xRange,yRange); % 

dxdt = Y;
dydt = sin(X) - c*Y;

figure
plot(t,x)
legend('inga pozícioja')
grid on

figure
plot(x(:,1),x(:,2))
grid on

figure
quiver(X,Y,dxdt,dydt)
title("inga átborul az also egyensúlyi állapotba");
grid on
hold on
plot(x(:,1),x(:,2))
hold off


c = 0.2;

tEnd = 20;

tGrid = 0:0.01:tEnd; %s

x0 = [0.1,0]; %Initial condition

[t,x] = ode45(@(t,x)oscillator(t,x),[0 tEnd],x0);

xRange = -1.5:0.1:1.5;
yRange = -1.5:0.1:1.5;

[X,Y] = meshgrid(xRange,yRange);

dxdt = Y+ X.*(1-X.^2 - Y.^2);
dydt = -X + Y.*(1-X.^2-Y.^2);

figure
plot(t,x)
grid on

figure
plot(x(:,1),x(:,2))
grid on

figure
quiver(X,Y,dxdt,dydt)
grid on
hold on
plot(x(:,1),x(:,2))
hold off

sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,a] = ode45(f,[0 100],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
plot3(a(:,1),a(:,2),a(:,3))
grid on
