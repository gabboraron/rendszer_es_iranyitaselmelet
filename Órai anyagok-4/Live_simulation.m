clc
clear
close all

%Stability test
tEnd = 3.5;

J = 0.01;
b = 0.1;
Ke = 0.01;
Kt = 0.01;
R = 1;
L = 0.5;

x0 = [0,0.5];

tGrid = 0:0.05:tEnd;
ut = tGrid;
u = ones(1,length(tGrid))*1;
%u = sin(tGrid)*0.5;

[t,x] = ode45(@(t,x)DCmotor(t,x,ut,u,J,b,Ke,Kt,R,L),[0,tEnd],x0);

A = [-b/J, Ke/J; -Ke/L, -R/L];
B = [0; 1/L];
C = [1,0];
D = 0;

sys = ss(A,B,C,D);
step(sys)
impulse(sys)

intSol = 0;
xSol(1,1:2) = x0;
for i=2:length(tGrid)
    for j=2:i
        intSol = intSol + (expm(A.*(tGrid(i)-tGrid(j)))*B*u(j) + expm(A.*(tGrid(i)-tGrid(j-1)))*B*u(j-1))*(0.05/2);
    end
    xSol(i,1:2) = expm(A.*tGrid(i))*x0' + intSol;
    intSol = 0;
end

figure()
plot(t,x(:,1),tGrid,xSol(:,1));

[T,D] = eig(A);

A
T*D*inv(T)

xRange = -1:0.1:1;
yRange = -2:0.1:2;

[X,Y] = meshgrid(xRange,yRange);

dxdt = 2*X-Y;
dydt = -X+2*Y;


figure()
plot(x(:,1),x(:,2))

figure()
quiver(X,Y,dxdt,dydt,2)
hold on
plot(x(:,1),x(:,2))
hold off

%Tanker stability
a1 = -0.6;
a2 = -0.3;
a3 = -5;
a4 = -2;
alpha = -2;

x1e = 0.075;

A = [a1 + 2*alpha*abs(x1e),a2;a3,a4]

eig(A)
