%Cruise control
Tm = 190;
beta = 0.4;
alphan = 16;
omegam = 420;
m = 1302;
g = 9.81;
Cr = 0.01;
rho = 1.3;
Cd = 0.32;
A = 2.4;

v0 = 35;
vref = 40;
angle = 4;

tEnd = 100;

tGrid = 0:0.1:tEnd; %s

ut = tGrid;
u = ones(1,length(ut))*0.5;

thetat = tGrid;
theta = zeros(1,length(thetat));

iStart = find(tGrid == 20);  %bukkanó megadása
iEnd = find(tGrid == 30);  %

theta(iStart:end) = ones(1,length(theta(iStart:end)))*deg2rad(angle); %
theta(iStart:iEnd) = interp1([thetat(iStart),thetat(iEnd)],[0,deg2rad(angle)],thetat(iStart:iEnd)); % lineáris interpoláció, hogy a bukkanó ne éles szögben váltson
                                                                                                    % ismerjük a két pont koordinátáit és az x koordinátákat a kett? között

theta = zeros(1,length(thetat));


x0 = [v0,0]; %Initial condition

[t,x] = ode15s(@(t,x)cruiseControl(t,x,Tm,beta,alphan,omegam,m,g,Cr,rho,Cd,A,ut,u,thetat,theta,vref),[0 tEnd],x0); %integrálja a korábbi hibákat, hogy összeadja ?ket, hogy emlékezzen rájuk.

figure
plot(t,x)
grid on

figure
plot(thetat,theta)
grid on