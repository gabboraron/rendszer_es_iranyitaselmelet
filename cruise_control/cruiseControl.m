function dxdt = cruiseControl(t,x,Tm,beta,alphan,omegam,m,g,Cr,rho,Cd,A,ut,u,thetat,theta,vref)
% u = interp1(ut,u,t);
theta = interp1(thetat,theta,t);

kp = 0.2;
ki = 0.1;

u = kp*(vref-x(1)) + ki*x(2); % imputjel legyen arányos a követési hibával

Torque = Tm*(1-beta*(alphan*x(1)/omegam-1)^2);

dxdt = zeros(2,1);

dxdt(1) = (1/m)*alphan*u*Torque - g*Cr*sign(x(1)) - ...
    0.5*(1/m)*rho*Cd*A*x(1)^2 - g*sin(theta);
dxdt(2) = vref - x(1);
end

