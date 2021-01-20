function dxdt = massSpringDamper(t,x,ut,u,m,k,c)
u = interp1(ut,u,t);

dxdt = zeros(2,1);

dxdt(1) = x(2);
dxdt(2) = -(c/m)*x(2) - (k/m)*x(1) + u/m;

end

