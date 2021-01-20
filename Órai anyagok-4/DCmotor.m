function dxdt = DCmotor(t,x,ut,u,J,b,Ke,Kt,R,L)
u = interp1(ut, u, t);

dxdt = zeros(2,1);

% dxdt(1) = (Ke/J)*x(2)-(b/J)*x(1);
% dxdt(2) = (1/L)*u-(Ke/L)*x(1)-(R/L)*x(2);

A = [-b/J, Ke/J; -Ke/L, -R/L];
B = [0; 1/L];

dxdt = A*x + B*u;


end

