function dxdt = balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma)
u = interp1(ut,u,t);

Mt = M + m;
Jt = J + m*l^2;

dxdt = zeros(4,1);

dxdt(1) = x(3);
dxdt(2) = x(4);
dxdt(3) = (-m*l*sin(x(2))*x(4) + m*g*((m*l^2)/Jt)*sin(x(2))*cos(x(2))...
    -c*x(3)+u)/(Mt - m*((m*l^2)/Jt)*cos(x(2))^2);
dxdt(4) = (-m*l^2*sin(x(2))*cos(x(2))*x(4)^2 + Mt*g*l*sin(x(2)) +...
    c*l*cos(x(2))*x(3) + gamma*x(4) + l*cos(x(2))*u)/(Jt*(Mt/m)-m*(l*cos(x(2)))^2);

end

