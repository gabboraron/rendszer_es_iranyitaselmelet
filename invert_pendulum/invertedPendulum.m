function dxdt = invertedPendulum(t,x,c)

dxdt = zeros(2,1);

dxdt(1) = x(2);
dxdt(2) = sin(x(1)) - c*x(2);

end

