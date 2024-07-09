function dx = CFFL1(t,x,par)

dx = zeros(6,1);

A = x(1); 
B = x(2);
X = x(3);
Y = x(4);
Z = x(5);
C = x(6);

dx(1) = par.mu + par.alpha_Z*Z - par.gamma_AC*A*C - par.delta_A*A; 
dx(2) = par.alpha_A*A - par.delta_B*B;
dx(3) = par.alpha_C*C - par.delta_X*X;
dx(4) = par.alpha_X*X - par.delta_Y*Y;
dx(5) = par.alpha_X*X + par.alpha_Y*Y - par.delta_Z*Z; 
dx(6) = par.alpha_B*B - par.gamma_AC*A*C - par.delta_C*C ; 


end
