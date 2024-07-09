function dx = IFFL_2x(t,x,par)

dx = zeros(9,1);

A = x(1); 
B = x(2);
X1 = x(3);
Y1 = x(4);
Z1 = x(5);
C = x(6);
X2 = x(7);
Y2 = x(8);
Z2 = x(9);

dx(1) = par.mu - par.gamma_AC*A*C - par.delta_A*A; 
dx(2) = par.alpha_Z*Z1 - par.delta_B*B;
dx(3) = par.alpha_A*A - par.delta_X*X1;
dx(4) = par.alpha_X*X1 - par.beta_Y*Y1*Z1 - par.delta_Y*Y1;
dx(5) = par.alpha_X*X1 - par.beta_Y*Y1*Z1 - par.delta_Z*Z1; 
dx(6) = par.alpha_Z*Z2 - par.gamma_AC*A*C - par.delta_C*C ; 
dx(7) = par.alpha_B*B - par.delta_X*X2; 
dx(8) = par.alpha_X*X2 - par.beta_Y*Y2*Z2 - par.delta_Y*Y2;
dx(9) = par.alpha_X*X2 - par.beta_Y*Y2*Z2 - par.delta_Z*Z2;


end