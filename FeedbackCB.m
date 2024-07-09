function dx = FeedbackCB(t,x,par)

dx = zeros(3,1);

A = x(1); 
B = x(2);
C = x(3);

dx(1) = par.mu - par.gamma_AC*A*C - par.delta_A*A; 
dx(2) = par.alpha_A*A - par.beta_CB*C*B- par.delta_B*B;
dx(3) = par.alpha_B*B - par.gamma_AC*A*C - par.beta_CB*C*B - par.delta_C*C;

end