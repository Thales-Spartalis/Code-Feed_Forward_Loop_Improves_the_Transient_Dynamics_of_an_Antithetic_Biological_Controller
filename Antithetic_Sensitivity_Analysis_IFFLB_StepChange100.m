clear all
close all
clc
LB = .1;
UB = 10;

par.mu = 1e-7;

par.alpha_A = .1; 
par.alpha_B = .1; 
par.alpha_C = .1;

par.gamma_AC = 5e5;

par.delta_A = 1e-3;
par.delta_B = 1e-3;
par.delta_C = 1e-3;

par.alpha_X = 0.001;
par.beta_Y = 1e7;
par.alpha_Z = 0.0001; 

par.delta_X = 1e-4;
par.delta_Y = 1e-4;
par.delta_Z = 1e-4;

   %%%%%%% Global SA Metric Calculation %%%%%% 

N = 1000;
Point = zeros(14,N);
for ij = 1:14
    Point(ij,:) = randperm(N);
end

Interval = linspace(LB,UB,N+1);
Sample = zeros(N,14);

for i = 1:N
    for j = 1:14
        Sample(i,j) = rand*(Interval(Point(j,i)+1) - Interval(Point(j,i))) + Interval(Point(j,i));
    end
end
t_rise = zeros(N,1);
Max = zeros(N,1);
RelativeE_ifflB = zeros(N,1);
Max_p_ifflB = zeros(N,1);
S_time_ifflB = zeros(N,1);
ss_diff = zeros(N,1);


tspan = 0:0.1:90000; %%% 33 hours simulation
t_cutoff = (length(tspan)-1)/(1.5); %%% 22 hour cutoff
ss_tolerance = 1e-5;

success_trial = []; %%% store trials that meet the criteria.
fail_trial = []; %%% store trials did not meet the criteria.

Gene_A = zeros(N,2);
Gene_B = zeros(N,2);
Gene_C = zeros(N,2);

Reference = par.mu/par.alpha_B;

for niter = 1:N
    
    Vpar = par;
    Vpar.mu = par.mu*Sample(niter,1);

    Vpar.alpha_A = par.alpha_A*Sample(niter,2);
    Vpar.alpha_B = par.alpha_B*Sample(niter,3);
    Vpar.alpha_C = par.alpha_C*Sample(niter,4);
    Vpar.gamma_AC = par.gamma_AC*Sample(niter,5);

    Vpar.delta_A = par.delta_A*Sample(niter,6);
    Vpar.delta_B = par.delta_B*Sample(niter,7);
    Vpar.delta_C = par.delta_C*Sample(niter,8);
    Vpar.alpha_X = par.alpha_X*Sample(niter,9);
    Vpar.beta_Y = par.beta_Y*Sample(niter,10);
    Vpar.alpha_Z = par.alpha_Z*Sample(niter,11);

    Vpar.delta_X = par.delta_X*Sample(niter,12);
    Vpar.delta_Y = par.delta_Y*Sample(niter,13);
    Vpar.delta_Z = par.delta_Z*Sample(niter,14);

    Reference(niter) = Vpar.mu/Vpar.alpha_B;

    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    x0 = [0 0 0 0 0 0 0 0 0];
    [t,x] = ode23s(@IFFL_2x,tspan,x0,options,Vpar);
    L = length(tspan);

    n_a = x(end,1);
    n_b = x(end,2) + x(end,2);
    n_x1 = x(end,3);
    n_y1 = x(end,4);
    n_z1 = x(end,5);
    n_c = x(end,6);
    n_x2 = x(end,7);
    n_y2 = x(end,8);
    n_z2 = x(end,9);

    x02 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];

    [t2,x2] = ode23s(@IFFL_2x,tspan,x02,options,Vpar);

    All_x = x2;
    All_t = t2./60; %%% convert to mins
    
    Gene_A(niter,1) = All_x(t_cutoff,1);
    Gene_B(niter,1) = All_x(t_cutoff,2);
    Gene_C(niter,1) = All_x(t_cutoff,6);

    Gene_A(niter,2) = All_x(end,1);
    Gene_B(niter,2) = All_x(end,2);
    Gene_C(niter,2) = All_x(end,6);
    
    
    label_list = ["A mRNA";"B mRNA"; "x RNA"; "Y RNA";"Z RNA"; "C mRNA"]; 
    
    delta_t = diff(All_t);
    t_increment = delta_t(1); %%%% time interval in min.
    Bact = All_x(:,2); %%% B output
    
    %%% process the extremely small and possibly negative value, by setting
    %%% them to 0
   
    for ij = t_cutoff:length(Bact)
        if Bact(ij) <= 1e-15        
           Bact(ij) = 0;
        end    
    end
    
    [max_Bact, locmax_Bact] = max(Bact);
    [min_Bact, locmin_BmRNA] = min(Bact);
    
    %%% Max
    Max_p_ifflB(niter) = (max_Bact - Bact(end))/Bact(end); %/par.Px_tot;
    
    if locmax_Bact <= t_cutoff 
  
            %%% Rise time
%             t90 = find(Bact(1:locmax_Bact) <= 0.9*max_Bact);
%             t10 = find(Bact(1:locmax_Bact) >= 0.1*max_Bact);
%             
%             rise = (t90(end) - t10(1))*t_increment; 

            %%% CHECK STEADY STATE
            t_end1 = find(Bact(locmax_Bact:end) > (1.05*Bact(end)));
            t_end2 = find(Bact(locmax_Bact:end) < (0.95*Bact(end)));
            
            if ~isempty(t_end1) && ~isempty(t_end2)
                S_time_ifflB(niter) = max(t_end1(end),t_end2(end))*delta_t(1); 
                ss = max(t_end1(end),t_end2(end));
            elseif isempty(t_end1) && isempty(t_end2) 
                S_time_ifflB(niter) = locmax_Bact*delta_t(1); 
                ss = locmax_Bact;
            elseif isempty(t_end1)
                S_time_ifflB(niter) = t_end2(end)*delta_t(1);
                ss = t_end2(end);
            elseif isempty(t_end2) 
                S_time_ifflB(niter) = t_end1(end)*delta_t(1);
                ss = t_end1(end);
            end
            de_avg = abs(Bact(ss) - Bact(end));
            
            %%% check if the simulation reaches steady state at 16.5 hours
            
            if de_avg <= ss_tolerance %%% mins value got from biomolecules paper
                
                success_trial = [success_trial;niter]; 

                %%% Steady State
                ss_diff(niter) = Bact(end); %%%% New definition


                %%% Steady state errorrrr
                RelativeE_ifflB(niter) = (Bact(end) - Reference(niter))/Reference(niter);
            
            else
                
                fail_trial = [fail_trial;niter];     
               
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                RelativeE_ifflB(niter) = nan;
                
            end
            
    else
                fail_trial = [fail_trial;niter];      
                S_time_ifflB(niter) = nan;
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                RelativeE_ifflB(niter) = nan;
               
    end
        
    
    clear All_x All_t
   
end

save('ifflB_Global_SA_StepChange100.mat','Sample','RelativeE_ifflB',...
  'Max_p_ifflB','S_time_ifflB','success_trial', 'fail_trial','ss_diff','Gene_A','Gene_B','Gene_C','-v7.3')
