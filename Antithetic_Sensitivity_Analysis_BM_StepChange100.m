clear all
close all
clc
LB = .1;
UB = 10;

par.mu = 1e-7;
par.alpha_B = .1; 

par.alpha_A = .1;  
par.alpha_C = .1;

par.gamma_AC = 5e5;

par.delta_A = 1e-4;
par.delta_B = 1e-4;
par.delta_C = 1e-4;

   %%%%%%% Global SA Metric Calculation %%%%%% 

N = 1000;
Point = zeros(8,N);
for ij = 1:8
    Point(ij,:) = randperm(N);
end

Interval = linspace(LB,UB,N+1);
Sample = zeros(N,8);

for i = 1:N
    for j = 1:8
        Sample(i,j) = rand*(Interval(Point(j,i)+1) - Interval(Point(j,i))) + Interval(Point(j,i));
    end
end

t_rise = zeros(N,1);
Max = zeros(N,1);
RelativeE_BMu = zeros(N,1);
Max_p_BMu = zeros(N,1);
S_time_BMu = zeros(N,1);
ss_diff = zeros(N,1);


tspan = 0:0.1:90000; 
t_cutoff = (length(tspan)-1)/(1.5);
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
    
    Reference(niter) = Vpar.mu/Vpar.alpha_B;
    
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    x0 = [0 0 0];
    [t,x] = ode23s(@Base_Model,tspan,x0,options,Vpar);
    L = length(tspan);

    n_a = x(end,1);
    n_b = x(end,2) + x(end,2);
    n_c = x(end,3);

    x02 = [n_a n_b n_c]; 

    [t2,x2] = ode23s(@Base_Model,tspan,x02,options,Vpar);

    All_x = x2;
    All_t = t2./60; %%% convert to mins
    
    Gene_A(niter,1) = All_x(t_cutoff,1);
    Gene_B(niter,1) = All_x(t_cutoff,2);
    Gene_C(niter,1) = All_x(t_cutoff,3);

    Gene_A(niter,2) = All_x(end,1);
    Gene_B(niter,2) = All_x(end,2);
    Gene_C(niter,2) = All_x(end,3);
    
    
    label_list = ["A mRNA"; "B mRNA"; "C mRNA"]; 
    
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
    Max_p_BMu(niter) = (max_Bact - Bact(end))/Bact(end); %/par.Px_tot;
    
    if locmax_Bact <= t_cutoff 
  
%             %%% Rise time
%             t90 = find(Bact(1:locmax_Bact) <= 0.9*max_Bact);
%             t10 = find(Bact(1:locmax_Bact) >= 0.1*max_Bact);
%             
%             rise = (t90(end) - t10(1))*t_increment; 

            %%% CHECK STEADY STATE
            t_end1 = find(Bact(locmax_Bact:end) > (1.05*Bact(end)));
            t_end2 = find(Bact(locmax_Bact:end) < (0.95*Bact(end)));
            
            if ~isempty(t_end1) && ~isempty(t_end2)
                S_time_BMu(niter) = max(t_end1(end),t_end2(end))*delta_t(1); 
                ss = max(t_end1(end),t_end2(end));
            elseif isempty(t_end1) && isempty(t_end2) 
                S_time_BMu(niter) = locmax_Bact*delta_t(1); 
                ss = locmax_Bact;
            elseif isempty(t_end1)
                S_time_BMu(niter) = t_end2(end)*delta_t(1);
                ss = t_end2(end);
            elseif isempty(t_end2) 
                S_time_BMu(niter) = t_end1(end)*delta_t(1);
                ss = t_end1(end);
            end
            de_avg = abs(Bact(ss) - Bact(end));
            
            %%% check if the simulation reaches steady state at 16.5 hours
            
            if de_avg <= ss_tolerance %%% mins value got from biomolecules paper
                
                success_trial = [success_trial;niter]; 

                %%% Steady State
                ss_diff(niter) = Bact(end); %%%% New definition

                %%% Steady state errorrrr
                RelativeE_BMu(niter) = (Bact(end) - Reference(niter))/Reference(niter);
            
            else
                
                fail_trial = [fail_trial;niter];     
               
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                RelativeE_BMu(niter) = nan;
                
            end
            
    else
                fail_trial = [fail_trial;niter];      
                S_time_BMu(niter) = nan;
                %%% Steady State
                ss_diff(niter) = nan; %%%% New definition
                RelativeE_BMu(niter) = nan;
               

    end
        
    
    clear All_x All_t
   
end

save('BM2_Global_SA_StepChange100.mat','Sample','t_rise','RelativeE_BMu',...
  'Max_p_BMu','S_time_BMu','success_trial', 'fail_trial','ss_diff','Gene_A','Gene_B','Gene_C','-v7.3')
