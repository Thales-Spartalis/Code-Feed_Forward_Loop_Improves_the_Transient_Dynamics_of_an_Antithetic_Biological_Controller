tspan = 0:.01:20000;
tspan2 = 0:.01:100000; %medium
tspan3 = 0:.01:200000;%long
options = odeset('RelTol',1e-10,'AbsTol', 1e-10, 'OutputFcn', @odeOutputFcn);
for ij = 1:72

    if ij == 1
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        x0 = [0 0 0];  
        [t3,y0]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);
        Reference3 = ones(1,length(t3))*par.mu/par.alpha_B;  

    elseif ij == 2
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        n_a = y0(end,1);
        n_b = y0(end,2) + 0.1*y0(end,2);
        n_c = y0(end,3);

        x0 = [n_a n_b n_c];  
        [t3,y0b]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);
        Reference3 = ones(1,length(t3))*par.mu/par.alpha_B; 

    elseif ij == 3
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        x0 = [0 0 0];  
        [t3,y1]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);

    elseif ij == 4
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        n_a = y1(end,1);
        n_b = y1(end,2) + 0.1*y1(end,2);
        n_c = y1(end,3);

        x0 = [n_a n_b n_c];  
        [t3,y1b]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);
    
    elseif ij == 5
       par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        x0 = [0 0 0];  
        [t3,y2]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);

    elseif ij == 6
        par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        n_a = y2(end,1);
        n_b = y2(end,2) + 0.1*y2(end,2);
        n_c = y2(end,3);

        x0 = [n_a n_b n_c]; 
        [t3,y2b]=ode23s(@(t,x)Base_Model(t,x,par),tspan3,x0,options);

    elseif ij == 7
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        x0 = [0 0 0];  
        [t2,y3]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);
    elseif ij == 8
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        n_a = y3(end,1);
        n_b = y3(end,2) + 0.1*y3(end,2);
        n_c = y3(end,3);

        x0 = [n_a n_b n_c];  
        [t2,y3b]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);

    elseif ij == 9
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        x0 = [0 0 0];  
        [t2,y4]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);
        Reference2 = ones(1,length(t2))*par.mu/par.alpha_B;  

    elseif ij == 10
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        n_a = y4(end,1);
        n_b = y4(end,2) + 0.1*y4(end,2);
        n_c = y4(end,3);

        x0 = [n_a n_b n_c];  
        [t2,y4b]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);
        Reference2 = ones(1,length(t2))*par.mu/par.alpha_B; 

    elseif ij == 11
        par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        x0 = [0 0 0];  
        [t2,y5]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);

    elseif ij == 12
       par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        n_a = y5(end,1);
        n_b = y5(end,2) + 0.1*y5(end,2);
        n_c = y5(end,3);

        x0 = [n_a n_b n_c];  
        [t2,y5b]=ode23s(@(t,x)Base_Model(t,x,par),tspan2,x0,options);
    
    elseif ij == 13
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        x0 = [0 0 0];  
        [t,y6]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);

    elseif ij == 14
        par.mu = 1e-7; 
        
        par.alpha_A = .01; 
        par.alpha_B = .1;
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        n_a = y6(end,1);
        n_b = y6(end,2) + 0.1*y6(end,2);
        n_c = y6(end,3);

        x0 = [n_a n_b n_c]; 
        [t,y6b]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);

    elseif ij == 15
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        x0 = [0 0 0];  
        [t,y7]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);

    elseif ij == 16
        par.mu = 1e-7; 
        
        par.alpha_A = .1; 
        par.alpha_B = .1;
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        n_a = y7(end,1);
        n_b = y7(end,2) + 0.1*y7(end,2);
        n_c = y7(end,3);

        x0 = [n_a n_b n_c];  
        [t,y7b]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);

    elseif ij == 17
        par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        x0 = [0 0 0];  
        [t,y8]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;  

    elseif ij == 18
        par.mu = 1e-7; 
        
        par.alpha_A = 1; 
        par.alpha_B = .1;
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        n_a = y8(end,1);
        n_b = y8(end,2) + 0.1*y8(end,2);
        n_c = y8(end,3);

        x0 = [n_a n_b n_c];  
        [t,y8b]=ode23s(@(t,x)Base_Model(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;
%%% FB %%%
     elseif ij == 19
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];  
        [t3,y9]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options); 
       
    elseif ij == 20
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;   

        n_a = y9(end,1);
        n_b = y9(end,2) + 0.1*y9(end,2);
        n_c = y9(end,3);

        x0 = [n_a n_b n_c];   
        [t3,y9b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options);

    elseif ij == 21
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];  
        [t3,y10]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options);
 
    elseif ij == 22
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;

        n_a = y10(end,1);
        n_b = y10(end,2) + 0.1*y10(end,2);
        n_c = y10(end,3);
       
         
        x0 = [n_a n_b n_c];   
        [t3,y10b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options);

     elseif ij == 23
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0]; 
        [t3,y11]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options);

    elseif ij == 24
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.beta_CB = 1e5;

        n_a = y11(end,1);
        n_b = y11(end,2) + 0.1*y11(end,2);
        n_c = y11(end,3);
               
        x0 = [n_a n_b n_c];   
        [t3,y11b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan3,x0,options);

    elseif ij == 25
        
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0]; 
        [t2,y12]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);
               
    elseif ij == 26
        
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;

        n_a = y12(end,1);
        n_b = y12(end,2) + 0.1*y12(end,2);
        n_c = y12(end,3);
                 
        x0 = [n_a n_b n_c];    
        [t2,y12b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);

    elseif ij == 27
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];
        [t2,y13]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);

    elseif ij == 28
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;   

        n_a = y13(end,1);
        n_b = y13(end,2) + 0.1*y13(end,2);
        n_c = y13(end,3);
                 
        x0 = [n_a n_b n_c];   
        [t2,y13b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);

    elseif ij == 29
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];
        [t2,y14]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);
 
    elseif ij == 30
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.beta_CB = 1e5;

        n_a = y14(end,1);
        n_b = y14(end,2) + 0.1*y14(end,2);
        n_c = y14(end,3);
                 
        x0 = [n_a n_b n_c];   
        [t2,y14b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan2,x0,options);

     elseif ij == 31
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];
        [t,y15]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);

    elseif ij == 32
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;

        n_a = y15(end,1);
        n_b = y15(end,2) + 0.1*y15(end,2);
        n_c = y15(end,3);
                 
        x0 = [n_a n_b n_c];   
        [t,y15b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);

    elseif ij == 33
       
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];
        [t,y16]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);
               
    elseif ij == 34
       
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;

        n_a = y16(end,1);
        n_b = y16(end,2) + 0.1*y16(end,2);
        n_c = y16(end,3);
                 
        x0 = [n_a n_b n_c] ;  
        [t,y16b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);

    elseif ij == 35
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;   

        x0 = [0 0 0];
        [t,y17]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);

    elseif ij == 36
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.beta_CB = 1e5;   

        n_a = y17(end,1);
        n_b = y17(end,2) + 0.1*y17(end,2);
        n_c = y17(end,3);
                
        x0 = [n_a n_b n_c];   
        [t,y17b]=ode23s(@(t,x)FeedbackCB(t,x,par),tspan,x0,options);
%%% IFFL B %%%
    elseif ij == 37
       par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];         
        [t3,y18]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);               
    elseif ij == 38
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y18(end,1);
        n_b = y18(end,2) + 0.1*y18(end,2);
        n_x1 = y18(end,3);
        n_y1 = y18(end,4);
        n_z1 = y18(end,5);
        n_c = y18(end,6);
        n_x2 = y18(end,7);
        n_y2 = y18(end,8);
        n_z2 = y18(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];  
        [t3,y18b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);

    elseif ij == 39
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];
        [t3,y19]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);

    elseif ij == 40
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y19(end,1);
        n_b = y19(end,2) + 0.1*y19(end,2);
        n_x1 = y19(end,3);
        n_y1 = y19(end,4);
        n_z1 = y19(end,5);
        n_c = y19(end,6);
        n_x2 = y19(end,7);
        n_y2 = y19(end,8);
        n_z2 = y19(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];  
        [t3,y19b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);
    
    elseif ij == 41
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];  
        [t3,y20]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);

    elseif ij == 42
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y20(end,1);
        n_b = y20(end,2) + 0.1*y20(end,2);
        n_x1 = y20(end,3);
        n_y1 = y20(end,4);
        n_z1 = y20(end,5);
        n_c = y20(end,6);
        n_x2 = y20(end,7);
        n_y2 = y20(end,8);
        n_z2 = y20(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];   
        [t3,y20b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan3,x0,options);

    elseif ij == 43
       par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.alpha_X = .00008;
        par.beta_Y = 1e7;
        par.alpha_Z = .0001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];  
        [t2,y21]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);

    elseif ij == 44
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.alpha_X = .00008;
        par.beta_Y = 1e7;
        par.alpha_Z = .0001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y21(end,1);
        n_b = y21(end,2) + 0.1*y21(end,2);
        n_x1 = y21(end,3);
        n_y1 = y21(end,4);
        n_z1 = y21(end,5);
        n_c = y21(end,6);
        n_x2 = y21(end,7);
        n_y2 = y21(end,8);
        n_z2 = y21(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];   
        [t2,y21b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);

    elseif ij == 45
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005;  
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;


        x0 = [0 0 0 0 0 0 0 0 0];         
        [t2,y22]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);               
    elseif ij == 46
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y22(end,1);
        n_b = y22(end,2) + 0.1*y22(end,2);
        n_x1 = y22(end,3);
        n_y1 = y22(end,4);
        n_z1 = y22(end,5);
        n_c = y22(end,6);
        n_x2 = y22(end,7);
        n_y2 = y22(end,8);
        n_z2 = y22(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];  
        [t2,y22b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);

    elseif ij == 47
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00003;
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0]; 
        [t2,y23]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);

    elseif ij == 48
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;
        
       par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00003; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y23(end,1);
        n_b = y23(end,2) + 0.1*y23(end,2);
        n_x1 = y23(end,3);
        n_y1 = y23(end,4);
        n_z1 = y23(end,5);
        n_c = y23(end,6);
        n_x2 = y23(end,7);
        n_y2 = y23(end,8);
        n_z2 = y23(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];  
        [t2,y23b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan2,x0,options);
    
    elseif ij == 49
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];  
        [t,y24]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);

    elseif ij == 50
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y24(end,1);
        n_b = y24(end,2) + 0.1*y24(end,2);
        n_x1 = y24(end,3);
        n_y1 = y24(end,4);
        n_z1 = y24(end,5);
        n_c = y24(end,6);
        n_x2 = y24(end,7);
        n_y2 = y24(end,8);
        n_z2 = y24(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];   
        [t,y24b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);

    elseif ij == 51
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .001;
        par.beta_Y = 1e7;
        par.alpha_Z = .0003; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];  
        [t,y25]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);

    elseif ij == 52
      par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .001;
        par.beta_Y = 1e7;
        par.alpha_Z = .0003; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y25(end,1);
        n_b = y25(end,2) + 0.1*y25(end,2);
        n_x1 = y25(end,3);
        n_y1 = y25(end,4);
        n_z1 = y25(end,5);
        n_c = y25(end,6);
        n_x2 = y25(end,7);
        n_y2 = y25(end,8);
        n_z2 = y25(end,9);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];   
        [t,y25b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);

    elseif ij == 53
       par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0 0 0 0];         
        [t,y26]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);               
    elseif ij == 54
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;
        
        par.alpha_X = .00005;
        par.beta_Y = 1e7;
        par.alpha_Z = .00005; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y26(end,1);
        n_b = y26(end,2) + 0.1*y26(end,2);
        n_x = y26(end,3);
        n_y = y26(end,4);
        n_z = y26(end,5);
        n_c = y26(end,6);
  
        x0 = [n_a n_b n_x1 n_y1 n_z1 n_c n_x2 n_y2 n_z2];  
        [t,y26b]=ode23s(@(t,x)IFFL_2x(t,x,par),tspan,x0,options);
%%% IFFL A %%%
     elseif ij == 55
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t3,y27]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);
        Reference3 = ones(1,length(t3))*par.mu/par.alpha_B;

    elseif ij == 56
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

       par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y27(end,1);
        n_b = y27(end,2) + 0.1*y27(end,2);
        n_x = y27(end,3);
        n_y = y27(end,4);
        n_z = y27(end,5);
        n_c = y27(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];  
        [t3,y27b]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);
        Reference3 = ones(1,length(t3))*par.mu/par.alpha_B;

    elseif ij == 57
       par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

       par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t3,y28]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);
 
    elseif ij == 58
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y28(end,1);
        n_b = y28(end,2) + 0.1*y28(end,2);
        n_x = y28(end,3);
        n_y = y28(end,4);
        n_z = y28(end,5);
        n_c = y28(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t3,y28b]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);

     elseif ij == 59
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t3,y29]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);

    elseif ij == 60
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-5;
        par.delta_B = 1e-5;
        par.delta_C = 1e-5;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y29(end,1);
        n_b = y29(end,2) + 0.1*y29(end,2);
        n_x = y29(end,3);
        n_y = y29(end,4);
        n_z = y29(end,5);
        n_c = y29(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t3,y29b]=ode23s(@(t,x)IFFL1(t,x,par),tspan3,x0,options);

    elseif ij == 61
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t2,y30]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);
               
    elseif ij == 62
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y30(end,1);
        n_b = y30(end,2) + 0.1*y30(end,2);
        n_x = y30(end,3);
        n_y = y30(end,4);
        n_z = y30(end,5);
        n_c = y30(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t2,y30b]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);

    elseif ij == 63
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t2,y31]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);
        Reference2 = ones(1,length(t2))*par.mu/par.alpha_B;

    elseif ij == 64
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y31(end,1);
        n_b = y31(end,2) + 0.1*y31(end,2);
        n_x = y31(end,3);
        n_y = y31(end,4);
        n_z = y31(end,5);
        n_c = y31(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];  
        [t2,y31b]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);
        Reference2 = ones(1,length(t2))*par.mu/par.alpha_B;

    elseif ij == 65
       par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t2,y32]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);
 
    elseif ij == 66
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-4;
        par.delta_B = 1e-4;
        par.delta_C = 1e-4;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y32(end,1);
        n_b = y32(end,2) + 0.1*y32(end,2);
        n_x = y32(end,3);
        n_y = y32(end,4);
        n_z = y32(end,5);
        n_c = y32(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t2,y32b]=ode23s(@(t,x)IFFL1(t,x,par),tspan2,x0,options);

     elseif ij == 67
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t,y33]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);

    elseif ij == 68
        par.mu = 1e-7;
        
        par.alpha_A = .01; 
        par.alpha_B = .1; 
        par.alpha_C = .01;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y33(end,1);
        n_b = y33(end,2) + 0.1*y33(end,2);
        n_x = y33(end,3);
        n_y = y33(end,4);
        n_z = y33(end,5);
        n_c = y33(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t,y33b]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);

    elseif ij == 69
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t,y34]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;
               
    elseif ij == 70
        par.mu = 1e-7;
        
        par.alpha_A = .1; 
        par.alpha_B = .1; 
        par.alpha_C = .1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y34(end,1);
        n_b = y34(end,2) + 0.1*y34(end,2);
        n_x = y34(end,3);
        n_y = y34(end,4);
        n_z = y34(end,5);
        n_c = y34(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];   
        [t,y34b]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;

    elseif ij == 71
       par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        x0 = [0 0 0 0 0 0];  
        [t,y35]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;

    elseif ij == 72
        par.mu = 1e-7;
        
        par.alpha_A = 1; 
        par.alpha_B = .1; 
        par.alpha_C = 1;
        
        par.gamma_AC = 5e5;
        
        par.delta_A = 1e-3;
        par.delta_B = 1e-3;
        par.delta_C = 1e-3;

        par.alpha_X = 0.001;
        par.beta_Y = 1e5;
        par.alpha_Z = 0.001; 
        
        par.delta_X = 1e-4;
        par.delta_Y = 1e-4;
        par.delta_Z = 1e-4;

        n_a = y35(end,1);
        n_b = y35(end,2) + 0.1*y35(end,2);
        n_x = y35(end,3);
        n_y = y35(end,4);
        n_z = y35(end,5);
        n_c = y35(end,6);
         
        x0 = [n_a n_b n_x n_y n_z n_c];  
        [t,y35b]=ode23s(@(t,x)IFFL1(t,x,par),tspan,x0,options);
        Reference = ones(1,length(t))*par.mu/par.alpha_B;
    end

end
y0_new = [y0; y0b];
y1_new = [y1; y1b];
y2_new = [y2; y2b];
y3_new = [y3; y3b];
y4_new = [y4; y4b];
y5_new = [y5; y5b];
y6_new = [y6; y6b];
y7_new = [y7; y7b];
y8_new = [y8; y8b];
y9_new = [y9; y9b];
y10_new = [y10; y10b];
y11_new = [y11; y11b];
y12_new = [y12; y12b];
y13_new = [y13; y13b];
y14_new = [y14; y14b];
y15_new = [y15; y15b];
y16_new = [y16; y16b];
y17_new = [y17; y17b];
y18_new = [y18; y18b];
y19_new = [y19; y19b];
y20_new = [y20; y20b];
y21_new = [y21; y21b];
y22_new = [y22; y22b];
y23_new = [y23; y23b];
y24_new = [y24; y24b];
y25_new = [y25; y25b];
y26_new = [y26; y26b];
y27_new = [y27; y27b];
y28_new = [y28; y28b];
y29_new = [y29; y29b];
y30_new = [y30; y30b];
y31_new = [y31; y31b];
y32_new = [y32; y32b];
y33_new = [y33; y33b];
y34_new = [y34; y34b];
y35_new = [y35; y35b];

t_total =linspace (0, 4000000, 4000002);
t_total2 =linspace (0, 20000000, 20000002);
t_total3 =linspace (0, 40000000, 40000002);

Reference_final = [Reference Reference];
Reference_final2 = [Reference2 Reference2];
Reference_final3 = [Reference3 Reference3];

figure
fig = tiledlayout(3,3,'TileSpacing','compact');
txt = xlabel(fig,'Time (mins)');
txt2 = ylabel(fig,'Output and Reference');
txt.FontSize = 32;
txt2.FontSize = 32;
txt.FontName = 'Times New Roman';
txt2.FontName = 'Times New Roman';

nexttile
plot(t_total3./3600, Reference_final3,'--k','LineWidth',2)
hold on
ln = plot(t_total3./3600,  y0_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total3./3600,  y9_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total3./3600,  y18_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total3./3600,  y27_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
ylim([0 inf])
% xlim([0 3000])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based ','Interpreter','latex')
hold off

nexttile
plot(t_total3./3600, Reference_final3,'--k','LineWidth',2)
hold on
ln = plot(t_total3./3600,  y1_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total3./3600,  y10_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total3./3600,  y19_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total3./3600,  y28_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
ylim([0 inf])
% xlim([0 3000])
% xlim([0 1.5e6])
hold off


nexttile
plot(t_total3./3600, Reference_final3,'--k','LineWidth',2)
hold on
ln = plot(t_total3./3600,  y2_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total3./3600,  y11_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total3./3600,  y20_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total3./3600,  y29_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
ylim([0 inf])
% xlim([0 3000])
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
hold off

nexttile
plot( t_total2./3600, Reference_final2,'--k','LineWidth',2)
hold on
ln = plot(t_total2./3600,  y3_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total2./3600,  y12_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total2./3600,  y21_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total2./3600,  y30_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
% xlim([0 500])
hold off

nexttile
plot(t_total2./3600, Reference_final2,'--k','LineWidth',2)
hold on
ln = plot(t_total2./3600,  y4_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total2./3600,  y13_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total2./3600,  y22_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total2./3600,  y31_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
ylim([0 inf])
% xlim([0 500])
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based ','Interpreter','latex')
hold off

nexttile
plot(t_total2./3600, Reference_final2,'--k','LineWidth',2)
hold on
ln = plot(t_total2./3600,  y5_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total2./3600,  y14_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total2./3600,  y23_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total2./3600,  y32_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
ylim([0 inf])
% xlim([0 500])
% xlim([0 1.5e6])
hold off


nexttile
plot(t_total./3600, Reference_final,'--k','LineWidth',2)
hold on
ln = plot(t_total./3600,  y6_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total./3600,  y15_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total./3600,  y24_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total./3600,  y33_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
% ylim([0 12])
% xlim([0 300])
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
hold off

nexttile
plot(t_total./3600, Reference_final,'--k','LineWidth',2)
hold on
ln = plot(t_total./3600,  y7_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total./3600,  y16_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total./3600,  y25_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total./3600,  y34_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
% xlim([0 300])
hold off

nexttile
plot(t_total./3600, Reference_final,'--k','LineWidth',2)
hold on
ln = plot(t_total./3600,  y8_new(:,2),'LineWidth',3)
ln.Color = [0.32 0.3 0.38];
hold on
ln2 = plot(t_total./3600,  y17_new(:,2),'LineWidth',2)
ln2.Color = [0.4660 0.6740 0.1880];
hold on
ln4 = plot(t_total./3600,  y26_new(:,2),'LineWidth',2)
ln4.Color = [0.6350 0.0780 0.1840];
hold on
ln6 = plot(t_total./3600,  y35_new(:,2),'LineWidth',2)
ln6.Color = [0.9290 0.6940 0.1250];
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
% xlim([0 1.5e6])
legend ('Reference = $\frac{\mu}{\alpha_{B}}$', 'Base Circuit','FB-based','2IFFL-based', '1IFFL-based','Interpreter','latex')
% xlim([0 300])
hold off