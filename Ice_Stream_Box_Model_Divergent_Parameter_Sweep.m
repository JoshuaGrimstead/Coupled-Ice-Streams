%% Compute 12D Poincare Sections and Switches (Grimstead et al., ,2026)
% This is the high resolution parameter sweep for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026.  
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run and the times and indices when the event 
% 'box 2 void ratio crossing 0.6 increasing' is triggered are recorded, 
% in addition to the 12D parameter space at that event. Events determining
% when the till enters a different state are also recorded. 
% 
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
% Base values and initial conditions recorded in p_base_values.mat
%
% Code written by Joshua Grimstead. 
%
%%
% set the number of cores for the HPC to use
parpool(48);

%%
% T_s2_vals from 9 to 25 with step 0.01
T_s2_vals = 9.0 : 0.01 : 25.0;
% L1_vals from 0 to 200 with step 0.4
L1_vals_raw = 0 : 0.4 : 200;
% adjusting for km -> m
L1_vals = L1_vals_raw .*1e3;

T_s2_len = numel(T_s2_vals);
L1_len = numel(L1_vals);

%Initialising empty cells

FailedPixelsMask = false(T_s2_len, L1_len);
n_points_mat = nan(T_s2_len, L1_len);
n_points_mat_e3 = nan(T_s2_len, L1_len);
PC_results_e2 = cell(T_s2_len, L1_len);
PC_results_e3 = cell(T_s2_len, L1_len);
peak_times_e2_mat = cell(T_s2_len, L1_len);
peak_times_e3_mat = cell(T_s2_len, L1_len);
crossing_times_e2_mat = cell(T_s2_len, L1_len);
crossing_times_e3_mat = cell(T_s2_len, L1_len);
ThresholdCrossingIndiciesMat = zeros(T_s2_len, L1_len,9);


%% Create a DataQueue for progress reporting
total_iters = T_s2_len * L1_len;

dq = parallel.pool.DataQueue;
t_start = tic;

afterEach(dq, @(iter_time) updateProgress(iter_time, total_iters, t_start));

dq_state = parallel.pool.DataQueue;
afterEach(dq_state, @(x) printIterationInfo4(x));

% load parameter file
load("p_base_values.mat")
p.L2 = 215e3;
p.L1f =   p.L1 ;
p.L2f =   p.L2 ;
p.L3f =   p.L3 ;
p.T_s1f = p.T_s1 ; 
p.T_s2f = p.T_s2 ;
p.T_s3f = p.T_s3 ; 
p.rate=0;


t_final = 20e4;


iter_start_total = tic;

parfor h = 1:T_s2_len
    
  
    p_loop = p;      % copy of conditions for the loop

    p_loop.T_s2 = T_s2_vals(h);
    p_loop.T_s2f = p_loop.T_s2 ;

    row_failed = false(1, L1_len);
    row_T_data = cell(1, L1_len);
    
    for i = 1:L1_len

        iter_start = tic;
        t_iter = tic;
        PC_cell = cell(30,1);
        p_loop.L1 = L1_vals(i);
        p_loop.L1f = p_loop.L1;
        

        global_iter = (h-1)*L1_len + i;

        send(dq_state, struct( ...
            'iter', global_iter, ...
            'T_s2', p_loop.T_s2, ...
            'L1', p_loop.L1 ...
        ));


        % 1. Spin up attractor run to equilibrium
        p_loop.ic = p.ic; 
        
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                 'Events', @(t,X) myEvents(t, X, p_loop));
        
        try
            % Attempt standard integration
            [~, Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p_loop), ...
                               [0, p_loop.year*t_final], p_loop.ic, options);
            
            % 2. Production run
            t_final_prod = t_final * 30; 
            p_loop.ic = Tinit(end,:); 
            % t: Time points
            % T: The full state trajectory
            % te: Time of events
            % ye: State at events
            % ie: Index of events
            tspan = linspace(0,p_loop.year*t_final,5000);

            [t, T, te, ye, ie] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p_loop), ...
                           tspan, p_loop.ic, options);
            
           
        catch ME
            row_failed(i) = true;
            iter_time = toc(t_iter);
            send(dq, iter_time);
            continue
            
        end
        
        % ie contains which event triggered (1 or 2), ye contains the state 
        % Filter results s.t not enough data doesn't cause problems
        e2_indices = (ie == 1);
        if sum(e2_indices) < 4
             iter_time = toc(t_iter);
            send(dq, iter_time);
            continue
        end

        e3_indices = (ie == 2);
         % te is a vector of times corresponding to the indices in ie
        e2_crossing_times = te(ie == 1); 
        e3_crossing_times = te(ie == 2); 

        dampening_indicies1 = (ie == 6);
        dampening_indicies2 = (ie == 7);
        dampening_indicies3 = (ie == 8);

        max_till_indicies1 = (ie == 3);
        max_till_indicies2 = (ie == 4); 
        max_till_indicies3 = (ie == 5);

        voidc_ratio_indicies1  = (ie == 9);
        voidc_ratio_indicies2  = (ie == 10);
        voidc_ratio_indicies3  = (ie == 11);
        
        
        [~, peak_locs_e2] = findpeaks(T(:,6));
        [~, peak_locs_e3] = findpeaks(T(:,10));
        
        e2_peak_times = t(peak_locs_e2);
        e3_peak_times = t(peak_locs_e3);

        

     

        temp_crossings_binary = zeros(1, 1, 9);

        if any(dampening_indicies1)
            temp_crossings_binary(1) = 1;
        end

        if any(dampening_indicies2)
            temp_crossings_binary(2) = 1;
        end

        if any(dampening_indicies3)
            temp_crossings_binary(3) = 1;
        end
        
        if any(max_till_indicies1)
            temp_crossings_binary(4) = 1;
        end

        if any(max_till_indicies2)
            temp_crossings_binary(5) = 1;
        end

        if any(max_till_indicies3)
            temp_crossings_binary(6) = 1;
        end

        if any(voidc_ratio_indicies1)
            temp_crossings_binary(7) = 1;
        end

        if any(voidc_ratio_indicies2)
            temp_crossings_binary(8) = 1;
        end

        if any(voidc_ratio_indicies3)
            temp_crossings_binary(9) = 1;
        end
            
       ThresholdCrossingIndiciesMat(h,i,:) = temp_crossings_binary; 

        V1 = (ye(e2_indices, 1)*p_loop.L1*p_loop.W1)/10^9;
        V2 = (ye(e2_indices, 5)*p_loop.L2*p_loop.W2)/10^9; 
        V3 = (ye(e2_indices, 9)*p_loop.L3*p_loop.W3)/10^9;
        

        V1_e3 = (ye(e3_indices, 1)*p_loop.L1*p_loop.W1)/10^9;
        V2_e3 = (ye(e3_indices, 5)*p_loop.L2*p_loop.W2)/10^9; 
        V3_e3 = (ye(e3_indices, 9)*p_loop.L3*p_loop.W3)/10^9;
        
        % Force column vectors to be 2-column matrices, always
        if ~isempty(e2_indices) && any(e2_indices)
            % Ensure these are column vectors
            col1 = (V1 + V2 + V3); 
            col2 = ye(e2_indices, 10);
            
            % Force them into a 2-column matrix [N x 2]
            PC_2D = [col1(:), col2(:)]; 
            PC_all = ye(e2_indices, :);
        else
            PC_2D = []; % Explicitly handle empty cases
            PC_all = [];
        end
        
        % Do the same for PC_2D_e3
        if ~isempty(e3_indices) && any(e3_indices)
            col1_e3 = (V1_e3 + V2_e3 + V3_e3);
            col2_e3 = ye(e3_indices, 6);
            PC_2D_e3 = [col1_e3(:), col2_e3(:)];
        else
            PC_2D_e3 = [];
        end
        
       
      
        % capturing points depending on the till height recorded
        n_points = size(PC_2D, 1);
        n_points_e3 = size(PC_2D_e3, 1);
       

       
  
        %save all parameters of interest
        % peak indecies e2
        % peak indecies e3
        %PC_2D
        %PC_2D_e3
        %n_points
        %n_points_e3
        % peak times e_2
        % peak times e_3
        %ThermalDampeningMat
        
        n_points_mat(h,i) = n_points;
        n_points_mat_e3(h,i) = n_points_e3;
        PC_results_e2_all{h,i} = PC_all;
        PC_results_e2{h,i} = PC_2D;
        PC_results_e3{h,i} = PC_2D_e3;
        peak_times_e2_mat{h,i} = e2_peak_times;
        peak_times_e3_mat{h,i} = e3_peak_times;
        crossing_times_e2_mat{h,i} = e2_crossing_times;
        crossing_times_e3_mat{h,i} = e3_crossing_times;

        iter_time = toc(t_iter);
        send(dq, iter_time);

    end
    FailedPixelsMask(h, :) = row_failed;
    
end
%% Save Data
% Saving number of points, 12D parameter state at void ratio crossings and
% associated times, indices where till switches are crossed, and where not
% enough data was available as NAN. 
save('raw_simulation2_output3_HPC.mat','-v7.3', 'n_points_mat', 'n_points_mat_e3', ...
    'PC_results_e2_all','PC_results_e2', 'PC_results_e3', 'peak_times_e2_mat', 'peak_times_e3_mat', ...
    "crossing_times_e2_mat", "crossing_times_e3_mat", ...
    'ThresholdCrossingIndiciesMat', 'FailedPixelsMask', 'T_s2_vals', 'L1_vals');


function [value, isterminal, direction] = myEvents(t, X, p_loop)
 
    thresh = 0.6;
    min_e = p_loop.e_c; % set as 0.3

    tol = 1e-6;  % choose arbitrary small tolerance to ensure detection
  
   
    
    value = [X(6) - thresh;       % Event 1: Cross e_2=0.6
             X(10) - thresh;      % Event 2: Cross e_3=0.6
           
             X(3) - p_loop.htmax1 + tol; % Event 3: Cross till_1 ~ 0
             X(7) - p_loop.htmax2 + tol; % Event 4: Cross till_2 ~ 0
             X(11) - p_loop.htmax3 + tol; % Event 5: Cross till_3 ~ 0
            
             X(3) - p_loop.h_t_min; % Event 6: Cross till_1 ~ max
             X(7) - p_loop.h_t_min; % Event 7: Cross till_2 ~ max
             X(11) - p_loop.h_t_min; % Event 8: Cross till_3 ~ max

             X(2) - min_e;       % Event 9: Cross e_1=0.3
             X(6) - min_e;      % Event 10: Cross e_2=0.3
             X(10) - min_e];       % Event 11: Cross e_3=0.3
   
     %choosing direction of event detection              
             
    isterminal = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
    direction  = [1; 1; 0; 0; 0; -1; -1; -1; -1; -1; -1]; 
end