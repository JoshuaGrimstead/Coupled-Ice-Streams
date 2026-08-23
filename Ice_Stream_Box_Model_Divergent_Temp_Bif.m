%% Divergent Topology Temperature Bifurcation Diagram  
%% (Grimstead et al., ,2026)
% This is the bifurcation diagram varying temperature for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026. This code is a
% revised version of the code used in Kypke et al., APS, 2026. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run and the magnitudes of the peaks in total ice
% volume calculated.
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
% See peak_Vice_Divergent.m for the calculation of peaks in volume.
%
% Base values and initial conditions recorded in p_base_values.mat
%
% Code written by Joshua Grimstead, based on the code of Kolja Kypke
%
%%
load("p_base_values.mat");

max_Vice = cell(541, 1);


t_final = 10e4;    %total time of integration

icinit =  p.ic;
total_tic = tic;
for i = 1:1021
    % choose T_s2
   % p.T_s1 = 0 + 30/1020*(i-1);
    p.T_s2 = 0 + 30/1020*(i-1);
    %p.T_s3 = 0 + 30/1020*(i-1);
    


    p.T_s2f = p.T_s2 ;
    % p.T_s1f = p.T_s1 ; 
     %p.T_s3f = p.T_s3 ;
    max_Vice_hires(i,1) =  {p.T_s2};
    
    
    % first, run for 20e4 years to get a point on the attractor
    p.tspan=[0,20e4*p.year];  
    
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
    [timeinit,Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p),p.tspan,icinit,options);  %integrate box model
    

    p.ic= Tinit(end,:);
    p.tspan=[0,p.year*t_final];     


    max_Vice_hires(i,2) = {peak_Vice_Divergent(p)};


    elapsed = toc(total_tic);
    avg_time = elapsed / i;
    est_remaining = avg_time * (1021 - i) / 60; % in minutes
    fprintf('Iter %d/1021 | T_s2: %.3f | Remaining: %.1f mins\n', i, p.T_s2, est_remaining);
    %pause(20)
end
save("max_Vice_hires_Divergent24.mat","max_Vice_hires");
%%
load("max_Vice_hires_Divergent24.mat","max_Vice_hires");


% Extract T_s2 values and peak velocities
T_s2_vals = zeros(1021,1);
Vice_vals = cell(1021,1);

for i = 1:1021
    T_s2_vals(i) = max_Vice_hires{i,1};
    Vice_vals{i} = max_Vice_hires{i,2};
end

% --- Plot 1: Bifurcation diagram (T_s2 vs peak velocity) ---
figure(); clf;
hold on;

for i = 1:1021
    v = Vice_vals{i};
    if ~isempty(v)
        % if multiple peaks (oscillatory), plot all of them
        plot(T_s2_vals(i) * ones(size(v)), v/1000, 'k.', 'MarkerSize', 1);
    end
end

xlabel('T_{s2} (-°C)', 'FontSize', 13);
ylabel('Peak of Ice Stream Volume (1000 km^3)', 'FontSize', 13);
%title('Bifurcation Diagram: Peak Ice Volume vs Surface Temperature T_{s2} (T_{s1}, T_{s3} = 15)', 'FontSize', 13);
%grid on;
%box on;
%ylim([0, 2e4]);

