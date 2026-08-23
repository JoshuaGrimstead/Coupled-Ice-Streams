%% Convergent Topology Temperature Bifurcation Diagram 
%% (Grimstead et al., ,2026)
% This is the bifurcation diagram varying temperature for the convergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026. This code is a
% revised version of the code used in Kypke et al., APS, 2026. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run and the magnitudes of the peaks in total ice
% volume calculated.
%
% See Ice_Stream_Box_Model_RHS_Convergent.m for implementation of model 
% equations.
%
% See peak_Vice_Convergent.m for the calculation of peaks in volume.
%
% Base values and initial conditions recorded in p_base_values_conv.mat
%
% Code written by Joshua Grimstead, based on the code of Kolja Kypke
%
% 
%
%
%% Set parameters
load("p_base_values_conv.mat",'p'); %load base values
% p.ic = [400.857202814214, 0.646369044549, 1.00002796940515, -3.79683468406686e-05, 494.189813758044, 0.324124717495546, 1.0000169094565, 7.35735153932896e-05, 484.69763013887, 0.299993534301783, 0.414729082679925, 6.08967772196069e-05 ];

 
%icinit =  [2.513464916211947e+02, 0.719606522947843e+04, 1.000027969405150, -3.796834684066860e-05, 2.675040993016955e+02, 0.683872730768198, 1.000043214013321, -3.326578858649950e-05, 4.067700935219721e+02, 1.642175807274792, 1.999885855326992, 2.683602508176260e-05];
    % 
    icinit = p.ic;
p.T_s2 = 15;      
p.L2 = 235e3; % for sufficently cold temperatures, accumulation dominates
%  unrealistically for small ice stream footprints,thus a sufficent length
%  must be chosen to match the accumulation
p.rate = 0;		% K/yr
p.start_time = 0; 
%%

%icinit = p.ic;
max_Vice = cell(541, 1);


t_final = 10e4;    %total time of integration


total_tic = tic;
for i = 1:221

    p.T_s2 = 0 + 30/220*(i-1);

    p.T_s2f = p.T_s2 ;
 
    max_Vice_hires(i,1) =  {p.T_s2};
    
    % first, run for 20e4 years to get a point on the attractor
    p.tspan=[0,20e4*p.year];  
   
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
    [timeinit,Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Convergent(t,X,p),p.tspan,icinit,options);  %integrate box model
    
    p.ic= Tinit(end,:);
    p.tspan=[0,p.year*t_final];     

    max_Vice_hires(i,2) = {peak_Vice_Convergent(p)};


    elapsed = toc(total_tic);
    avg_time = elapsed / i;
    est_remaining = avg_time * (221 - i) / 60; % in minutes
    fprintf('Iter %d/221 | T_s2: %.3f | Remaining: %.1f mins\n', i, p.T_s2, est_remaining);
    %pause(20)
end
%%
save("max_Vice_hires_Shared6.mat","max_Vice_hires");
%%
load("max_Vice_hires_Shared6.mat","max_Vice_hires");


% Extract T_s2 values and peak velocities
T_s2_vals = zeros(221,1);
Vice_vals = cell(221,1);

for i = 1:221
    T_s2_vals(i) = max_Vice_hires{i,1};
    Vice_vals{i} = max_Vice_hires{i,2};
end

% --- Plot 1: Bifurcation diagram (T_s2 vs peak velocity) ---
figure(); clf;
hold on;

for i = 1:221
    v = Vice_vals{i};
    if ~isempty(v)
        % if multiple peaks (oscillatory), plot all of them
        plot(T_s2_vals(i) * ones(size(v)), v, 'k.', 'MarkerSize', 1);
    end
end

xlabel('T_{s2} (-°C)', 'FontSize', 13);
ylabel('Peak Total Ice Volume (km^3)', 'FontSize', 13);
%title('Bifurcation Diagram: Peak Ice Volume vs Surface Temperature T_{s2} (T_{s1}, T_{s3} = 15)', 'FontSize', 13);
%grid on;
box off;
xlim([0, 30]);

