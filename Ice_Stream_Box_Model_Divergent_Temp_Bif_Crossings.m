%% Divergent Topology Temperature Bifurcation Diagram with Switches
%% (Grimstead et al., ,2026)
% This is the bifurcation diagram varying temperature for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026. This code is a
% revised version of the code used in Kypke et al., APS, 2026. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run while varying Box 2 temperature and the 
% magnitudes of the peaks in total ice volume calculated.
% Crossings between thawed, partially frozen and frozen are recorded.
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
% See peak_Vice_Divergent.m for the calculation of peaks in volume.
%
% Base values and initial conditions recorded in p_base_values.mat
%
% Code written by Joshua Grimstead, based on the code of Kolja Kypke

%%
load("p_base_values.mat");

t_final = 10e4;    %total time of integration

icinit =  p.ic;
total_tic = tic;


% Pre-allocate columns: 
% col 1: current T_s2
% col 2: peak_Vice_Divergent
% col 3: T_s2 value if X(6) crossed (or NaN)
% col 4: T_s2 value if X(7) crossed (or NaN)
max_Vice_hires = cell(1021, 4); 

for i = 1:1021
    p.T_s2 = 0 + 30/1020*(i-1);
    p.T_s2f = p.T_s2;
    max_Vice_hires{i,1} = p.T_s2;
    
    % 1. Burn in run
    options_spin = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);      
    [~, Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p), [0, 20e4*p.year], icinit, options_spin);  
    
    p.ic = Tinit(end,:);
    p.tspan = [0, p.year * t_final];     
    
    % 2. Production Run with Events 
    options_event = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @myEvents);
    [time, T, te, ye, ie] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p), p.tspan, p.ic, options_event);

    % Initialize as NaN 
    X6_crossing_temp = NaN;
    X7_crossing_temp = NaN;
    
    % Check for Event 1 void ratio of box 2 meets consolidation threshold
    idx1 = find(ie == 1, 1, 'first');
    if ~isempty(idx1)
        
        X6_crossing_temp = p.T_s2; 
    end
    
    % Check for Event 2 till height reaches zero
    idx2 = find(ie == 2, 1, 'first');
    if ~isempty(idx2)
        X7_crossing_temp = p.T_s2;
    end
    
    % Save data into the matrix
    max_Vice_hires{i,2} = peak_Vice_Divergent(p); 
    max_Vice_hires{i,3} = X6_crossing_temp;   
    max_Vice_hires{i,4} = X7_crossing_temp;   
    
    % Timing metrics
    elapsed = toc(total_tic);
    avg_time = elapsed / i;
    est_remaining = avg_time * (1021 - i) / 60; 
    fprintf('Iter %d/1021 | T_s2: %.3f | Remaining: %.1f mins\n', i, p.T_s2, est_remaining);
end

save("max_Vice_hires_Divergent21.mat", "max_Vice_hires");
%%
load("max_Vice_hires_Divergent21.mat","max_Vice_hires");

% Extract T_s2 values, peak velocities, and event crossings
T_s2_vals = zeros(1021,1);
Vice_vals = cell(1021,1);
weak_to_strong = zeros(1021,1);
full_freeze = zeros(1021,1);

for i = 1:1021
    T_s2_vals(i) = max_Vice_hires{i,1};
    Vice_vals{i} = max_Vice_hires{i,2};
    
    % Store crossing states
    weak_to_strong(i) = max_Vice_hires{i,3}; 
    full_freeze(i)    = max_Vice_hires{i,4};
end

% Plot Bifurcation diagrams
figure(); clf;
hold on;

y_limits = [5, 14]; 

% 2. Plot translucent overlays
% We identify non-NaN indices where the events occurred
idx_weak = ~isnan(weak_to_strong);
idx_freeze = ~isnan(full_freeze);

% Plot red overlay for X(6) crossing (Weak-to-strong binge-purge)
if any(idx_weak)
    T_weak = T_s2_vals(idx_weak);
    % Create a shaded patch spanning the X(6) event region
    x_patch = [min(T_weak), max(T_weak), max(T_weak), min(T_weak)];
    y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
    fill(x_patch, y_patch, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Box 2 Partially Frozen');
end

% Plot blue overlay for X(7) crossing (Full till freeze on transition)
if any(idx_freeze)
    T_freeze = T_s2_vals(idx_freeze);
    % Create a shaded patch spanning the X(7) event region
    x_patch2 = [min(T_freeze), max(T_freeze), max(T_freeze), min(T_freeze)];
    y_patch2 = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
    fill(x_patch2, y_patch2, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Box 2 Frozen');
end

% Plot the Bifurcation points on top of the overlays
for i = 1:1021
    v = Vice_vals{i};
    if ~isempty(v)
        % Plotting bifurcation points in black
        plot(T_s2_vals(i) * ones(size(v)), v/1000, 'k.', 'MarkerSize', 1, 'HandleVisibility', 'off');
    end
end

% --- Formatting ---
ylim(y_limits);
xlim([min(T_s2_vals), max(T_s2_vals)]);
xlabel('T_{s2} (-°C)', 'FontSize', 13);
ylabel('Peak of Ice Stream Volume (1000 km^3)', 'FontSize', 13);
legend('show', 'Location', 'best');
box on;
%title('Bifurcation Diagram: Peak Ice Volume vs Surface Temperature T_{s2} (T_{s1}, T_{s3} = 15)', 'FontSize', 13);
%grid on;
%box on;
%ylim([0, 2e4]);


function [value, isterminal, direction] = myEvents(t, X)
   
    value = [X(6) - 0.3;     % Event 1: X(6) crosses 0.3
             X(7) - 0.01];    % Event 2: X(7) crosses 0.01
    isterminal = [0; 0];    
    direction  = [0; 0];     
end