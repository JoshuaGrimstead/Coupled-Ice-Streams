%% Poincare Section Plots (Grimstead et al., ,2026)
% This is the Poincare section code for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026.  
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run and void ratio crossings over 0.6 are recorded
% using the Events function.
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
% Base values and initial conditions recorded in p_base_values.mat
%
% Code written by Joshua Grimstead, based on the code of Kolja Kypke,
% revisions made using Claude Opus 5. 
%
%% Set up
load("p_base_values.mat", "p"); 

Ts2val_combined_v2 = [14.02, 10.2705, 10.27055, 10.2706];%[10.149, 10.150, 10.15058, 10.152];%[19.01,13.68,14.8,19.51];%[13.251, 13.252, 13.562, 13.63];[19.01,13.68,14.8,19.51];%[13.251, 13.252, 13.562, 13.63];%
L1vals_combined = [55e3,30e3,30e3,30e3];%[55e3,55e3,55e3,55e3]; [58e3,71.6e3,54.8e3,35.5e3];

nRuns = numel(Ts2val_combined_v2);

% Preallocate storage so each run's data survives the loop
V_total_all = cell(1, nRuns);
e3_all      = cell(1, nRuns);

for i = 1:nRuns
    p.e_c     = 0.3;      % till consolidation void ratio
p.hmin = 10;        % minimum till height (renamed from "hmin" so it matches
                       % what myEventsLive actually reads: p_live.h_t_min)
p.c       = 21.7;      % empirical till exponent
p.tau0    = 9.44e8;    % empirical till coefficient
p.L2  = 215e3;
p.L2f = p.L2;

p_live = p;               % loading parameters so each iteration has same ics
p_live.L1  = L1vals_combined(i);
p_live.L1f = p_live.L1;

t_final = 20e4;

    p_live.T_s2  = Ts2val_combined_v2(i);
    p_live.T_s2f = p_live.T_s2;   % FIXED: was p_live.Ts2 (typo, undefined field)

    % Rebuild options each iteration so the events fn closes over current p_live
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                      'Events', @(t,X) myEventsLive(t, X, p_live));

    fprintf('Run %d/%d - Step 1/2: burn-in at T_s2 = %.4f ...\n', i, nRuns, p_live.T_s2);
    [~, Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p_live), ...
                       [0, p_live.year*t_final], p_live.ic, options);

    p_live.ic = Tinit(end,:);
    t_final_prod = t_final * 20; % long enough for dense Poincare points

    fprintf('Run %d/%d - Step 2/2: production run...\n', i, nRuns);
    tspan = [0, p_live.year*t_final_prod];
    [t, T, te, ye, ie] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p_live), ...
                               tspan, p_live.ic, options);

    % Poincare section using e2 threshold crossings (event 1)
    e2_indices = (ie == 1);
    if any(e2_indices)
        V1_e = (ye(e2_indices, 1) * p_live.L1 * p_live.W1) / 1e9;
        V2_e = (ye(e2_indices, 5) * p_live.L2 * p_live.W2) / 1e9;
        V3_e = (ye(e2_indices, 9) * p_live.L3 * p_live.W3) / 1e9;
        V_total_all{i} = V1_e + V2_e + V3_e;
        e3_all{i}      = ye(e2_indices, 10);
    else
        V_total_all{i} = [];
        e3_all{i}      = [];
    end
end

save('poincare_results_10.mat', 'V_total_all', 'e3_all', 'Ts2val_combined_v2', 'L1vals_combined', 'p_live');



%% reload and plot
S = load('poincare_results_10.mat', 'V_total_all', 'e3_all', 'Ts2val_combined_v2', 'L1vals_combined', 'p_live');
plot_poincare_2x2(S.V_total_all, S.e3_all, S.Ts2val_combined_v2, S.L1vals_combined);

%% functions
function [value, isterminal, direction] = myEventsLive(t, X, p_live)
    thresh = 0.6;
    value = [X(6)  - thresh;           % Event 1: e2 cross 0.6
             X(10) - thresh;           % Event 2: e3 cross 0.6
             X(7)  - p_live.h_t_min];  % Event 3: till height at minimum
    isterminal = [0; 0; 0];
    direction  = [1; 1; -1];           % only positive crossings for 1,2;
    %  negative for 3
end


function plot_poincare_2x2(V_total_all, e3_all, Ts2vals, L1vals)
    % axis limits are optionally global to ensure comparable plots
    xAll = []; yAll = [];
    for m = 1:numel(V_total_all)
        xAll = [xAll; V_total_all{m}(:) ./ 1000]; 
        yAll = [yAll; e3_all{m}(:)];              
    end
 
    if ~isempty(xAll)
        xPad = 0.05 * range(xAll); if xPad == 0, xPad = 1; end
        yPad = 0.05 * range(yAll); if yPad == 0, yPad = 0.01; end
        xLims = [min(xAll)-xPad, max(xAll)+xPad];
        yLims = [min(yAll)-yPad, max(yAll)+yPad];
    else
        xLims = [0 1]; yLims = [0 1]; % fallback if every run was empty
    end
 
    figure('Color', 'w', 'Name', 'Poincare Sections', 'Position', [1150 100 600 500]);
    tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
    axHandles = gobjects(1, numel(Ts2vals));
    for m = 1:numel(Ts2vals)
        axHandles(m) = nexttile;
        if ~isempty(V_total_all{m})
            scatter(V_total_all{m} ./ 1000, e3_all{m}, 20, 'filled', 'MarkerFaceAlpha', 0.6);
        else
            text(0.5, 0.5, 'No crossings', 'HorizontalAlignment', 'center');
        end
        title(sprintf('T_{s2} = %.3f, L_1 = %.1f km', Ts2vals(m), L1vals(m)/1000));
       % xlim(xLims);
        %ylim(yLims);
        ax = gca;
        ax.FontSize = 12;
        grid off;
    end
    %allows zooming to be simultaneous across plots
    %linkaxes(axHandles, 'xy');
 
    xlabel(tl, 'Total Ice Stream Volume (1000 km^3)', 'FontSize', 16);
    ylabel(tl, 'Box 3 Void Ratio', 'FontSize', 16);
end
 