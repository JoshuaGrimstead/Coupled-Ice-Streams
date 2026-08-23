%% Transient Timeseries (Grimstead et al., ,2026)
% This is the transient chaos timeseries for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026.  
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run and ice stream length or temperature is ramped
% incrementally.
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
% Base values and initial conditions recorded in p_base_values.mat
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5. 
%
% Method works as follows:
%   1. burn in at the starting value
%   2. ramp linearly over RAMP_YEARS
%   3. hold at the final value and watch the transient
%
% Time is plotted with t = 0 at RAMP START, so the ramp occupies
% [0, RAMP_YEARS/1000] kyr and the burn-in runs backwards from zero.

%% Ramping choice
% 'L1', 'L2', 'L3' ,'T_s2' 
RAMP_VAR   = 'L2';
RAMP_FROM  = 243.2e3;
RAMP_TO    = 244.2e3;
RAMP_YEARS = 50e3;


% RAMP_VAR   = 'L3';
% RAMP_FROM  = 250e3;
% RAMP_TO    = 249e3;
% RAMP_YEARS = 50e3;

%% ---- fixed parameters for this point -------------------------------
FIX_L1   = 55e3;
FIX_T_S2 = 14.02;
% L2 is the ramped variable, so its value comes from RAMP_FROM/RAMP_TO
% above rather than being fixed here.

BURN_YEARS = 5e4;    % equilibration before the ramp
HOLD_YEARS = 30e4;    % after the ramp -- long enough to see the transient end
SHOW_BURN_KYR = 50;   % how much pre-ramp history to display

N_OUT = 40000;        % output samples across the plotted window

%% ---- set up p ------------------------------------------------------
load('p_base_values.mat', 'p');
p.L1 = FIX_L1;   p.L1f = p.L1;
%p.L2 = 243.2e3;
p.L2f = p.L2;
p.L3f = p.L3;
p.T_s1f = p.T_s1;
p.T_s2  = FIX_T_S2;  p.T_s2f = p.T_s2;
p.T_s3f = p.T_s3;
p.rate  = 0;

p = setRampVar(p, RAMP_VAR, RAMP_FROM);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);   % same as every other run

%% Phase 1: burn in 
fprintf('Burn-in (%g kyr) ...\n', BURN_YEARS/1000);
tb = linspace(0, p.year*BURN_YEARS, round(N_OUT * SHOW_BURN_KYR*1000 / BURN_YEARS));
[t_b, X_b] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p), ...
                   tb, p.ic, opts);

%% Phase 2: ramp 
fprintf('Ramp %s: %g -> %g over %g kyr ...\n', RAMP_VAR, RAMP_FROM, RAMP_TO, RAMP_YEARS/1000);
ramp_sec = p.year * RAMP_YEARS;
rate = (RAMP_TO - RAMP_FROM) / ramp_sec;

tr = linspace(0, ramp_sec, round(N_OUT * RAMP_YEARS / (RAMP_YEARS + HOLD_YEARS)));
[t_r, X_r] = ode45(@(t,X) rampedRHS(t, X, p, RAMP_VAR, RAMP_FROM, rate, ramp_sec), ...
                   tr, X_b(end,:), opts);

%% Phase 3: hold 
fprintf('Hold (%g kyr) ...\n', HOLD_YEARS/1000);
p_hold = setRampVar(p, RAMP_VAR, RAMP_TO);
th = linspace(0, p.year*HOLD_YEARS, round(N_OUT * HOLD_YEARS / (RAMP_YEARS + HOLD_YEARS)));
[t_h, X_h] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p_hold), ...
                   th, X_r(end,:), opts);

%% create continuous timeseries
% t = 0 at ramp start.
t_kyr = [ t_b/p.year - BURN_YEARS ; ...
          t_r/p.year              ; ...
          t_h/p.year + RAMP_YEARS ] / 1000;
X_all = [X_b; X_r; X_h];

% The ramped variable's value at every output time 
val = [ repmat(RAMP_FROM, numel(t_b), 1) ; ...
        RAMP_FROM + rate * t_r           ; ...
        repmat(RAMP_TO,   numel(t_h), 1) ];

L1v = repmat(p.L1, numel(t_kyr), 1);
L2v = repmat(p.L2, numel(t_kyr), 1);
L3v = repmat(p.L3, numel(t_kyr), 1);
switch upper(RAMP_VAR)
    case 'L1'; L1v = val;
    case 'L2'; L2v = val;
    case 'L3'; L3v = val;
end

V_tot = ( X_all(:,1).*L1v.*p.W1 + ...
          X_all(:,5).*L2v.*p.W2 + ...
          X_all(:,9).*L3v.*p.W3 ) / 1e9;   % km^3

%% Plot
figure('Color','w','Position',[100 100 1100 450]);
hold on

ramp_end_kyr = RAMP_YEARS/1000;

% shade the ramp interval
yl = [min(V_tot) max(V_tot)];
pad = 0.04*range(yl);
yl = [yl(1)-pad, yl(2)+pad];
patch([0 ramp_end_kyr ramp_end_kyr 0], [yl(1) yl(1) yl(2) yl(2)], ...
      [0.95 0.9 0.75], 'EdgeColor','none', 'FaceAlpha', 0.7);

plot(t_kyr, V_tot/1000, 'k', 'LineWidth', 0.9);

xline(0,            'r-', 'Ramp Start', 'LineWidth', 1.8, ...
      'LabelVerticalAlignment','bottom', 'FontSize', 12);
xline(ramp_end_kyr, 'b-', 'Ramp End',   'LineWidth', 1.8, ...
      'LabelVerticalAlignment','bottom', 'FontSize', 12);

xlim([-SHOW_BURN_KYR, (RAMP_YEARS + HOLD_YEARS)/1000]);
ylim(yl/1000);
xlabel('Time Relative to Ramp Start (kyr)', 'FontSize', 13);
ylabel('Total Ice Stream Volume (1000 km^3)', 'FontSize', 13);
%title(sprintf('%s ramped %g \\rightarrow %g km over %g kyr   (L1 = %g km, T_{s2} = %g)', ...
     % RAMP_VAR, RAMP_FROM/1e3, RAMP_TO/1e3, RAMP_YEARS/1000, ...
      %FIX_L1/1e3, FIX_T_S2), 'FontSize', 13);
set(gca,'FontSize', 12, 'Layer', 'top');
box on


%% Functions
function dXdt = rampedRHS(t, X, p_base, var_name, v0, rate, ramp_sec)
% p_base is captured by value, so the ramped value cannot leak anywhere
% outside this call.
    t_c = min(max(t, 0), ramp_sec);
    p_local = setRampVar(p_base, var_name, v0 + rate*t_c);
    dXdt = Ice_Stream_Box_Model_RHS_Divergent(t, X, p_local);
end


function p = setRampVar(p, var_name, v)
% Sets both the parameter and its 'f' partner. The model interpolates
% from X toward Xf at p.rate; with p.rate = 0 that mechanism is inert,
% but both are set anyway so the ramp cannot be partially applied.
    switch upper(var_name)
        case 'L1';   p.L1 = v;   p.L1f = v;
        case 'L2';   p.L2 = v;   p.L2f = v;
        case 'L3';   p.L3 = v;   p.L3f = v;
        case 'T_S2'; p.T_s2 = v; p.T_s2f = v;
        otherwise
            error('setRampVar:unknown', 'Unsupported ramp variable "%s".', var_name);
    end
end