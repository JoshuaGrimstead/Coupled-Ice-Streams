%% Arnold tongues for the generic circle map
% This is the generic arnold tongue code used in the MSc disseration of 
% Joshua Grimstead. 
%
% Code by Joshua Grimstead, revisions made using Claude Opus 5
%
% Parameter space: (a, r), with critical threshold a = 1.
%
% Lift:   F(x) = x + (a/(2*pi))*sin(2*pi*x) + r
% Circle: f(x) = mod(F(x), 1)
% Deriv:  F'(x) = 1 + a*cos(2*pi*x)  => min = 1 - a
% So: invertible for a < 1, critical at 1, non-invertible for a > 1.

%clear; close all; clc;


%% Settings
a_min = 0.0;
a_max = 5.0;     % max coupling strength
Na    = 1001;     % resolution in a

r_min = 0.0; % range of frequency detuning (Omega in diss)
r_max = 1.0;
Nr    = 1001;     % resolution in r

x0     = 0.1234; %initial condition
Nburn  = 20000; %burn in iterations
Nkeep  = 2000; %kept iterations

Qmax   = 32;        % max denominator for tongue labelling
rhoTol = 2.5e-3;    % tolerance to declare rho locked to p/q



%% Map
F  = @(x,a,r) x + (a/(2*pi)).*sin(2*pi*x) + r;
dF = @(x,a,~) 1 + a.*cos(2*pi*x);

%% Grid
avec = linspace(a_min, a_max, Na); %coupling strength vector
rvec = linspace(r_min, r_max, Nr); %frequency detuning vector

rhoGrid   = nan(Na, Nr);
lyapGrid  = nan(Na, Nr);


%% MAIN LOOP
for ia = 1:Na
    a = avec(ia);
% compute derivatives
    for ir = 1:Nr
        r = rvec(ir);
        %x0 = rand;
        if a > 1
            x0 = (acos(-1/a))/(2*pi);
        end
  
        [rho_hat, lyap_hat] = rotation_and_lyapunov(F, dF, x0, a, r, Nburn, Nkeep);
        rhoGrid(ia,ir)  = rho_hat;
        lyapGrid(ia,ir) = lyap_hat;

        
    end

    if mod(ia, max(1,round(Na/10)))==0
        fprintf('Progress: %d / %d a rows (%.0f%%)\n', ia, Na, 100*ia/Na);
    end
end

%% Fig: Arnold tongues
Qmax = 32;

C       = nan(size(rhoGrid));
bestErr = inf(size(rhoGrid));

for q = 1:Qmax
    p      = round(rhoGrid*q);
    err    = abs(rhoGrid - p./q);
    better = err < bestErr;                 % strict: lowest q wins ties

    bestErr(better) = err(better);
    C(better)       = p(better)./q;
end

figure('Name','Arnold tongues');
imagesc(rvec, avec, C);
set(gca,'YDir','normal');
xlabel('\Omega'); ylabel('a');
cb = colorbar; ylabel(cb, '\rho');
clim([min(C(:)) max(C(:))]);
hold on;
plot([r_min r_max], [1 1], 'k--', 'LineWidth', 1.2);

%% Fig: Lyapunov exponent

% colormap with a hard break at zero
n = 256;
lo = -1; hi = 1;                      % must match clim
frac = (0.005 - lo) / (hi - lo);          % fraction of range that is negative
nneg = round(n * frac);
npos = n - nneg;

%define both colour maps in the split
neg = [linspace(0.24, 0.05, nneg)', ...
       linspace(0.15, 0.75, nneg)', ...
       linspace(0.66, 0.95, nneg)'];


pos = [linspace(0.98, 0.20, npos)', ...
       linspace(0.95, 0.75, npos)', ...
       linspace(0.10, 0.65, npos)'];

lyapmap = [neg; pos];
% make figure
figure('Name','Lyapunov exponent');
imagesc(rvec, avec, lyapGrid);
set(gca,'YDir','normal');
xlabel('\Omega'); ylabel('a');
%title('\lambda): -ive contracting, +ive expanding');
ylabel(colorbar, 'Lyapunov Exponent Clipped \lambda = [-1,1]');
colormap(gca, lyapmap);
clim([lo hi]);
hold on;
plot([r_min r_max], [1 1], 'k--', 'LineWidth', 1.2);

%% rotation and lyapunov functions 

function [rho_hat, lyap_hat] = rotation_and_lyapunov(F, dF, x0, a, r, Nburn, Nkeep)
    th = mod(x0,1);

    for n = 1:Nburn
        th = mod(F(th,a,r), 1);
    end

    total  = 0;
    logder = zeros(Nkeep,1);

    for n = 1:Nkeep
        xn    = F(th,a,r);
        total = total + (xn - th);      % lift displacement

        d = dF(th,a,r);
        if abs(d) < 1e-14, d = 1e-14; end
        logder(n) = log(abs(d));

        th = mod(xn,1);
    end

    rho_hat  = total / Nkeep;           % no mod
    lyap_hat = mean(logder);
end







