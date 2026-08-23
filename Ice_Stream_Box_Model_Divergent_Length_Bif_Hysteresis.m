%% Divergent Topology Length Bifurcation Diagram (Grimstead et al., ,2026)
% This is the bifurcation diagram varying temperature for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026. This code is a
% revised version of the code used in Kypke et al., APS, 2026. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). Simulations are run whilst box 2 length is varied, with initial
% conditions carried forward to test hysteresis. The magnitudes of the 
% peaks in total ice volume calculated and the initial conditions at the
% end of the run recorded. 
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



% Pre-allocate cell array
bifurc_L2_inc = cell(201, 2);

%% Increasing Run
tic
for i = 1:201
    iter_start = tic;

    p.L2 =  160e3 + 140e3/200*(i-1);


    bifurc_L2_inc(i,1) =  {p.L2};
 
    % run for 20e4 years to get a point on the attractor
    p.tspan=[0,20e4*p.year];  
   
   % run ODE
    options = odeset('RelTol',1e-6,'AbsTol',1e-6, 'OutputFcn', []);     %set ode integration settings
    [timeinit,Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p),[0,20e4*p.year],p.ic,options);  %integrate box model
    %

    p.ic= Tinit(end,:);
    p.tspan=[0,p.year*t_final];     


    bifurc_L2_inc(i,2) = {peak_Vice_Divergent(p)};
    %
    iter_time = toc(iter_start);
    total_time = toc;
    avg_time = total_time / i;
    remaining = avg_time * (201 - i);
    fprintf('Iter %d/201 | L2=%.1fkm | This: %.1fmin | Avg: %.1fmin | Left: %.1fhr\n', ...
            i, p.L2/1000, iter_time/60, avg_time/60, remaining/3600);
    icinit = Tinit(end,:);
end
ic_from_inc = Tinit(end,:);
save("bifurc_L2_inc_vol_fix7.mat","bifurc_L2_inc");
save("inc_ic.mat", "ic_from_inc")
%%
load("bifurc_L2_inc_vol_fix7.mat","bifurc_L2_inc");

figure()
clf
hold on
for i =1:201
   scatter( bifurc_L2_inc{i,1}/1000*ones(1,length(bifurc_L2_inc{i,2})), bifurc_L2_inc{i,2} , 1, 'filled', 'MarkerFaceColor', 'black'  ) 
end
%xlim([180,240])
%ylim([1.15e4,1.45e4])
ax = gca;
ax.FontSize = 16; 
xlabel("Box 2 length (km)")
ylabel("Peak total ice volume (km^3)")
drawnow  
%% Decreasing Run
load("p_base_values.mat");
load("inc_ic.mat", "ic_from_inc")
icinit = ic_from_inc;

% Pre-allocate cell array
%%
bifurc_L2_dec = cell(201, 2);

%% 
tic
for i = 1:201
    iter_start = tic;

    p.L2 =  300e3 - 140e3/200*(i-1);


    bifurc_L2_dec(i,1) =  {p.L2};
    
    p.tspan=[0,20e4*p.year];  
   


    options = odeset('RelTol',1e-6,'AbsTol',1e-6, 'OutputFcn', []);     %set ode integration settings
    [timeinit,Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p),[0,20e4*p.year],icinit,options);  %integrate box model
   

    p.ic= Tinit(end,:);
    p.tspan=[0,p.year*t_final];     


    bifurc_L2_dec(i,2) = {peak_Vice_Divergent(p)};
   
    iter_time = toc(iter_start);
    total_time = toc;
    avg_time = total_time / i;
    remaining = avg_time * (201 - i);
    fprintf('Iter %d/201 | L2=%.1fkm | This: %.1fmin | Avg: %.1fmin | Left: %.1fhr\n', ...
            i, p.L2/1000, iter_time/60, avg_time/60, remaining/3600);
    icinit = Tinit(end,:);
 end
save("bifurc_L2_dec_vol_fix7_inc_ic.mat","bifurc_L2_dec");
%%
load("bifurc_L2_dec_vol_fix7_inc_ic.mat","bifurc_L2_dec");

figure()
clf
hold on
for i =1:201
   scatter( bifurc_L2_dec{i,1}/1000*ones(1,length(bifurc_L2_dec{i,2})), bifurc_L2_dec{i,2} , 1, 'filled', 'MarkerFaceColor', 'black'  ) 
end
%xlim([180,240])
%ylim([1.15e4,1.45e4])
ax = gca;
ax.FontSize = 16; 
xlabel("Box 2 length (km)")
ylabel("Peak total ice volume (km^3)")
drawnow  


%%
load("bifurc_L2_dec_vol_fix6_inc_ic.mat","bifurc_L2_dec");
load("bifurc_L2_inc_vol_fix6.mat","bifurc_L2_inc");

lightBlueRegions = [ 187  200;    % region 1 -> light blue
                     230 242];  % region 2 -> light blue

darkBlueRegions  = [ 168 187;    % region 3 -> dark blue
                     225 230;
                     212 218];  % region 4 -> dark blue

lightBlueColor = [0.65 0.82 0.95];  % light blue
darkBlueColor  = [0.05 0.20 0.55];  % darker blue
regionAlpha    = 0.35;              % transparency 

% generate figure
figure()
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');



% branch
axTop = nexttile; hold(axTop,'on')
for i = 1:201
    h1 = scatter(axTop, bifurc_L2_inc{i,1}/1000*ones(1,length(bifurc_L2_inc{i,2})), ...
        bifurc_L2_inc{i,2}/1000, 1, 'filled', 'MarkerFaceColor', 'black');
end
axTop.FontSize = 16;


title(axTop, 'Increasing Branch')

% Decreasing branch
axBottom = nexttile; hold(axBottom,'on')
for i = 1:201
    h2 = scatter(axBottom, bifurc_L2_dec{i,1}/1000*ones(1,length(bifurc_L2_dec{i,2})), ...
        bifurc_L2_dec{i,2}/1000, 1, 'filled', 'MarkerFaceColor', 'black');
end
axBottom.FontSize = 16;
xlabel(axBottom, "Box 2 length (km)")
%ylabel(axBottom, "Peak total ice volume (km^3)")
title(axBottom, 'Decreasing Branch')

%sgtitle(t, 'Bifurcation Diagram of Ice Volume vs. Box 2 Length')

% overlay both plots
axesList = [axBottom, axTop];
pfLight = gobjects(1);
pfDark  = gobjects(1);

for a = 1:numel(axesList)
    ax = axesList(a);
    yl = ylim(ax);   % lock in the data-driven y-range before adding patches

    for r = 1:size(lightBlueRegions,1)
        x0 = lightBlueRegions(r,1); x1 = lightBlueRegions(r,2);
        p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], lightBlueColor, ...
            'FaceAlpha', regionAlpha, 'EdgeColor', 'none');
        uistack(p,'top')      % render behind data points
        if a == 1, pfLight = p; end   % keep one handle for the legend
    end

    for r = 1:size(darkBlueRegions,1)
        x0 = darkBlueRegions(r,1); x1 = darkBlueRegions(r,2);
        p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], darkBlueColor, ...
            'FaceAlpha', regionAlpha, 'EdgeColor', 'none');
        uistack(p,'top')
        if a == 1, pfDark = p; end
    end

    ylim(ax, yl);  % restore in case adding patches auto-rescaled the axis
end
ylabel("            Peak Total Ice Volume (1000 km^3)")
legend(axBottom,    [pfLight, pfDark], ...
    {'Chaotic in Increasing','Chaotic in Both'}, 'Location','best')
%legend(axTop, h1, 'Increasing Branch', 'Location','best')

linkaxes(axesList, 'x')   % zooming/panning one subplot moves the other too
drawnow