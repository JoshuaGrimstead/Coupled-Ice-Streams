%% Convergent Topology Length Bifurcation Diagram (Grimstead et al., ,2026)
% This is the bifurcation diagram varying temperature for the convergent
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
% See Ice_Stream_Box_Model_RHS_Convergent.m for implementation of model 
% equations.
%
% See peak_Vice_Convergent.m for the calculation of peaks in volume.
%
% Base values and initial conditions recorded in p_base_values_conv.mat
%
% Code written by Joshua Grimstead, based on the code of Kolja Kypke


%% Set parameters

load("p_base_values_conv.mat",'p');

p.T_s2 = 15;      


p.rate = 0;		% K/yr
p.start_time = 0; 
icinit = p.ic;
bifurc_L2_inc = cell(201, 2);

%% generate bif data
tic
for i = 1:201
    iter_start = tic;

    p.L2 =  180e3 + 60e3/200*(i-1);


    bifurc_L2_inc{i,1} =  {p.L2};
%     
    % first, run for 20e4 years to get a point on the attractor
    p.tspan=[0,20e4*p.year];  
   % icinit = Tinit(end,:);
   
  
    options = odeset('RelTol',1e-6,'AbsTol',1e-6, 'OutputFcn', []);     %set ode integration settings
     [time, T] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Convergent(t,X,p), p.tspan, icinit, options);

   
    
    p.ic= T(end,:);
    p.tspan=[0,p.year*t_final];  

    % Save data 
    bifurc_L2_inc{i,2} = {peak_Vice_Convergent(p)};
 
    iter_time = toc(iter_start);
    total_time = toc;
    avg_time = total_time / i;
    remaining = avg_time * (201 - i);
    fprintf('Iter %d/201 | L2=%.1fkm | This: %.1fmin | Avg: %.1fmin | Left: %.1fhr\n', ...
            i, p.L2/1000, iter_time/60, avg_time/60, remaining/3600);
    icinit = T(end,:);
end

ic_from_inc_st = icinit;
save("bifurc_terminus_L2_inc_3.mat","bifurc_L2_inc");
save("inc_ic_st3.mat", "ic_from_inc_st")
%%
load("bifurc_terminus_L2_inc_3.mat","bifurc_L2_inc");
figure()
clf
hold on
for i =1:201
   scatter( bifurc_L2_inc{i,1}{1}/1000*ones(1,length(bifurc_L2_inc{i,2}{1})), bifurc_L2_inc{i,2}{1} , 1, 'filled', 'MarkerFaceColor', 'black'  ) 
end
%xlim([180,240])
%ylim([1.15e4,1.45e4])
ax = gca;
ax.FontSize = 16; 
xlabel("Box 2 length (km)")
ylabel("Peak total ice volume (km^3)")
drawnow  % Force MATLAB to render the figure NOW


%% decreasing 
max_Vice = cell(541, 1);


t_final = 10e4;    %total time of integration

bifurc_L2_inc = cell(201, 2);
load("inc_ic_st.mat", "ic_from_inc_st")
icinit = ic_from_inc_st;
%% generate points on bifurc diagram
tic
for i = 1:201
    iter_start = tic;
% choose decreasing L2 values
    p.L2 =  240e3 - 60e3/200*(i-1);
    bifurc_L2_inc(i,1) =  {p.L2};

    %length of run
    p.tspan=[0,20e4*p.year];  
   

  %run model
    options = odeset('RelTol',1e-6,'AbsTol',1e-6, 'OutputFcn', []);     %set ode integration settings
    [timeinit,Tinit] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Convergent(t,X,p),[0,20e4*p.year],icinit,options);  %integrate box model
   

    p.ic= Tinit(end,:);
    p.tspan=[0,p.year*t_final];     


    bifurc_L2_inc(i,2) = {peak_Vice_Convergent(p)};
   
    iter_time = toc(iter_start);
    total_time = toc;
    avg_time = total_time / i;
    remaining = avg_time * (201 - i);
    fprintf('Iter %d/201 | L2=%.1fkm | This: %.1fmin | Avg: %.1fmin | Left: %.1fhr\n', ...
            i, p.L2/1000, iter_time/60, avg_time/60, remaining/3600);
    icinit = Tinit(end,:);
end

ic_from_inc_st = icinit;
save("bifurc_terminus_L2_dec_3.mat","bifurc_L2_inc");

%%
load("bifurc_terminus_L2_dec_3.mat","bifurc_L2_inc");
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
drawnow  % Force MATLAB to render the figure NOW
