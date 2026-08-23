%% Convergent Topology Timeseries (Grimstead et al., ,2026)
% This is the dynamics function for the convergent ice stream box model
% described in the MSc dissertation of Joshua Grimstead and subsequently 
% Grimstead et al., ,2026. This code is an adapted version of the code used
% in Kypke et al., APS, 2026 which is based on the original box model 
% described in Robel et al., JGR, 2013. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well)
%
% See Ice_Stream_Box_Model_RHS_Convergent.m for implementation of model 
% equations.
%
% Base values and initial conditions recorded in p_base_values_conv.mat
%
% Original code written by Alex Robel
% Debugging help from Joshua Grimstead, Christian Schoof, Eli Tziperman,
% Eric DeGiuli and Elisa Mantelli
%
% Coupling initially added by Kolja Kypke
%
% Convergent Coupling Model by Joshua Grimstead, last updated August, 2026

%% Set initial conditions and integration time

load("p_base_values_conv.mat",'p');

p.L2 = 190.8e3;
p_test.ic = p.ic;
   % p.T_s2 = 15;
    p.tspan = [0, 3e5*p.year];
    p.ic = [400.857202814214, 0.646369044549, 1.00002796940515, -3.79683468406686e-05, 494.189813758044, 0.324124717495546, 1.0000169094565, 7.35735153932896e-05, 484.69763013887, 0.299993534301783, 0.414729082679925, 6.08967772196069e-05 ];


%% Run Model
%p.tspan=[0,p.year*t_final];     %time steps

options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
[time,T] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Convergent(t,X,p),p.tspan,p.ic,options);  %integrate box model

%% Warning for unrealistic parameter choice

if any(min(T(:,[1 5 9])) < p.hmin)
    warning('Thickness fell below hmin — box %d. Check parameters.', ...
        find(min(T(:,[1 5 9])) < p.hmin));
end

% Keep only the last 1e5 years 
keep = time <= 2e5*p.year;
time = time(keep);
T = T(keep,:);

%% Diagnostic 1
h1 = T(:,1);
e1 = T(:,2);
h_till1 = T(:,3);
T_b1 = T(:,4);
deltaT1 = p.T_s1 - T_b1;               %difference between basal and surface ice temperature
e1(e1<=p.e_c)=p.e_c;               %make sure till void ratio doesn't go below threshold
h_till1(h_till1<=p.h_t_min) = p.h_t_min; h_till1(h_till1>=p.htill_init) = p.htill_init;  
                                     %make sure unfrozen till thickness stays within bounds
T_b1(T_b1<=0)=0;                     %make sure basal ice temperature never goes above zero

%% Diagnostic 2
h2 = T(:,5);
e2 = T(:,6);
h_till2 = T(:,7);
T_b2 = T(:,8);
deltaT2 = p.T_s2 - T_b2;               %difference between basal and surface ice temperature
e2(e2<=p.e_c)=p.e_c;               %make sure till void ratio doesn't go below threshold
h_till2(h_till2<=p.h_t_min) = p.h_t_min; h_till2(h_till2>=p.htill_init) = p.htill_init; 
                                     %make sure unfrozen till thickness stays within bounds
T_b2(T_b2<=0)=0;                     %make sure basal ice temperature never goes above zero



%% Diagnostic 3
h3 = T(:,9);
e3 = T(:,10);
h_till3 = T(:,11);
T_b3 = T(:,12);
deltaT3 = p.T_s3 - T_b3;               %difference between basal and surface ice temperature
e3(e3<=p.e_c)=p.e_c;               %make sure till void ratio doesn't go below threshold
h_till3(h_till3<=p.h_t_min) = p.h_t_min; h_till3(h_till3>=p.htill_init) = p.htill_init; 
                                     %make sure unfrozen till thickness stays within bounds
T_b3(T_b3<=0)=0;                     %make sure basal ice temperature never goes above zero

%%
% mean slopes
% if h1 >= h2
%     slope13 = (h1-h3)/p.L1;
%     slope2 = h2/p.L2;
% else
%     slope23 = (h2-h3)/p.L2; 
%     slope1 = h1/p.L1;
% end
% slope13 = (h1-h3)/p.L1;
% slope23 = (h2-h3)/p.L2; 
% slope3 = h3/p.L3;

%slope1 = (h1-max(h2,h3))/p.L1;

slope13 = (h1-(p.W1/p.W3)*(h3))/p.L1;
slope23 = (h2 -(p.W2/p.W3)*(h3))/p.L2;
slope3 = h3/p.L3;
% driving, frictional stresses and basal velocities
tau_d1 = p.rho_i*p.g*h1.*(slope13);  
%calculate driving stress
tau_f1 = p.tau0*exp(-p.c*e1);         %calculate basal shear stress
U1 = (p.A_f/256)*((p.W1)^(p.n+1))*(((tau_d1-tau_f1)./h1).^p.n); U1 = max(U1,0); %calculate centerline ice stream velocity
tau_d2 = p.rho_i*p.g*h2.*(slope23);        %calculate driving stress
tau_f2 = p.tau0*exp(-p.c*e2);         %calculate basal shear stress
U2 = (p.A_f/256)*((p.W2)^(p.n+1))*(((tau_d2-tau_f2)./h2).^p.n); U2 = max(U2,0); %calculate centerline ice stream velocity
tau_d3 = p.rho_i*p.g*h3.*(slope3);        %calculate driving stress
tau_f3 = p.tau0*exp(-p.c*e3);         %calculate basal shear stress
U3 = (p.A_f/256)*((p.W3)^(p.n+1))*(((tau_d3-tau_f3)./h3).^p.n); U3 = max(U3,0); %calculate centerline ice stream velocity

V1 = (T(:,1)*p.L1*p.W1)/10^9;
V2 = (T(:,5)*p.L2*p.W2)/10^9;
V3 = (T(:,9)*p.L3*p.W3)/10^9;
V = V1+V2+V3;


%% Make some plots % %

f = figure();
%set(3,'units','pixels','position',[0 0 1002 1202])
t1=tiledlayout(4,1, "TileSpacing", "tight", "Padding","tight");
nexttile
plot(time./p.year/1000,V/1000,'linewidth',2)
title("Total Volume")
ylabel('1000 km^3');
ax = gca;
ax.FontSize = 16; 
t2 = tiledlayout(t1,2,1, "TileSpacing", "tight", "Padding","tight");
t2.Layout.Tile = 2;
t2.Layout.TileSpan = [1 1];
nexttile(t2);
plot(time./p.year/1000,(U2*p.year+1),'b','linewidth',2);
ylim([300,450])
ax = gca;
set(gca,'xtick',[])
set(gca,'ytick',[300,400])
ax.FontSize = 16; 
nexttile(t2)
hold on
plot(time./p.year/1000,(U3*p.year+1),'r','linewidth',2);
plot(time./p.year/1000,(U1*p.year+1),'k','linewidth',2);
plot(time./p.year/1000,(U2*p.year+1),'b','linewidth',2);
ylim([0,100])
title(t2,"Basal Sliding Velocity",'Fontsize', 16, 'FontWeight', 'bold');
ylabel(t2,"m/yr", 'Fontsize', 16);
ax = gca;
ax.FontSize = 16; 
nexttile(t1);
%subplot(4,1,3)
hold on;
plot(time./p.year/1000, T(:,3),'k','linewidth',2);
plot(time./p.year/1000, T(:,7),'b','linewidth',2);
plot(time./p.year/1000, T(:,11),'r','linewidth',2);
title('Unfrozen Till Thickness')
ylabel('m')
ylim([-0.1,2.1])
ax = gca;
ax.FontSize = 16; 
nexttile(t1);
hold on;
plot(time./p.year/1000, T(:,2),'k','linewidth',2);
plot(time./p.year/1000, T(:,6),'b','linewidth',2);
plot(time./p.year/1000, T(:,10),'r','linewidth',2);
title('Void Ratio');
xlabel('Time (kyr)')
ax = gca;
ax.FontSize = 16; 

