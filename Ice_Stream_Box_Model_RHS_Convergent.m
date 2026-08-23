function rhs=Ice_Stream_Box_Model_RHS_Convergent(t,X,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convergent Topology Ice Stream Box Model (Grimstead et al., ,2026)
% This is the dynamics function for the convergent ice stream box model 
% described in the MSc dissertation of Joshua Grimstead and subsequently 
% Grimstead et al., ,2026. This code is an adapted version of the 
% diveregent code used in Kypke et al., APS, 2026 which is based on the 
% original box model described in Robel et al., JGR, 2013. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well)
%
% See Ice_Stream_Box_Model_Convergent_Timeseries.m for model driver.
%
% Original code written by Alex Robel
% Debugging help from Joshua Grimstead, Christian Schoof, Eli Tziperman, 
% Eric DeGiuli and Elisa Mantelli
%
% Coupling initially added by Kolja Kypke
%
% Convergent Coupling Model by Joshua Grimstead, last updated August, 2026


%% Read in variables
h1 = X(1);
e1 = X(2); 
h_till1 = X(3);
T_b1 = X(4);

h2 = X(5);
e2 = X(6); 
h_till2 = X(7);
T_b2 = X(8);

h3 = X(9);
e3 = X(10); 
h_till3 = X(11);
T_b3 = X(12);


%% Set parameters
n=p.n;

rho_i = p.rho_i;
L_f = p.L_f;
K_i = p.K_i;
A_f = p.A_f;
g = p.g;

e_c = p.e_c;
c = p.c;
C_ice = p.C_ice;

L1 = p.L1;
W1 = p.W1;

L2 = p.L2;
W2 = p.W2;

L3 = p.L3;
W3 = p.W3;

eta_b = p.eta_b;
h_t_min = p.h_t_min;


T_s1=p.T_s1;
T_s2=p.T_s2;
T_s3=p.T_s3;


%% Diagnostic model equations 1
q_g1 = p.q_g1;
deltaT1 = T_s1-T_b1;               %difference between basal and surface ice temperature
e1 = max([e1,e_c]);            %make sure till void ratio doesn't go below threshold
h_till1 = max([h_till1,h_t_min]);h_till1=min([h_till1 p.htmax1]);   %make sure unfrozen till thickness stays within bounds
T_b1 = max([T_b1 0]);             %make sure basal ice temperature never goes above zero

%set min for h
h1 = max([h1, p.hmin]);
h2 = max([h2, p.hmin]);
h3 = max([h3, p.hmin]);

% slope12 = h2/L2;  
% slope13 = h3/L3;  
% slope1 = (W2/W1*(h1-h2) + W3/W1*(h1-h3))/L1;


slope13 = (h1-(W1/W3)*(h3))/L1;
slope23 = (h2 -(W2/W3)*(h3))/L2;
slope3 = h3/L3;

tau_d1 = rho_i*g*h1*(slope13);
tau_f1 = p.tau0*exp(-c*e1);         %calculate basal shear stress
U1 = (A_f/256)*((W1)^(n+1))*((max([(tau_d1-tau_f1)/h1,0]))^n); %calculate centerline ice stream velocity

a1 = p.a1 ;
%% Diagnostic model equations 2
q_g2 = p.q_g2;
deltaT2 = T_s2-T_b2;               %difference between basal and surface ice temperature
e2 = max([e2,e_c]);               %make sure till void ratio doesn't go below threshold
h_till2 = max([h_till2,h_t_min]);h_till2=min([h_till2 p.htmax2]);   %make sure unfrozen till thickness stays within bounds
T_b2 = max([T_b2 0]);             %make sure basal ice temperature never goes above zero

tau_d2 = rho_i*g*h2*(slope23);       %calculate driving stress
tau_f2 = p.tau0*exp(-c*e2);         %calculate basal shear stress

U2 = (A_f/256)*((W2)^(n+1))*((max([(tau_d2-tau_f2)/h2,0]))^n); %calculate centerline ice stream velocity
a2 = p.a2 ;
%% Diagnostic model equations 3
q_g3 = p.q_g3;
deltaT3 = T_s3-T_b3;               %difference between basal and surface ice temperature
e3 = max([e3,e_c]);               %make sure till void ratio doesn't go below threshold
h_till3 = max([h_till3,h_t_min]);h_till3=min([h_till3 p.htmax3]);   %make sure unfrozen till thickness stays within bounds
T_b3 = max([T_b3 0]);             %make sure basal ice temperature never goes above zero


tau_d3 = rho_i*g*h3*(slope3);        %calculate driving stress
tau_f3 = p.tau0*exp(-c*e3);         %calculate basal shear stress

U3 = (A_f/256)*((W3)^(n+1))*((max([(tau_d3-tau_f3)/h3,0]))^n); %calculate centerline ice stream velocity
a3 = p.a3 ;
%% Prognostic model equations (see Robel et al., JGR, 2013 for details) 1

if((h_till1==h_t_min && T_b1==0 && ((tau_f1*U1) + q_g1 - (K_i*deltaT1/h1))<0) || (h_till1==h_t_min && T_b1>0))   %till is frozen
    U1=0;
    dedt1 = 0;
    dhtdt1 = 0;
    dTbdt1 = -1*(1/(eta_b*C_ice))*((tau_f1*U1) + q_g1 - (K_i*deltaT1/h1));
    
    else if((e1 == e_c && h_till1==p.htmax1 && ((tau_f1*U1) + q_g1 - (K_i*deltaT1/h1))<0) ||  (e1 == e_c && h_till1<p.htmax1))   %basal ice is temperate, till is partially frozen
        U1=0;
        dedt1 = 0;
        dhtdt1 = ((tau_f1*U1) + q_g1 - (K_i*deltaT1/h1))/(e_c.*L_f*rho_i);
        dTbdt1 = 0;
        
        else
            dedt1 = ((tau_f1*U1) + q_g1 - (K_i*deltaT1/h1))/(h_till1*L_f*rho_i);
            dhtdt1 = 0;
            dTbdt1 = 0;
    end
end

dhdt1 = a1 - h1*U1/L1;

%% Prognostic model equations (see Robel et al., JGR, 2013 for details) 2
% if h2<=h1 && dhdt1<0
%     Hflux2 = -(W2/W1)*(L1*W1*a1 - W1*(h1-h2)*U1)/(L2*W2);
% else 
%     Hflux2 = 0;
% end

if((h_till2==h_t_min && T_b2==0 && ((tau_f2*U2) + q_g2 - (K_i*deltaT2/h2))<0) || (h_till2==h_t_min && T_b2>0))   %till is frozen
    U2=0;
    dedt2 = 0;
    dhtdt2 = 0;
    dTbdt2 = -1*(1/(eta_b*C_ice))*((tau_f2*U2) + q_g2 - (K_i*deltaT2/h2));
    
    else if((e2 == e_c && h_till2==p.htmax2 && ((tau_f2*U2) + q_g2 - (K_i*deltaT2/h2))<0) ||  (e2 == e_c && h_till2<p.htmax2))   %basal ice is temperate, till is partially frozen
        U2=0;
        dedt2 = 0;
        dhtdt2 = ((tau_f2*U2) + q_g2 - (K_i*deltaT2/h2))/(L_f*rho_i*e_c);
        dTbdt2 = 0;
        
        else %basal ice is temperate, till is thawed, void ratio varies
            dedt2 = ((tau_f2*U2) + q_g2 - (K_i*deltaT2/h2))/(h_till2*L_f*rho_i) ;
            dhtdt2 = 0;
            dTbdt2 = 0;
    end
end


dhdt2 = a2 - h2*U2/L2;

%% Prognostic model equations (see Robel et al., JGR, 2013 for details) 3

% if h3<=h1 && h3<=h2 && dhdt1<0 && dhdt2<0
%     Hflux3 = -(W1/W3)*(L3*W3*a3 - W3*(h1-h3)*U3)/(W1*L1)-(W2/W3)*(L3*W3*a3 - W3*(h2-h3)*U3)/(W2*L2);
% else 
%     Hflux3 = 0;
% end


if((h_till3==h_t_min && T_b3==0 && ((tau_f3*U3) + q_g3 - (K_i*deltaT3/h3))<0) || (h_till3==h_t_min && T_b3>0))   %till is frozen
    U3=0;
    dedt3 = 0;
    dhtdt3 = 0;
    dTbdt3 = -1*(1/(eta_b*C_ice))*((tau_f3*U3) + q_g3 - (K_i*deltaT3/h3));
    
    else if((e3 == e_c && h_till3==p.htmax3 && ((tau_f3*U3) + q_g3 - (K_i*deltaT3/h3))<0) ||  (e3 == e_c && h_till3<p.htmax3))   %basal ice is temperate, till is partially frozen
        U3=0;
        dedt3 = 0;
        dhtdt3 = ((tau_f3*U3) + q_g3 - (K_i*deltaT3/h3))/(L_f*rho_i*e_c);
        dTbdt3 = 0;
        
        else  %basal ice is temperate, till is thawed, void ratio varies
            dedt3 = ((tau_f3*U3) + q_g3 - (K_i*deltaT3/h3))/(h_till3*L_f*rho_i);
            dhtdt3 = 0;
            dTbdt3 = 0;
    end
end

dhdt3 = a3 - (h3*U3/L3) + W1/W3*(h1*U1/L3) + W2/W3*(h2*U2/L3);

%%
%disp(['% of integration finished: ' num2str(100*(t/p.tspan(2))) '%']);

%%
rhs= [dhdt1 ;
    dedt1
    dhtdt1
    dTbdt1
    dhdt2
    dedt2
    dhtdt2
    dTbdt2
    dhdt3
    dedt3
    dhtdt3
    dTbdt3];