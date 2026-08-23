
function peaks = peak_Vice_Divergent(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Divergent Topology Ice Volume Peaks (Grimstead et al., ,2026)
% This is the volume peaks function for the divergent
% ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026. This code is a
% revised version of the code used in Kypke et al., APS, 2026. 
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well). 
%
% See Ice_Stream_Box_Model_RHS_Divergent.m for implementation of model 
% equations.
%
%
% Code adapted by Joshua Grimstead, based on the code of Kolja Kypke

options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
[time,T] = ode45(@(t,X) Ice_Stream_Box_Model_RHS_Divergent(t,X,p),p.tspan,p.ic,options);  %integrate box model

V1 = (T(:,1)*p.L1*p.W1)/10^9;
V2 = (T(:,5)*p.L2*p.W2)/10^9;
V3 = (T(:,9)*p.L3*p.W3)/10^9;
V = V1+V2+V3;

Vmax = [];
k=1;
for i=2:(length(time)-1)
    if V(i) >= V(i-1) && V(i)>= V(i+1) 
        % find time of local max (peaks) of total ice volume
       Vmax(k) = V(i);
       k=k+1;
    end 
end

peaks = Vmax;

end