function Q = merge(geometry,P,rho_0)
% =========================================================================
% This function is used to compute the flows (or -OPTIONNAL- the densities)
% at a 2-to-1 merge (say a junction with 2 incoming roads and simply 
% 1 outgoing road)
% -------------------------------------------------------------------------
% Inputs:
% (1) key network functions encoded in the *geometry* structure:
%     (a) demand functions *Demand_i* on each incoming road
%     (b) supply function *Supply_j* on the outgoing road
%     [OPTIONAL]
%     (c) flow functions for each road
%     (d) maximal density value for each road
% (2) priority parameter *P* between both incoming roads (0<P<1)
% (3) initial condition *rho_0*: densities on each road at initial time, as
%     a vector (discrete in space)
% -------------------------------------------------------------------------
% Outputs:
% (1) the vector of all the flows *Q* defined as Q=[q1 q2 q3]
% OR
% (optional) densities *rho* on each road at next step, as a vector
% =========================================================================

% I. Initialization
Demand_1 = geometry(1).Demand;
Demand_2 = geometry(2).Demand;
Supply_3 = geometry(3).Supply;

% ------------------ Optional ------------------
% (only if one wants to compute the densities and not the flows)
% rho = NaN(size(rho_0));
% flow_1 = geometry(1).flow_function;
% flow_2 = geometry(2).flow_function;
% flow_3 = geometry(3).flow_function;
% rho_max_1 = geometry(1).rho_max;
% rho_max_2 = geometry(2).rho_max;
% rho_max_3 = geometry(3).rho_max;
% ----------------------------------------------

% Demand on the first incoming road
D1 = Demand_1(rho_0(1));
% Demand on the second incoming road
D2 = Demand_2(rho_0(2));
% Supply on the outgoing road
S3 = Supply_3(rho_0(3));

% II. Distinction of the cases + Compute the densities from the flows
if S3 >= D1 + D2
    % Demand-contrained regime
    q1 = D1;
    q2 = D2;
    q3 = D1+D2;
    
%     ------------------ Optional ------------------
%     % Free-flow densities
%     [rho(1),~] = density_from_flow(flow_1,rho_max_1,q1);
%     [rho(2),~] = density_from_flow(flow_2,rho_max_2,q2);
%     [rho(3),~] = density_from_flow(flow_3,rho_max_3,q3);    
%     ----------------------------------------------
    
else
    % Supply-contrained regime
    q1 = min(D1, max(P*S3, S3-D2) );
    q2 = min(D2, max((1-P)*S3, S3-D1) );
    q3 = q1+q2;
    
%     ------------------ Optional ------------------
%     % Congested densities
%     [~,rho(1)] = density_from_flow(flow_1,rho_max_1,q1);
%     [~,rho(2)] = density_from_flow(flow_2,rho_max_2,q2);
%     [~,rho(3)] = density_from_flow(flow_3,rho_max_3,q3);  
%     ----------------------------------------------
end

Q = [q1 q2 q3];