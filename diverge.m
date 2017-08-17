function Q = diverge(geometry,A,rho_0)
% =========================================================================
% This function is used to compute the flows (or -OPTIONNAL- the densities)
% at a 1-to-2 diverge (say a junction with 1 incoming road and simply 
% 2 outgoing roads)
% -------------------------------------------------------------------------
% Inputs:
% (1) key network functions encoded in the *geometry* structure:
%     (a) demand function *Demand_i* on the incoming road
%     (b) supply functions *Supply_j* on each outgoing road
%     [OPTIONAL]
%     (c) flow functions for each road
%     (d) maximal density value for each road
% (2) the assignment matrix *A* that gives the drivers preferences between 
%     both outgoing roads (q2 = A q1 and q3 = (1-A) q1, with 0<A<1)
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
Supply_2 = geometry(2).Supply;
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

% Demand on the incoming road
D1 = Demand_1(rho_0(1));
% Supply on the first outgoing road
S2 = Supply_2(rho_0(2));
% Supply on the second outgoing road
S3 = Supply_3(rho_0(3));

q1 = min( D1, min(S2/A, S3/(1-A)) ) ;
q2 = A*q1 ;
q3 = (1-A)*q1;

Q = [q1 q2 q3];

% ------------------ Optional ------------------
% if q1 == D1
%     % Free-flow densities
%     [rho(1),~] = density_from_flow(flow_1,rho_max_1,q1);
%     [rho(2),~] = density_from_flow(flow_2,rho_max_2,q2);
%     [rho(3),~] = density_from_flow(flow_3,rho_max_3,q3);
% else
% % Congested densities
%     [~,rho(1)] = density_from_flow(flow_1,rho_max_1,q1);
%     [~,rho(2)] = density_from_flow(flow_2,rho_max_2,q2);
%     [~,rho(3)] = density_from_flow(flow_3,rho_max_3,q3);
% end
% ----------------------------------------------