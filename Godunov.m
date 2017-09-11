function [rho,Delta_t]=Godunov(geometry,rho_0,...
    Demand_upstream,Supply_downstream,...
    Delta_x,T)
% =========================================================================
% This function aims to solve the LWR model (1D scalar conservation law) 
% thanks to the Godunov scheme (1D finite volume scheme)
% -------------------------------------------------------------------------
% Inputs:
% (1) all the informations about the network for each homogeneous section:
%       (a) length of the section
%       (b) demand function
%       (c) supply function
%       (d) maximal speed
% (2) initial condition: densities on each road at initial time
%     /!\ WARNING : here it is a smooth function of space
%     but with real data, it should be discrete in space
% (3) demand at the upstream boundary of the network (entry)
%     /!\ WARNING : here it is a smooth function of time 
%     but with real data, it should be discrete in time
% (4) supply at the downstream boundary of the network (exit)
%     /!\ WARNING : here it is a smooth function of time 
%     but with real data, it should be discrete in time
% (5) discrete spatial step *Delta_x*
% (6) time horizon *T*
% -------------------------------------------------------------------------
% Outputs:
% (1) densities *rho* on each road at next step, as a vector
% (2) discrete time step *Delta_t* that satisfies the CFL condition
% =========================================================================


number_section = length(geometry);
V_max = 0;
X = [];
L_section = NaN(1,number_section);
L_partial = 0;

% %--------------------- Optional ---------------------
% figure; hold on
for i = 1:number_section
    V_max = max( V_max, geometry(i).Vmax) ;
    
    % Partial computation for the total length of the network
    L_section(i) = geometry(i).length;
    dim_prev_sec = floor(L_partial/Delta_x)+1 ;
    dim_curr_sec = floor((L_partial+L_section(i))/Delta_x) ;
    X(2,dim_prev_sec:dim_curr_sec) = ...
        i.*ones(1,dim_curr_sec-dim_prev_sec+1) ;
    L_partial = L_partial + L_section(i);
    
%     %--------------------- Optional ---------------------
%     ezplot(geometry(i).Demand,[0 200])
%     ezplot(geometry(i).Supply,[0 200])
end
% hold off
% axis tight
% xlabel('Density (veh/km)','Fontsize',14)
% ylabel('Flow (veh/hr)','Fontsize',14)
% title('Demand and supply functions')

% Total length of the network
L = sum(L_section);

% Discrete vector of the network cells
X(1,:) = Delta_x/2:Delta_x:L-Delta_x/2;

% CFL condition
security_factor = 1.5; %should be greater or equal to 1
Delta_t = Delta_x / (security_factor*V_max);

% Initialization
n_T = length(0:Delta_t:T);
n_X = length(X(1,:));
rho = NaN(n_T,n_X);

rho(1,:) = rho_0(X(1,:));

lambda = Delta_t / Delta_x ;

% Loop in time
for i=1:n_T-1
    
    % Loop in space
    for j = 1
        % Index of the road to which the cell belongs
        index_j = X(2,j);
        index_j_plus = X(2,j+1);
        Supply = geometry(index_j).Supply ;
        Demand = geometry(index_j).Demand ;
        Supply_plus = geometry(index_j_plus).Supply;
        inflow = min( Demand_upstream(i*Delta_t), Supply(rho(i,j)) );
        outflow = min( Demand(rho(i,j)), Supply_plus(rho(i,j+1)) );
        rho(i+1,j) = rho(i,j) + lambda * (inflow-outflow);
        inflow = outflow;
    end
    
    for j= 2:n_X-1
        % Index of the road to which the cell belongs
        index_j = X(2,j);
        index_j_plus = X(2,j+1);
        Demand = geometry(index_j).Demand ;
        Supply_plus = geometry(index_j_plus).Supply ;
        outflow = min( Demand(rho(i,j)), Supply_plus(rho(i,j+1)) );
        rho(i+1,j) = rho(i,j) + lambda * (inflow-outflow);
        inflow = outflow;
    end
    
    for j = n_X
        % Index of the road to which the cell belongs
        index_j = X(2,j);
        Demand = geometry(index_j).Demand ;
        outflow = min( Demand(rho(i,j)), Supply_downstream(i*Delta_t) );
        rho(i+1,j) = rho(i,j) + lambda * (inflow-outflow);
    end
end

clear X

%{
% Plotting 
% -----------------------------------------
% --> could be a different script function:
% plot_density(L,T,Delta_x,Delta_t,rho)
% -----------------------------------------

figure
[X,Y]=meshgrid(0:Delta_t:T,Delta_x/2:Delta_x:L-Delta_x/2);
surf(X,Y,rho','EdgeColor','none');
view(2)
colormap(jet)
c = colorbar('Limits',[0 max(max(rho))]) ;
c.Label.String = 'Density (veh/km)'; % OPTIONNAL
xlabel('Time (hr)','Fontsize',14)
ylabel('Space (km)','Fontsize',14)
%}