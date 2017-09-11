% ------------------------------------
% Minimal working example
% by Guillaume COSTESEQUE - July 2017
% ------------------------------------
clear allvariables
clc

% (1) Specify the geometry characteristics
geometry.length(1)=5; %in km
geometry.Demand(1)=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry.Supply(1)=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry.Vmax(1)=90; %in km/hr

% --------------------- OPTIONNAL -----------------------
% Graphical representation of the demand/supply functions
figure
hold on
ezplot(geometry.Demand,[0 200])
ezplot(geometry.Supply,[0 200])
hold off
% -------------------------------------------------------

% (2) Enter the initial densities
rho_0=@(x) 20.*(x<=0.5) + 100.*(x>0.5) ; %in veh/km

% (3) Enter the upstream demand
Demand_upstream=@(t) 1200; %in veh/hr

% (4) Enter the downstream supply
Supply_downstream=@(t) 2000; %in veh/hr

% (5) Specify the discrete step in space and the time horizon
Delta_x = 0.2; %in km
T = 0.5;        %in hour

% Run the Godunov solver
[rho,Delta_t]=Godunov(geometry,rho_0,...
    Demand_upstream,Supply_downstream,...
    Delta_x,T) ;
