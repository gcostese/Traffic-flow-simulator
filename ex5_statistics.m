% -----------------------------------------
% Minimal working example
% by Guillaume COSTESEQUE - September 2017
%
% Compute a posteriori some statistics on 
% the traffic state
% -----------------------------------------
clear allvariables
clc

% ==========================================
% (I) Definition of the network
% ==========================================

geometry(1).length = 10 ;
geometry(1).Vmax = 90 ;
geometry(1).rho_crit = 30 ;
geometry(1).Demand = @(k) min( geometry(1).Vmax * geometry(1).rho_crit ,...
    geometry(1).Vmax.*k );
geometry(1).Supply = @(k) min( geometry(1).Vmax * geometry(1).rho_crit ,...
    geometry(1).Vmax * geometry(1).rho_crit + ...
    15.*(geometry(1).rho_crit - k) );
geometry(1).Flux = @(k) min( geometry(1).Vmax.*k ,...
    geometry(1).Vmax * geometry(1).rho_crit +...
    15.*(geometry(1).rho_crit - k) );

geometry(2).length = 5 ;
geometry(2).Vmax = 70 ;
geometry(2).rho_crit = 20 ;
geometry(2).Demand = @(k) min( geometry(2).Vmax * geometry(2).rho_crit ,...
    geometry(2).Vmax.*k );
geometry(2).Supply = @(k) min( geometry(2).Vmax * geometry(2).rho_crit ,...
    geometry(2).Vmax * geometry(2).rho_crit + ...
    15.*(geometry(2).rho_crit - k) );

Delta_x = 0.25 ;

density(1).rho = 100.* rand(20,floor(geometry(1).length/Delta_x)+1) ;
density(2).rho = 80.* rand(20,floor(geometry(2).length/Delta_x)+1) ;

% ==========================================
% (II) Computations of some traffic indicators
% ==========================================

[TTT,TD,QL,FC] = statistics(geometry,density,Delta_x) ;