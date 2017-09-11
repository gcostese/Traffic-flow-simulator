function [TTT,TD,QL,FC] = statistics(geometry,density,Delta_x)
% =========================================================================
% This function is used to compute some traffic indicators for a given set
% of links on the network
% -------------------------------------------------------------------------
% Inputs:
% (1) key network functions encoded in the *geometry* structure:
%     (a) length of each link
%     (b) maximal speed for each link
%     (c) critical density for each link
%     (d) the flow-density function
% (2) the *density* structure containing the matrix of all discrete
%     densities at different time steps
% (3) the value of the space step *Delta_x*
% -------------------------------------------------------------------------
% Outputs:
% (1) the Total Travel Time *TTT*
% (2) the Total Delay *TD*
% (3) the queue length *QL*
% =========================================================================

%Number of links on the network
nb_link = length(geometry) ;
L = 0;
weight=NaN(1,nb_link);
for link=1:nb_link
    %total length of the network
    L = L + geometry(link).length;
    %weight of the link length on the total network length
    weight(link) = geometry(link).length;
    
    %definition of the flow-density function
    if ~isfield(geometry,'Flux') || isempty(geometry(link).Flux)
        geometry(link).Flux = @(rho) min( geometry(link).Demand(rho), ...
            geometry(link).Supply(rho) ) ;
    end
end
weight=weight/L;

%Number of time steps to be considered
nb_time_step = +Inf;
for link=1:nb_link
    nb_time_step = min ( nb_time_step, length( density(link).rho(:,1) ) ) ;
end

% ==========================================
% (I) Total Travel Time
%     and Average Speed
% ==========================================

%Initialization of the vector
TTT = zeros(1,nb_time_step) ;
average_speed = NaN(1,nb_time_step) ;
for t=1:nb_time_step
    speed_link = NaN(1,nb_link) ;
    for link=1:nb_link
        K = density(link).rho(t,:) ; %vector of all the densities on the link at time step *t*
        f = geometry(link).Flux ;    %flux function for the considered link
        TTT(t) = TTT(t) + Delta_x * sum( K ./ f(K) ) ; %the total travel time is given by the sum of all the travel times on each link
        
        speed_link(link) = mean( f(K)./K ) ; %mean speed on the link
    end
    average_speed(t) = mean( weight.*speed_link ) ; %mean speed on the network at time step *t*
end

% ==========================================
% (II) Total Delay
% ==========================================

% Normal Travel Time
NTT=0;
for link=1:nb_link
    NTT = NTT + geometry(link).length / geometry(link).Vmax ;
end
NTT = NTT.*ones(1,nb_time_step) ;
TD = TTT - NTT;

% ==========================================
% (III) Total queue length
% ==========================================
QL = NaN(1,nb_time_step) ;
for t=1:nb_time_step
    QL(t) = 0;
    if TD(t) > 0
        for link=1:nb_link
            QL(t) = QL(t) + length( find( ...
                density(link).rho(t,:) > geometry(link).rho_crit) ) ...
                * Delta_x ;
        end
    end
end

% Fraction of the network that is congested
FC = 100*QL/L ;

% ==========================================
% (IV) Graphical representation
% ==========================================

T = 1:nb_time_step ;

figure

%Total Travel Time
subplot(2,2,1)
plot(T,TTT)
ylabel('Total Travel Time (hr)')
axis tight

%Total Delay
subplot(2,2,2)
plot(T,TD)
ylabel('Total Delay (hr)')
axis tight

%Queue length
subplot(2,2,3)
plot(T,QL)
ylabel('Queue length (km)')
axis tight

%Fraction of congested network
subplot(2,2,4)
plot(T,FC)
ylabel('Portion of congested network (%)')
axis tight

%{ 
subplot(2,2,[3 4])
[hAx,~,~] = plotyy(T,QL,T,FC) ;
ylabel(hAx(1),'Queue length (km)') % left y-axis 
ylabel(hAx(2),'Portion of congested network (%)') % right y-axis
    %{
    %For more recent MATLAB versions
    yyaxis left
    plot(QL)
    ylabel('Queue length (km)')
    yyaxis right
    plot(FC)
    ylabel('Portion of congested network (%)')
    %}
%}

%Average Speed
figure
hold on
plot(average_speed)
plot(L./TTT)
ylabel('Speed (km/hr)')
legend('Average','Length / TTT')
axis tight
hold off