% -----------------------------------------
% Minimal working example for Javier ORTIZ
% by Guillaume COSTESEQUE - August 2017
%
% Simulate traffic densities on a simple
% network with one diverge and one merge
% -----------------------------------------
clear allvariables
clc

% ==========================================
% (I) Definition of the network
% ==========================================

% 6 sections: -<>-
%  -<: 1-2 diverge (3 links) and >-: 2-1 merge (3 links)

% (1) Specify the geometry characteristics

% Link 1
geometry(1).length=5; %in km
geometry(1).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(1).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(1).Vmax=90; %in km/hr
geometry(1).rho_crit=30; %in veh/km

% Link 2
geometry(2).length=5; %in km
geometry(2).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(2).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(2).Vmax=90; %in km/hr
geometry(2).rho_crit=30; %in veh/km

% Link 3
geometry(3).length=5; %in km
geometry(3).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(3).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(3).Vmax=90; %in km/hr
geometry(3).rho_crit=30; %in veh/km

% Link 4
geometry(4).length=5; %in km
geometry(4).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(4).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(4).Vmax=90; %in km/hr
geometry(4).rho_crit=30; %in veh/km

% Link 5
geometry(5).length=5; %in km
geometry(5).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(5).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(5).Vmax=90; %in km/hr
geometry(5).rho_crit=30; %in veh/km

% Link 6
geometry(6).length=5; %in km
geometry(6).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(6).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(6).Vmax=90; %in km/hr
geometry(6).rho_crit=30; %in veh/km


nb_link = length(geometry) ;


% ==========================================
% (II) Initial and boundary conditions
% ==========================================

% (2) Enter the initial densities (constant on each link)
% /!\ here, it is a vector and not a function
Rho_0=[125 100 50 20 80 120] ; %in veh/km

% (3) Enter the upstream demand
Demand_upstream=@(t) 1000; %in veh/hr

% (4) Enter the downstream supply
Supply_downstream=@(t) 1500; %in veh/hr


% ==========================================
% (III) Numerical scheme
% ==========================================

% (5) Specify the discrete step in space and the time horizon
Delta_x = 0.2; %in km
T = 0.6;       %in hour

% (6) CFL condition
V_max = -inf;
for i = 1:nb_link
    V_max = max(V_max, geometry(i).Vmax) ;
end
k = 1.5; %security factor
Delta_t = Delta_x / (k*V_max) ;


% (7) Time loop
tic

Density = struct;
for link = 1:nb_link
    Density(link).rho = Rho_0(link).*ones(1, ...
        length(Delta_x/2:Delta_x:geometry(link).length-Delta_x/2));
end

Rho(1,:) = Rho_0;
i = 1;
for t=Delta_t:Delta_t:T
    
    % Run the diverge solver (links 1 to 3)
    rho_0 = Rho(i,1:3) ;
    A = 1/2 ;
    Q = diverge(geometry,A,rho_0) ;
    outflow_1 = Q(1);
    inflow_2 = Q(2);
    inflow_3 = Q(3);
    
    % Run the merge solver (links 4 to 6)
    rho_0 = Rho(i,4:6) ;
    P = 1/2 ;
    Q = merge(geometry,P,rho_0) ;
    outflow_4 = Q(1);
    outflow_5 = Q(2);
    inflow_6 = Q(3);
    
    demand_4 = geometry(2).Supply(Density(2).rho(i,end));
    demand_5 = geometry(3).Supply(Density(3).rho(i,end));
    supply_2 = geometry(4).Supply(Density(4).rho(i,1)) ;
    supply_3 = geometry(5).Supply(Density(5).rho(i,1)) ;
    
    Demand_upstream_global=@(link,t) Demand_upstream(t).*(link==1) +...
        inflow_2.*(link==2) + inflow_3.*(link==3) +...
        demand_4.*(link==4) + demand_5.*(link==5) +...
        inflow_6.*(link==6);
    
    Supply_downstream_global=@(link,t) outflow_1.*(link==1) +...
        supply_2.*(link==2) + supply_3.*(link==3) +...
        outflow_4.*(link==4) + outflow_5.*(link==5) +...
        Supply_downstream(t).*(link==6) ;

    % Run the Godunov solver for each link for a single time step
    for link = 1:nb_link
        rho_0=@(x) Density(link).rho(end,floor(x/Delta_x)+1);
        
        [rho,~]=Godunov(geometry(link),rho_0,...
            @(t) Demand_upstream_global(link,t), ...
            @(t) Supply_downstream_global(link,t),...
            Delta_x,Delta_t) ;
        Rho(i+1,link) = (rho(2,end)).*(link==1 || link==4 || link ==5) ...
            + (rho(2,1)).*(link==2 || link==3 || link ==6) ;   
        
        Density(link).rho = [Density(link).rho; rho(2,:)];
        
        close all
    end
    
    i = i+1;
    
    clc
    fprintf('Computations: %0.3g%% done \n',(t/T)*100)
end

toc

clearvars -except geometry nb_link Delta_x Delta_t T Rho Density

% ==========================================
% (IV) Graphical representation
% ==========================================

figure;

rho_max = -inf;
for link = 1:nb_link
    rho_max = max(rho_max , max(max(Density(link).rho)) );
end

pos=@(link) ([1,5]).*(link==1) + (2).*(link==2) + (6).*(link==3) + ...
  (3).*(link==4) + (7).*(link==5) + ([4,8]).*(link==6) ;

hAx = NaN(1,nb_link);
for link = 1:nb_link
    hAx(link) = subplot(2,4,pos(link)) ;
    [X,Y] = meshgrid(Delta_x/2:Delta_x:geometry(link).length-Delta_x/2, ...
        0:Delta_t:T) ;
    surf(X,Y,Density(link).rho,'EdgeColor','none')
    view(2)
    axis tight
    colormap(jet)
    caxis manual
    caxis([0 rho_max]);
end

clear X Y

h=colorbar;
top = 0.93;
bottom = .11;
width = .02;
heigth = .8150;
set(h, 'Position', [top bottom width heigth])

p1=get(hAx(1),'position');  % position of UL-most axes
p2=get(hAx(nb_link),'position');  % ditto for LR-most
axes('position',[p1(1) p2(2) p2(1)+p2(3)-p1(1)  ...
    p1(2)+p1(4)-p2(2)], ...
    'color','none','visible','off'); % an outer axis for global p
text(-0.08,0.5,'Time (hr)','rotation',90,'fontsize',14, ...
    'horizontalalignment','center','verticalalignment','bottom');
text(0.5,-0.07,'Space (km)','rotation',0,'fontsize',14, ...
    'horizontalalignment','center','verticalalignment','top');

clear p1 p2 h top bottom width heigth

% ==========================================
% (V) Computations of some traffic indicators
% ==========================================

[TTT,TD,QL,FC] = statistics(geometry,Density,Delta_x) ;