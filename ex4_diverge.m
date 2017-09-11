% -----------------------------------------
% Minimal working example for Javier ORTIZ
% by Guillaume COSTESEQUE - September 2017
% 
% Comparison between two FIFO and non-FIFO
% diverge solvers
% -----------------------------------------
clear allvariables
clc

% ==========================================
% (I) Definition of the network
% ==========================================

%  -<: 1-2 diverge (3 links)

% (1) Specify the geometry characteristics

% Link 1
geometry(1).length=5; %in km
geometry(1).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(1).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(1).Vmax=90; %in km/hr

% Link 2
geometry(2).length=5; %in km
geometry(2).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(2).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(2).Vmax=90; %in km/hr

% Link 3
geometry(3).length=5; %in km
geometry(3).Demand=@(rho) (90.*rho).*(rho<=30) + (2700).*(rho>30);
geometry(3).Supply=@(rho) (2700).*(rho<=30) + (15.*(30-rho)+2700).*(rho>30);
geometry(3).Vmax=90; %in km/hr


nb_link = length(geometry) ;


% ==========================================
% (II) Initial and boundary conditions
% ==========================================

% (2) Enter the initial densities (constant on each link)
% /!\ here, it is a vector and not a function
Rho_0=[30 150 100] ; %in veh/km

% (3) Enter the upstream demand
Demand_upstream=@(t) 2000; %in veh/hr

% (4) Enter the downstream supply
Supply_downstream_2=@(t) 200; %in veh/hr

Supply_downstream_3=@(t) 800; %in veh/hr

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
    Density(link).rho_nonFIFO = Rho_0(link).*ones(1, ...
        length(Delta_x/2:Delta_x:geometry(link).length-Delta_x/2));
end

Rho(1,:) = Rho_0;
Rho_nonFIFO(1,:) = Rho_0;
i = 1;
for t=Delta_t:Delta_t:T
    
    % Run the FIFO diverge solver (links 1 to 3)
    rho_0 = Rho(i,1:3) ;
    A = 1/2 ;
    Q = diverge(geometry,A,rho_0) ;
    outflow_1 = Q(1);
    inflow_2 = Q(2);
    inflow_3 = Q(3);
    
    Demand_upstream_global=@(link,t) Demand_upstream(t).*(link==1) +...
        inflow_2.*(link==2) + inflow_3.*(link==3) ;
    
    Supply_downstream_global=@(link,t) outflow_1.*(link==1) +...
        Supply_downstream_2(t).*(link==2) +...
        Supply_downstream_3(t).*(link==3) ;
    
    
    % Run the non-FIFO diverge solver (links 1 to 3)
    rho_0 = Rho(i,1:3) ;
    A = 1/2 ;
    Q = diverge_nonFIFO(geometry,A,rho_0) ;
    outflow_1_nonFIFO = Q(1);
    inflow_2_nonFIFO = Q(2);
    inflow_3_nonFIFO = Q(3);
    
    Demand_upstream_global_nonFIFO=@(link,t) ...
        Demand_upstream(t).*(link==1) +...
        inflow_2_nonFIFO.*(link==2) + inflow_3_nonFIFO.*(link==3) ;
    
    Supply_downstream_global_nonFIFO=@(link,t) ...
        outflow_1_nonFIFO.*(link==1) +...
        Supply_downstream_2(t).*(link==2) +...
        Supply_downstream_3(t).*(link==3) ;

    
    % Run the Godunov solver for each link for a single time step
    for link = 1:nb_link
        
        % FIFO
        rho_0=@(x) Density(link).rho(end,floor(x/Delta_x)+1);
        
        [rho,~]=Godunov(geometry(link),rho_0,...
            @(t) Demand_upstream_global(link,t), ...
            @(t) Supply_downstream_global(link,t),...
            Delta_x,Delta_t) ;
        Rho(i+1,link) = (rho(2,end)).*(link==1) ...
            + (rho(2,1)).*(link==2 || link==3) ;   
        
        Density(link).rho = [Density(link).rho; rho(2,:)];
        
        % Non-FIFO
        rho_0_nonFIFO=@(x) Density(link).rho_nonFIFO(end,floor(x/Delta_x)+1);
        
        [rho_nonFIFO,~]=Godunov(geometry(link),rho_0_nonFIFO,...
            @(t) Demand_upstream_global_nonFIFO(link,t), ...
            @(t) Supply_downstream_global_nonFIFO(link,t),...
            Delta_x,Delta_t) ;
        Rho_nonFIFO(i+1,link) = (rho_nonFIFO(2,end)).*(link==1) ...
            + (rho_nonFIFO(2,1)).*(link==2 || link==3) ;   
        
        Density(link).rho_nonFIFO = [Density(link).rho_nonFIFO; ...
            rho_nonFIFO(2,:)];
        
        close all
    end
    
    i = i+1;
    
    clc
    fprintf('Computations: %0.3g%% done \n',(t/T)*100)
end

toc

% ==========================================
% (IV) Graphical representation
% ==========================================

pos=@(link) ([1,3]).*(link==1) + (2).*(link==2) + (4).*(link==3) ;

rho_max = -inf;
for link = 1:nb_link
    rho_max = max( rho_max , max ( ...
        max(max(Density(link).rho)), ...
        max(max(Density(link).rho_nonFIFO)) ) );
end

% FIFO
figure;

hAx = NaN(1,nb_link);
for link = 1:nb_link
    hAx(link) = subplot(2,2,pos(link)) ;
    [X,Y] = meshgrid(Delta_x/2:Delta_x:geometry(link).length-Delta_x/2, ...
        0:Delta_t:T) ;
    surf(X,Y,Density(link).rho,'EdgeColor','none')
    view(2)
    axis tight
    colormap(jet)
    caxis manual
    caxis([0 rho_max]);
end

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


% Non FIFO
figure;

hAx = NaN(1,nb_link);
for link = 1:nb_link
    hAx(link) = subplot(2,2,pos(link)) ;
    [X,Y] = meshgrid(Delta_x/2:Delta_x:geometry(link).length-Delta_x/2, ...
        0:Delta_t:T) ;
    surf(X,Y,Density(link).rho_nonFIFO,'EdgeColor','none')
    view(2)
    axis tight
    colormap(jet)
    caxis manual
    caxis([0 rho_max]);
end

h=colorbar;
top = 0.93;
bottom = .11;
width = .02;
heigth = .8150;
set(h, 'Position', [top bottom width heigth])

p1=get(hAx(1),'position');  % position of UL-most axes
p2=get(hAx(nb_link),'position');  % ditto for LR-most
hAxOuter=axes('position',[p1(1) p2(2) p2(1)+p2(3)-p1(1)  ...
    p1(2)+p1(4)-p2(2)], ...
    'color','none','visible','off'); % an outer axis for global p
hY=text(-0.08,0.5,'Time (hr)','rotation',90,'fontsize',14, ...
    'horizontalalignment','center','verticalalignment','bottom');
hX=text(0.5,-0.07,'Space (km)','rotation',0,'fontsize',14, ...
    'horizontalalignment','center','verticalalignment','top');