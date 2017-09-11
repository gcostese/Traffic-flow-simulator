function plot_density(L,T,Delta_x,Delta_t,rho)
% =========================================================================
% This function is used to plot a colormap of densities 
% on a 2D time-space mesh
% -------------------------------------------------------------------------
% Inputs:
% (1) link length *L*
% (2) time horizon *T*
% (3) discrete spatial step *Delta_x*
% (4) discrete time step *Delta_t*
% (5) matrix of all the densities *rho*
% -------------------------------------------------------------------------
% Output:
% (1) figure with a colormap of the densities depending on 
%     time on x-axis and space on y-axis
% =========================================================================

figure
[X,Y]=meshgrid(0:Delta_t:T,Delta_x/2:Delta_x:L-Delta_x/2);
surf(X,Y,rho','EdgeColor','none');
view(2)
axis tight
colormap(jet)
c = colorbar('Limits',[0 max(max(rho))]) ;
c.Label.String = 'Density (veh/km)'; % OPTIONNAL
xlabel('Time (hr)','Fontsize',14)
ylabel('Space (km)','Fontsize',14)