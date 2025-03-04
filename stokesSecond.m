clc; clear all; close all;

% this is a simple finite difference solution to Stoke's second problem
% written by Ebenezer P. Gnanamanickam for AE 521 Spring 2025
% Depending on the year you may need to modify the code.

fontSize=16;
set(0,'DefaultTextInterpreter','latex','DefaultAxesLineWidth',0.5,...
    'DefaultAxesFontSize',fontSize,'DefaultFigureInvertHardCopy',...
    'on','DefaultAxesFontName','Times','DefaultLineMarkerSize',6,...
    'DefaultLineLineWidth',1.5)
set(0,'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',fontSize,...
    'DefaultLegendOrientation','horizontal',...
    'DefaultLegendLocation','NorthOutside',...
    'DefaultLegendNumColumns',2)

% set up the wall-normal grid
y=linspace(0,1,100);
% wall-normal grid spacing
dy=y(2);

% frequncy of oscillation
omega=2*pi*100;

% kinematic viscosity
nu=10;

% time step size based on the criteria that the diffusion
% number d<=0.5
dt=0.5*dy^2/nu;
% time vector
t=0:dt:0.1;

%% initialize the velocity field
uInit=y.*0;
% set the bottom wall to be 1 at time t=0
uInit(1)=1;
uInit(end)=0;

% set the counter to plot at specific time steps
plotStep=round(linspace(2,length(t),3000));
plotCounter=1;


figure
for i=2:1:length(t)    
    uFiniteDiff=y.*0;
    % advance the solution in time using the finite difference method
    for j=2:1:length(y)-1
        uFiniteDiff(j)=(uInit(j+1)-2*uInit(j)+uInit(j-1))*...
            nu*dt/(dy^2)+uInit(j);
    end
    % enforce wall BC
    uFiniteDiff(end)=0; % boundary condition
    uFiniteDiff(1)=cos(omega*t(i)); % boundary condition
    
    % reinitialize u for next time step
    uInit=uFiniteDiff;
    
    % plot the solution
    if i==plotStep(plotCounter)
        plotCounter=plotCounter+1;
        plot(uFiniteDiff,y,'-b')
    end
    
    xlim([-1 1])
    ylim([0 1])
    xlabel('$u/U$')
    ylabel('$y$')
    drawnow
end


