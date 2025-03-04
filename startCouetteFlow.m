clc; clear all; close all;

% This script sets up the starting flow conditions for Couette flow.
% The flow is initialized to be zero everywhere except at the top wall  
% where the velocity is set to 1. The flow is then allowed to evolve
% according to the diffusion equation. The flow is considered to be 
% steady when the L2 error between the current and previous time step
% is less than 1e-2. The flow is visualized every 50 time steps.
% written by Ebenezer Gnanamanickam for AE 521 
% This meant as a teaching tool and nothing beyond that.
% Students in AE521 are free to modify and use this code for their assigments
% This was commented with a lot of help from copilot

% set some default figure properties
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
y=linspace(0,1,50);
% wall-normal grid spacing
dy=y(2);

% kinematic viscosity
nu=.1;
% time step to satisfy stability
dt=0.5*dy^2/nu;

% set up the initial conditio
uInit=y.*0;
% set the top wall to be 1 at time t=0
uInit(end)=1;

% Loss function to calculate deviation from steady state. initialize to a large value
L2Loss=100;
% set up a counter to keep track of time steps
% loop until the flow is steady
cnt=1;
figure
while L2Loss>1e-2
    % increment counter
    cnt=cnt+1;

    % update the flow field using finite difference
    uFiniteDiff=uInit.*0;
    for j=2:1:length(y)-1
        uFiniteDiff(j)=(uInit(j+1)-2*uInit(j)+uInit(j-1))*...
            nu*dt/(dy^2)+uInit(j);
    end
    % set the top wall to be 1 - BC
    uFiniteDiff(end)=1; % boundary condition

    % calculate the L2 error - Note the steady state solution is u/U=y/H
    L2Loss=sqrt(sum((uFiniteDiff-y).^2));

    uInit=uFiniteDiff;

    if rem(cnt,50)==0
        clf;
        plot(uFiniteDiff,y,'-b',y,y,'-.r')
        title(['$t = $' num2str((cnt-1)*dt)])

        xlim([0 1])
        ylim([0 1])
        xlabel('$u/U$')
        ylabel('$y$')
        text(0.3,0.8,['$L_2$ error = ' num2str(L2Loss)],...
            'FontSize',fontSize)
        legend('Starting','Steady state')
        drawnow
    end
end

tEnd=(cnt-1)*dt;

clf;
plot(uFiniteDiff,y,'-b',y,y,'-.r')
title(['$t = $' num2str((cnt-1)*dt)])

xlim([0 1])
ylim([0 1])
xlabel('$u/U$')
ylabel('$y/H$')
text(0.3,0.8,['$L_2$ error = ' num2str(L2Loss)],...
    'FontSize',fontSize)
legend('Starting','Unsteady')
drawnow