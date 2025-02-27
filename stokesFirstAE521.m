clc; clear all; close all;

% This code compares a finite difference solution of Stoke's first problem
% with the analytical solution. 
% The discretization is forward differnce (first order) in time and
% central difference (second order) in space
% written by Ebenezer P. Gnanamanickam for AE 521 Spring
% This is a simple 1D diffusion problem. You will need to modofy this code for your homework
% Some Questions to ask?
% Physics 
% what happens when the viscosity is increased or decreased?
% Do you see how this problem has no length or time scale i.e., its self similar?
% Where does all the vorticity come from or rather at what instant?
% Numerics
% We need one initial condition and two boundary conditions - what type of PDE is this?
% Change the time step and what happens to the error?
% What happens to the error if you change the grid spacing?
% What happens to the error if you change the viscosity?
% Somethings to try
% can you flip the problem i.e., the wall and fluid id moving at U and the wall is suddenly brought to rest at t=0?
% can you make this a Coutte flow starting problem?


fontSize=14;
set(0,'DefaultTextInterpreter','latex','DefaultAxesLineWidth',0.5,...
    'DefaultAxesFontSize',fontSize,'DefaultFigureInvertHardCopy',...
    'on','DefaultAxesFontName','Times','DefaultLineMarkerSize',6,...
    'DefaultLineLineWidth',1.5)
set(0,'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',fontSize,...
    'DefaultLegendOrientation','horizontal',...
    'DefaultLegendLocation','NorthOutside',...
    'DefaultLegendNumColumns',2)

%% Inputs
% wall normal grid - make sure that max(y) is large enough such that the
% infinity boundary condition is maintained.
y=linspace(0,0.5,150);
% kinematic viscosity
nu=1;
% wall boundary condition
uWall=10; % wall boundary condition 
% max time
tMax=0.04;

%% Rest of the code
dy=y(2);
% calculate the time step size based on the criteria that the diffusion
% number d<0.5
dt=0.5*dy^2/nu;
t=0:dt:tMax;

% initialize the velocity field 
uInit=y.*0;
% apply the boundary condition
uInit(1)=uWall; 

figure('Position',[440   105   973   735])
for i=2:1:length(t)
    % analytical u and omega
    uAnal=uWall*(1-erf(y./(2*sqrt(nu*t(i))))); 
    omegaAnal=uWall/sqrt(nu*pi*t(i)).*exp(-(y.^2./(4*nu*t(i))));
    
    % use omega at t=dt to non-dimensionalize plots
    if i==2
        omegaT0=omegaAnal(1);
    end
    
    % implement the FDE solution
    uFiniteDiff=uInit.*0;
    for j=2:1:length(y)-1
        uFiniteDiff(j)=(uInit(j+1)-2*uInit(j)+uInit(j-1))*...
            nu*dt/(dy^2)+uInit(j);
    end
    % enforce wall BC
    uFiniteDiff(1)=uWall; 
    
    
    % reinitialize u for next time step
    uInit=uFiniteDiff;

    % plot every 10 time steps
    if rem(i,10)==0
        clf;
        subplot(1,2,1)
        plot(uAnal(1:5:end)/uWall,y(1:5:end),'xr',...
            uFiniteDiff/uWall,y,'-.b')
        xlim([0 1])
        ylim([0 0.5])
        xlabel('$u/U$')
        ylabel('$y$')
        legend('Exact solution','Finite difference')
        
        subplot(1,2,2)
        plot(omegaAnal/omegaT0,y,'-.r')
        xlim([-0.7 0.7])
        ylim([0 0.5])
        xlabel('$\omega_z/\omega_{z(t=\delta t)}$')
        ylabel('$y$')

        drawnow
    end
end


