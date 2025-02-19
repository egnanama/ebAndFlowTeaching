
% This is just a simple plotting utility to visualize Coutte flow. Both the velocity and temperature profiles are plotted
clc; close all; clear all;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters
% mu: dynamic viscosity (Pa.s)
mu=[1.8e-3 1e-3 0.29];
% k: thermal conductivity (W/m.K)
k=[0.26 0.6 0.145];
% Br: Brinkman number
Br=[.007 .017 20];
% upper wall temperature (K)
Tu=350;
% lower wall temperature (K)
Tl=300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temperature difference
dT=Tu-Tl;

% y: normalized distance from the lower wall, -1 is the bottom wall and 1 is the top wall
y=linspace(-1,1,100);

% non-dimensionalized velocity profile
u=1/2.*(y+1);

% plot the velocity profile
figure
plot(u,y,'r')
xlabel('$u/U$')
ylabel('$y/h$')
xlim([0 1.1]);
ylim([-1 1]);

% non-dimensionalized temperature profile
T=zeros(length(Br),length(y));
for i=1:1:length(Br)
    T(i,:)=(Tu+Tl)/2+(Tu-Tl)/2.*y+Br(i)/8*(Tu-Tl).*(1-y.^2);
end

% plot the temperature profile
figure
plot(T(1,:)./dT,y,'k',T(2,:)./dT,y,'-.r',T(3,:)./dT,y,'--b')
legend(['$Br=$' num2str(Br(1))],['$Br=$' num2str(Br(2))],...
    ['$Br=$' num2str(Br(3))])
xlabel('$T/(T_u-T_l)$')
ylabel('$y/h$')

% Nusselt number
NusseltNo=1+Br/2
NusseltNo=1-Br/2

