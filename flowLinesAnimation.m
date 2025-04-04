% Ebenezer P Gnanamanickam
% last updated February 2025
% inspired by the gif from the wikipedia article on streaklines,
% streamlines and pathlines

% This script generates an animation of streamlines, pathlines and streaklines
% for a flow field with a constant velocity in the x direction and a velocity
% in the y direction that increases with time. The flow is steady for n=0 and
% unsteady for larger values of n. The flow is visualized in a 2D domain with
% a grid of points. The pathline starts at a point (xi,yi) and the animation
% starts at time tSt and stops at time tStp. The animation is generated by
% plotting the streamlines, pathlines and streaklines at each time step.
% I shamelessly used copilot to generate some of the coments

% set default values for the plots. Big font as I am getting old
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

clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% power closer to zero the flow is more steady and larger the number more
% unsteady
n=0.5;

% start and stop time for the animation
% start time
tSt=0.5;
% stop time
tStp=2.0;

% initial position of the flow lines
xi=0.1;
yi=.1;

% number of time steps
tSteps=200;

% x and y grid points
nx=15;
ny=40;

% setup the grid
x=linspace(0,2,nx);
y=linspace(-4,4,ny);
[x,y]=meshgrid(x,y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% starting position for the streamlines
sx=x(:,1);
sy=y(:,1);

% time steps
t=linspace(0,tStp,tSteps);

% initialize the velocity fields
u=zeros(size(x,1),size(x,2),length(t));
v=u;

% calculate the velocity fields
for i=1:1:length(t)
    u(:,:,i)=x.*0+1;
    v(:,:,i)=x.*t(i).^(n);
end

% initialize the pathline
xPath=t.*0;
yPath=xPath;
% initialize the streakline
xStreak=xPath+xi;
yStreak=xPath+yi;

[~,iSt]=min(abs(t-tSt));

% Initialize the starting position of the pathline
xPath(iSt)=xi;
yPath(iSt)=yi;

% Create figure window with specified position and size
figure('pos',[1843 151 771 526])
for i=iSt:1:length(t)
    clf; % Clear current figure
    % Get current velocity field
    uCurr=u(:,:,i);
    vCurr=v(:,:,i);
    
    hold on;
    % Plot streamlines starting from all points along left boundary
    h1=streamline(x,y,uCurr,vCurr,sx,sy);
    set(h1,'Color','k','LineStyle','-.','LineWidth',.75)
    % Plot streamline starting from initial point
    hs=streamline(x,y,uCurr,vCurr,xi,yi);
    set(hs,'Color','k','LineStyle','-.','LineWidth',2)
    % Plot velocity vectors
    hq=quiver(x,y,uCurr,vCurr,'g');
    
    % Calculate pathline position using velocity interpolation
    if i~=iSt
        uPath=interp2(x,y,uCurr,xCurr,yCurr);
        vPath=interp2(x,y,vCurr,xCurr,yCurr);
        xPath(i)=xCurr+uPath*(t(i)-t(i-1));
        yPath(i)=yCurr+vPath*(t(i)-t(i-1));
    end
    
    % Update streakline starting position
    xStreak(i)=xPath(i);
    yStreak(i)=yPath(i);
    
    % Calculate streakline positions by tracking particles released earlier
    if i>iSt+1
        for j=i-1:-1:iSt
            uSt=interp2(x,y,uCurr,xStreak(j-1),yStreak(j-1));
            vSt=interp2(x,y,vCurr,xStreak(j-1),yStreak(j-1));
            xStreak(j)=xStreak(j-1)+uSt*(t(j)-t(j-1));
            yStreak(j)=yStreak(j-1)+vSt*(t(j)-t(j-1));
        end
    end
    
    % Plot streakline
    hSt=plot(xStreak(iSt:i),yStreak(iSt:i),'-b','LineWidth',2);
    
    % Plot pathline and current particle position
    h2=plot(xPath(iSt:i),yPath(iSt:i),'--r','LineWidth',2);
    plot(xPath(i),yPath(i),'or')
    
    % Update current position
    xCurr=xPath(i);
    yCurr=yPath(i);
    
    % Set plot properties
    axis([-0 1.8  -0 1.8])
    legend([h1(1),h2(1),hq(1),hSt(1)],'Streamlines','Pathline',...
        'Velocity vectors','Streaklines','Location',...
        'NorthOutside','Interpreter','latex','NumColumns',2)
    xlabel('$x$')
    ylabel('$y$')
    box('on')
    drawnow
    hold off;
end
