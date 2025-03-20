
% coded by  Ebenezer Gnanamanickam for AE 521.
% This code calculates the Falkner-Skan boundary layer profile for a range of m values
% and plots the results. It also determines the value of m when the flow separates.
% This code is for teaching purposes only.

clc; clear all; close all;
fontSize=14;
set(0,'DefaultTextInterpreter','latex','DefaultAxesLineWidth',0.5,...
    'DefaultAxesFontSize',fontSize,'DefaultFigureInvertHardCopy',...
    'on','DefaultAxesFontName','Times','DefaultLineMarkerSize',6,...
    'DefaultLineLineWidth',1.5)
set(0,'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',fontSize,...
    'DefaultLegendOrientation','horizontal',...
    'DefaultLegendLocation','NorthOutside',...
    'DefaultLegendNumColumns',3,...
    'DefaultLegendNumColumnsMode','manual')

opts = odeset('RelTol',1e-2,'AbsTol',1e-4);

% determining exact value of m when the flow separates i.e., when f''(0)=0
% know from trial and error that it is about .091 but code will blow up for
% all m<-.0904. So calculate for a few values of m and extrapolate
mFit=(linspace(-.0904,-0.08,30));
initGuess=0.0;
etaMax=10;
y3_0=mFit.*0;


% plot the Falkner-Skan profiles for a range of m values as flow moves away from separation
figure
for i=1:1:length(mFit)
    [eta,y,err]=calFalknerSkan(mFit(i),etaMax,initGuess,opts);
    disp(err)

    initGuess=y(1,3);

    y3_0(i)=y(1,3);

    clf
    plot(y(:,2),eta,'r')
    ylabel('$\eta=y\sqrt{\frac{m+1}{2}\frac{U_e(x)}{\nu x}}$')
    xlabel('$\frac{u}{U_e(x)}=f^{\prime}$')
    ylim([0 4])
    drawnow
end

% use a second order fit to extrapolate to determine when f''(0)=0
PFit=polyfit(y3_0,mFit,2);

mZero=polyval(PFit,0);
disp(mZero)

figure
plot(mFit,y3_0,'xb',mZero,0,'xr',...
    polyval(PFit,[0 y3_0]),[0 y3_0],'-.k')

xlabel('$m$')
ylabel('$f^{\prime\prime}(0)$')

mPlot=[min(mFit) -.06 0 0.1 0.3 0.5 1 2 4];
etaMax=[10 9 8 7 6 5 4 3.5 3]; % need to adjust etaMax to avoid unrealistic solutions

lCol=colormap(['hot(' num2str(length(mPlot)*2) ')']);

initGuess=0;
figure
hold on;
for i=1:1:length(mPlot)
    [eta,y,err]=calFalknerSkan(mPlot(i),etaMax(i),initGuess,opts);
    disp(err)

    initGuess=y(1,3);

    y3_0(i)=y(1,3);

    plot(y(:,2),eta,'Color',lCol(i,:),...
        'DisplayName',['$m=$' num2str(mPlot(i))])
    drawnow
end

ylabel('$\eta=y\sqrt{\frac{m+1}{2}\frac{U_e(x)}{\nu x}}$')
xlabel('$\frac{u}{U_e(x)}=f^{\prime}$')
ylim([0 4])
box on
legend show

[eta,y,err]=calFalknerSkan(.3546,5,.5,opts);


figure
plot(y(:,2),eta,'-k')
ylabel('$\eta=y\sqrt{\frac{m+1}{2}\frac{U_e(x)}{\nu x}}$')
xlabel('$\frac{u}{U_e(x)}=f^{\prime}$')


return;

function [eta,y,err]=calFalknerSkan(m,etaMax,initGuess,opts)

clc
disp(m)
% need f'(inf)=1

err=10;
dEta=0.0001;
cnt=0;
while err>1e-6
    cnt=cnt+1;
    [eta,y]=ode45(@(eta,y) FalknerSkanEqn(eta,y,m),...
        [0 etaMax],[0; 0; initGuess],opts);
    [~,temp]=ode45(@(eta,y) FalknerSkanEqn(eta,y,m),...
        [0 etaMax],[0; 0; initGuess+dEta],opts);
    df=(temp(end,2)-y(end,2))/dEta;
    currGuess=initGuess;
    initGuess=initGuess+(1-y(end,2))/df;
    err=sqrt(((currGuess-initGuess)/currGuess)^2);
end

end

function fBL=FalknerSkanEqn(eta,y,m)

fBL=[y(2); y(3); -y(1)*y(3)-2*m/(m+1)*(1-(y(2))^2)];

end
