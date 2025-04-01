clc; clear all; close all;

% Set default properties for the plots
fontSize=16;
set(0,'DefaultTextInterpreter','latex','DefaultAxesLineWidth',0.5,...
    'DefaultAxesFontSize',fontSize,'DefaultFigureInvertHardCopy',...
    'on','DefaultAxesFontName','Times','DefaultLineMarkerSize',6,...
    'DefaultLineLineWidth',1.5)
set(0,'DefaultLegendInterpreter','latex',...
    'DefaultLegendFontSize',fontSize,...
    'DefaultLegendOrientation','horizontal',...
    'DefaultLegendLocation','NorthOutside',...
    'DefaultLegendNumColumns',2,...
    'DefaultLegendNumColumnsMode','manual')

% Prandtl number
Pr=[0.3 0.5 0.7 1.5 2 5 10 100];

% initialize figure
figure
lCol=colormap(['hot(' num2str(length(Pr)*2) ')']);
hold on;
% loop over Prandtl number
for i=1:1:length(Pr)
    [eta,uByUInf,theta]=calBoundaryLayerZPG(Pr(i));
    
    % plot the velocity profile (Pr=1)
    if i==1
        plot(eta,uByUInf,'-b',...
            'DisplayName','$u/U_\infty,~Pr=1$')
    end
    
    % plot the temperature profile
    plot(eta,theta,'Color',lCol(i,:),'LineStyle','--',...
        'DisplayName',['$Pr=$' num2str(Pr(i))])
end

xlabel('$\eta=y\sqrt{\frac{U_\infty}{\nu x}}$')
ylabel('{$u/U_\infty$},~{${(T-T_w)}/{(T_\infty-T_w)}$}')
xlim([0 8])
box on
legend show

% asses the validity of the empirical correlations
% Prandtl numbers
Pr=[0.005 0.01 .05 0.1 0.3 0.5 0.7 1 1.25 1.5 1.75 2 3 5 10 30 70 100 150];
Nu_xSqrtRe_x=Pr.*0;

for i=1:1:length(Pr)
    [~,~,~,Nu_xSqrtRe_x(i)]=calBoundaryLayerZPG(Pr(i));
end

PrPlot1=linspace(0.005,0.5,100);
PrPlot2=linspace(0.005,150,1000);
PrPlot3=linspace(.01,150,100);

figure
semilogy(Nu_xSqrtRe_x,Pr,'og',PrPlot1.^(1/2)*0.565,PrPlot1,'--k',...
    PrPlot2.^(1/3)*0.332,PrPlot2,'--r',...
    PrPlot3.^(1/3)*0.339,PrPlot3,'--b')
xlabel('$Nu_xRe_x^{-1/2}$')
ylabel('$Pr$')
legend('Similarity solution',...
    '$Nu_x=0.565 Pr^{1/2} Re_x^{1/2}$',...
    '$Nu_x=0.332 Pr^{1/3} Re_x^{1/2}$',...
    '$Nu_x=0.339 Pr^{1/3} Re_x^{1/2}$')

function [eta,uByUInf,theta,Nu_xSqrtRe_x]=calBoundaryLayerZPG(Pr)

opts=odeset('RelTol',1e-2,'AbsTol',1e-4);

etaMax=8;
cntMax=100;
nPts=200;

% need f'(inf)=1
initGuess=0.0;
err=100;
dEta=0.0001;
cnt=0;
while err>1e-5
    cnt=cnt+1;
    [eta,y]=ode45(@BlasiusEqn,linspace(0,etaMax,200),...
        [0; 0; initGuess]);
    [~,temp]=ode45(@BlasiusEqn,linspace(0,etaMax,200),...
        [0; 0; initGuess+dEta]);
    df=(temp(end,2)-y(end,2))/dEta;
    currGuess=initGuess;
    initGuess=initGuess+(1-y(end,2))/df;
    err=sqrt(((currGuess-initGuess)/currGuess)^2);
    if cnt>cntMax
        error=1e-8;
    end
end

theta=eta*0;
temp1=cumtrapz(eta,y(:,1));

for i=2:1:length(eta)
    theta(i)=trapz(eta(1:i),exp(-Pr/2.*temp1(1:i)))/...
        trapz(eta,exp(-Pr/2.*temp1));
end

uByUInf=y(:,2);
Nu_xSqrtRe_x=1./trapz(eta,exp(-Pr/2.*temp1));

end

function fBL=BlasiusEqn(eta,y)

fBL=[y(2); y(3); -y(1)*y(3)/2];

end