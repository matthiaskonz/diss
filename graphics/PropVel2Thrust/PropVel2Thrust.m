clear;

addpath(genpath('/home/matthias/MatlabFcns'));

%% load
load experimentThrustAngVel;

FR = erg(:,1);
w1 = erg(:,2);
w2 = erg(:,3);
w3 = erg(:,4);

%% least squares approximation
kappaF = [w1.^2; w2.^2; w3.^2] \ [FR; FR; FR]

%% figure
fig = figure(1); clf;
sp = subplot(1,1,1);
line(w1, FR, 'Color', 'r', 'LineStyle', 'none', 'Marker', 's', 'Markersize', 8);
%line(w2, FR, 'Color', 'g', 'LineStyle', 'none', 'Marker', '.', 'Markersize', 8);
line(w3, FR, 'Color', 'b', 'LineStyle', 'none', 'Marker', 'o', 'Markersize', 6);

wR = linspace(1, 800, 100);
line(wR, kappaF*wR.^2 , 'Color', 'k');
grid on;

leg = legend( ...
  'clockwise propeller', ...
  'counterclockwise propeller', ...
  'model', ...
'Location', 'NorthWest');

xlabel('angular velocity in RAD/s');
ylabel('thrust in N');
xlim([0 750]);
ylim([0 8]);

%% save figure
% set(fig1, 'position', [100 100 600 400]);
% export_fig('PropVel2Thrust.pdf', fig1, '-transparent');
%%
saveFigure = 1;
figSize = [12 7];

if (saveFigure > 0)
  
  sp.UserData.LatexYLabel = '$\PropForce$ in $\unit{N}$';
  sp.UserData.LatexXLabel = '$\PropVel$ in $\sfrac{\unit{RAD}}{\unit{s}}$';
% 	spProp(3).UserData.LatexYLabel = '$\BLDCCurr$ in $\unit{A}$';
%   spProp(4).UserData.LatexYLabel = '$\hat{\tau}_{\mathsf{MB}}$ in $\unit{Nm}$';
%   spProp(4).UserData.LatexXLabel = '$t$ in $\unit{s}$';
%   
%  leg.UserData.LatexString = {'clockwise propeller', 'counterclockwise propeller', 'LQ ident'};
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'PropVel2Thrust', 'pdf', 0 );
  
end
