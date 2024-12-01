clear;
filename = 'ServoTestStand4';
AgentNr = 1;

saveFigure = 1;
figSize = [15.9 17];

tW = 52.4;
TW = 2;

%% param
Ts = 5e-3;
JP = 3.7e-5;
kappaF = 1.467e-05;
iBLDC2A = 0.0177;

% input param
t0 = 1;
tEnd = 1000;

%% read measuremnt
load(filename);
meas = meas{AgentNr};
kk = find((meas.vSerial(1,:) > t0) & (meas.vSerial(1,:) <= tEnd));

tM = meas.vSerial(1,kk);
t = [tM(5):Ts:tM(end-5)];

FR = interp1(tM, meas.vSerial(8,kk), t);
FD = interp1(tM, meas.vSerial(9,kk), t);
FO = interp1(tM, meas.vSerial(10,kk), t+Ts);
iM = interp1(tM, meas.vSerial(3,kk), t);
p3O = interp1(tM, meas.vSerial(5,kk), t);
w  = interp1(tM, meas.vSerial(7,kk), t+1.5*Ts);
F = kappaF * w.^2;

aR = interp1(tM, meas.vSerial(18,kk), t); 
aU = interp1(tM, meas.vSerial(19,kk), t); 
aS = interp1(tM, meas.vSerial(20,kk), t); 
a  = interp1(tM, meas.vSerial(21,kk), t+Ts); % measured servo angle

t = t - tW;


%%
fig = figure(1); clf;
spProp(1) = subplot(5,1,1:2); grid on; hold on;
l11 = line(t, F, 'Color', 0.7*[1 1 1], 'Linewidth', 3);
l12 = stairs(t, FR, 'Color', 'k', 'Linewidth', 1);
l13 = stairs(t, FD, 'Color', 'b', 'Linewidth', 1);
l14 = stairs(t, FO, 'Color', 'r', 'Linewidth', 1);
ylim([0 7]);
ylabel('F in N');
spProp(2) = subplot(5,1,3); grid on; hold on;
stairs(t, F-FO, 'Color', 'r', 'Linewidth', 1);
stairs(t, FD-FO, 'Color', 'b', 'Linewidth', 1);
%stairs(t, FO-FR, 'Color', 'k', 'Linewidth', 1);
ylim(0.3*[-1 1]);
ylabel('error in N');
spProp(3) = subplot(5,1,4); grid on; hold on;
stairs(t, iBLDC2A*iM, 'Color', 'b', 'Linewidth', 1);
ylabel('iM in A');
ylim(iBLDC2A*[-400 800]);
%spProp(3).YTick = iBLDC2A*[-400 0 400 800];
spProp(4) = subplot(5,1,5); grid on; hold on;
stairs(t, JP*p3O, 'Color', 'b', 'Linewidth', 1);
%ylim([0 300]);
ylabel('bias torque');
xlabel('time in s');
linkaxes(spProp, 'x');

xlim([0 TW]);

axes(spProp(1));
l1 = legend([l11 l14 l12 l13], 'measured', 'observed', 'desired', 'filtered des.', 'location', 'southwest');

axes(spProp(2));
l2 = legend('eObs', 'eCtrl', 'location', 'southwest');

%%
if (saveFigure > 0)
  
  spProp(1).UserData.LatexYLabel = '$\PropForce$ in $\unit{N}$';
  spProp(2).UserData.LatexYLabel = 'error in $\unit{N}$';
	spProp(3).UserData.LatexYLabel = '$\BLDCCurr$ in $\unit{A}$';
  spProp(4).UserData.LatexYLabel = '$\hat{\tau}_{\mathsf{MB}}$ in $\unit{Nm}$';
  spProp(4).UserData.LatexXLabel = '$t$ in $\unit{s}$';
  
  l2.UserData.LatexString = {'$e_{\rm{Obs}}$', '$e_{\rm{Ctrl}}$'};
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'PropCtrlRes', 'pdf', 0 );
  
end