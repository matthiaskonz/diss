clear;
filename = 'ServoTestStand5';
AgentNr = 1;

saveFigure = 1;
figSize = [15.9 10];

tW = 43.1;
TW = 2;

%% param
Ts = 5e-3;
kappaF = 1.467e-05;

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
spServo(1) = subplot(3,1,1:2); grid on; hold on;
l1 = line(t, 180/pi*a, 'Color', 0.7*[1 1 1], 'Linewidth', 3);
l2 = stairs(t, 180/pi*aR, 'Color', 'k', 'Linewidth', 1);
l3 = stairs(t, 180/pi*aS, 'Color', 'b', 'Linewidth', 1);
spServo(1).YTick = 45*[-1 0 1];
%stairs(t, 180/pi*aU, 'Color', 'r', 'Linewidth', 1);
legend([l2 l1 l3], 'input', 'measured', 'simulated', 'location', 'southwest');
ylabel('angle in DEG');
spServo(2) = subplot(3,1,3); grid on; hold on;
stairs(t, 180/pi*(a-aS), 'Color', 'b', 'Linewidth', 1);
%line(t, 180/pi*(a-aR), 'Color', 'r', 'Linewidth', 1);
ylabel('error in DEG');
xlabel('time in s');
linkaxes(spServo, 'x');
xlim([0 TW])


%%
if (saveFigure > 0)
  
  spServo(1).UserData.LatexYLabel = '$\aServo$ in $\unit{RAD}$';
  spServo(2).UserData.LatexYLabel = 'error in $\unit{RAD}$';
	spservo(2).UserData.LatexXLabel = '$t$ in $\unit{s}$';
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'ServoSimRes', 'pdf', 0 );
  
end