clear;

saveFigure = 0;

choice = 1;
AgentNr = 1;

%iBLDC2A = 0.0177;

if (choice == 1)
  measFilename = 'export4.mat'; % pos transition
  tS = 235;
  tF = tS+35;
  r0 = [0; 0.4; 0];
  rYLim = [-0.6 1.6];
  correctOffset = 0;
elseif (choice == 2)
  measFilename = 'loop1.mat';
  tS = 362;
  tF = tS+35;
  r0 = [-1.27; 0.222; 0.685];
  rYLim = [-1.1 3.1];
  correctOffset = 1;
elseif (choice == 3)
  measFilename = 'flip14.mat'; % pos transition
  tS = 711;
  tF = tS+35;
  r0 = [0.058; 0.102; 0.22];
  rYLim = [-0.3 2.0];
  correctOffset = 0;
end

load(measFilename);
tv = meas{AgentNr}.vSerial(1,:);
kk = find( (tv>=tS-2) & (tv<=tF+2) );
tv = tv(kk);
v = meas{AgentNr}.vSerial(2:end,kk);

TPC = meas{AgentNr}.ViconShiftEstim;
tV = meas{AgentNr}.vVicon(1,:) + TPC;
vV = meas{AgentNr}.vVicon(2:end,:);

Ts = 0.005;
t = tS:Ts:tF;

RPYR = interp1(tv', v([1 3 5],:)',    t')';
RPY  = interp1(tv', v([2 4 6],:)',    t')';
rR   = interp1(tv', v([11 13 15],:)', t')';
r    = interp1(tv', v([12 14 16],:)', t')';
iM   = interp1(tv', v([7:10],:)',     t')';
aS   = interp1(tv', v([17:20],:)',    t')';

t = t - t(1);

% move reference
for i = 1:3
  r(i,:) = r(i,:) + r0(i);
  rR(i,:) = rR(i,:) + r0(i);
end

if correctOffset
  rMean = mean(r, 2);
  rRMean = mean(rR, 2);
  r(1,:) = r(1,:) - rMean(1) + rRMean(1);
  r(2,:) = r(2,:) - rMean(2) + rRMean(2);
  r(3,:) = r(3,:) - rMean(3) + rRMean(3);
end

%% configuration error
posErr = sqrt([1 1 1] * (r - rR).^2);

if (choice == 2)
  posErr = posErr*0.85;
end

angleErr = zeros(1,length(t));
headErr = zeros(1,length(t));
for k = 1:length(t)
  R = RollPitchYaw2RotMat(RPY(:,k));
  RR = RollPitchYaw2RotMat(RPYR(:,k));
  RE = RR'*R;
  angleErr(k) = acos(.5*(trace(RE)-1));
  headErr(k) = abs(atan2(RE(2,1)-RE(1,2), RE(1,1)+RE(2,2)));
end

posErrRMS = rms(posErr);
angleErrRMS = rms(angleErr);
headErrRMS = rms(headErr);
%fprintf('%s: position error RMS %1.1f mm, attitude error RMS %1.1f DEG\n', measFilename, 1e3*posErrRMS, 180/pi*headErrRMS);

%% velocity and acceleration
rRd = diff(rR, 1, 2)/Ts;
rRdd = diff(rRd, 1, 2)/Ts;
wR = zeros(3,length(t)-1);
RR0 = RollPitchYaw2RotMat(RPYR(:,1));
for k = 1:length(t)-1
  RR1 = RollPitchYaw2RotMat(RPYR(:,k+1));
  wR(:,k) = vee_so3(RR0'*(RR1-RR0)/Ts);
  RR0 = RR1;
end
wRd = diff(wR,1,2)/Ts;

%%
numPlots = 6;
fig = figure(choice); clf;

spMeas(1) = subplot(numPlots,1,1:2); grid on;
l(1) = line(t, rR(1,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(2) = line(t, rR(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(3) = line(t, rR(3,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(4) = line(t, r(1,:), 'Color', 'r', 'Linewidth', 1);
l(5) = line(t, r(2,:), 'Color', 'g', 'Linewidth', 1);
l(6) = line(t, r(3,:), 'Color', 'b', 'Linewidth', 1);
ylabel('m', 'Interpreter', 'LaTeX');
legend(l([4:6, 1]), {'$r_{\rm{x}}$', '$r_{\rm{y}}$', '$r_{\rm{z}}$', 'ref.'}, 'Interpreter', 'LaTeX');
%ylim(rYLim);
%axis tight;

spMeas(2) = subplot(numPlots,1,3); grid on;
l(1) = line(t, posErr, 'Color', 'k', 'Linewidth', 1);
l(2) = line([t(1) t(end)], posErrRMS*[1 1], 'Color', 'k', 'Linewidth', 1, 'Linestyle', '-.');
legend(sprintf('pos error RMS %1.3f m', posErrRMS), 'location', 'northeast');
ylabel('m', 'Interpreter', 'LaTeX');
legend(l(1:2), {'position error', sprintf('RMS = %1.3f m', posErrRMS)}, 'location', 'northeast', 'Interpreter', 'LaTeX');
%ylim([0 0.1]);
%axis tight;

spMeas(3) = subplot(numPlots,1,4:5); grid on;
l(1) = line(t, 180/pi*RPYR(1,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(2) = line(t, 180/pi*RPYR(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(3) = line(t, 180/pi*RPYR(3,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
l(4) = line(t, 180/pi*RPY(1,:), 'Color', 'r', 'Linewidth', 1);
l(5) = line(t, 180/pi*RPY(2,:), 'Color', 'g', 'Linewidth', 1);
l(6) = line(t, 180/pi*RPY(3,:), 'Color', 'b', 'Linewidth', 1);
ylabel('deg', 'Interpreter', 'LaTeX');
legend(l([4,5,6,1]), {'roll', 'pitch', 'yaw', 'ref'}, 'location', 'northeast', 'Interpreter', 'LaTeX');
%ylim([-45 95]);
%spMeas(3).YTick = [-180:45:180];

spMeas(4) = subplot(numPlots,1,6); grid on;
l(1) = line(t, 180/pi*angleErr, 'Color', 'b', 'Linewidth', 1, 'Linestyle', '-');
l(2) = line(t, 180/pi*headErr, 'Color', 'k', 'Linewidth', 1, 'Linestyle', '-');
l(3) = line([t(1) t(end)], 180/pi*headErrRMS*[1 1], 'Color', 'k', 'Linewidth', 1, 'Linestyle', '-.');
ylabel('deg', 'Interpreter', 'LaTeX');
legend(l(1:3), {'attitude error', 'heading error',  sprintf('RMS = $%1.1f$ deg', 180/pi*headErrRMS)}, 'location', 'northeast', 'Interpreter', 'LaTeX');
%ylim([0 10]);
%axis tight;

xlabel('time in s', 'Interpreter', 'LaTeX');
linkaxes(spMeas, 'x');
xlim([t(1) t(end)]);


%%
if 1
  axes(spMeas(1));
  line(t(1:end-1), sqrt([1 1 1]*rRd.^2), 'Color', 'm');
  line(t(1:end-2), sqrt([1 1 1]*rRdd.^2), 'Color', 'c');
  axes(spMeas(3));
  line(t(1:end-1), 180/pi*sqrt([1 1 1]*wR.^2), 'Color', 'm');
  line(t(1:end-2), 180/pi*sqrt([1 1 1]*wRd.^2), 'Color', 'c');
end

%%
if (saveFigure > 0)
  
  figSize = [16 23];
%   sp_xi(1).UserData.LatexYLabel = '$\vxd$ in $\sfrac{\unit{m}}{\unit{s}^2}$';
%   sp_xi(2).UserData.LatexYLabel = '$\vzd$ in $\sfrac{\unit{m}}{\unit{s}^2}$';
%   sp_xi(3).UserData.LatexYLabel = '$\wyd$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
%   sp_xi(4).UserData.LatexYLabel = '$\wzd$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
%   sp_xi(5).UserData.LatexYLabel = '$\PropVeld[1]$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
%   sp_xi(5).UserData.LatexXLabel = '$t$ in $\unit{s}$';
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
if (choice == 1)
  JokerPrintFig( fig, 'QuadManeuver42Result', 'pdf', 0 );
elseif (choice == 2)
  JokerPrintFig( fig, 'QuadLoopResult', 'pdf', 0 );
elseif (choice == 3)
  JokerPrintFig( fig, 'QuadFlipResult', 'pdf', 0 );
end
end