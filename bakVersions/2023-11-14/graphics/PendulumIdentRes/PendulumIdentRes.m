clear;

AgentNr = 1;

AX = 1;

tS = 0;
tF = 9999;

Ts = 5e-3;
m = 1.001;
g = 9.81;

measFilename = 'Quad_Jx.mat';
l = 0.982;
d = 0.760;

load(measFilename);
meas = meas{AgentNr};

kk = find((meas.vSerial(1,:) >= tS) & (meas.vSerial(1,:) <= tF));
tM = meas.vSerial(1,kk);
t = [tM(5) : Ts : tM(end-5)];
q = interp1(tM', meas.vSerial([2; 4; 6],kk)', t')';
w = interp1(tM', meas.vSerial([3 5 7],kk)', t'+Ts)';
rd = interp1(tM', meas.vSerial([12 14 16],kk)', t')';
a = interp1(tM', meas.vSerial([13 15 17],kk)', t'+Ts)';
iM = interp1(tM', meas.vSerial([8:11],kk)', t')';
wP = interp1(tM', meas.vSerial([18:21],kk)', t'+1.5*Ts)';

q = [sqrt(1 - [1 1 1]*(q.^2)); q];

N = length(t);
phi = zeros(1, N);
for k = 1:N
  R = Quaternion2RotMat(q(:,k));
  phi(k) = -asin(R(1,2));
end

phi = phi - mean(phi);

%%
SG_ORDER = 2;
SG_WINDOW = 31;
[tSG, phiSG] = sgolayDiffEstim(t, phi, SG_ORDER, SG_WINDOW, Ts, 2, 0);
[tSG, wSG] = sgolayDiffEstim(t, w(AX,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);


%%
% figMeas = figure(1); clf;
% spMeas(1) = subplot(3,1,1); grid on;
% line(t, 180/pi*phi, 'Color', 'k', 'Linewidth', 1);
% line(tSG, 180/pi*phiSG(1,:), 'Color', 'b', 'Linewidth', 1);
% spMeas(2) = subplot(3,1,2); grid on;
% line(t, w(AX,:), 'Color', 'k', 'Linewidth', 1);
% line(tSG, phiSG(2,:), 'Color', 'b', 'Linewidth', 1);
% spMeas(3) = subplot(3,1,3); grid on;
% line(tSG, wSG(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3);
% %line(tSG, phiSG(3,:), 'Color', 'b', 'Linewidth', 1);
% xlabel('time in s');
% 
% linkaxes(spMeas, 'x');

%%
A = wSG(2,:)';
b = m*g*d^2/4 * sin(phiSG(1,:)') ./ sqrt(l^2 - d^2/2*(1-cos(phiSG(1,:)')));

J = -A \ b

phidd = -m*g*d^2/4/J * sin(phiSG(1,:)) ./ sqrt(l^2 - d^2/2*(1-cos(phiSG(1,:))));


%%
fig = figure(1); clf;
sp = subplot(1,1,1);
line(tSG-100, wSG(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(tSG-100, phidd, 'Color', 'b', 'Linewidth', 1, 'Linestyle', '-');
legend('measurement', 'model');
xlim([0 3.5]);
grid on;
ylabel('add in RAD/s^2', 'Interpreter', 'none');
xlabel('time in s');


saveFigure = 1;
figSize = [9 6.3];

if (saveFigure > 0)
  
  sp.UserData.LatexYLabel = '$\ddot{\alpha}$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
  sp.UserData.LatexXLabel = '$t$ in $\unit{s}$';
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'PendulumIdentRes', 'pdf', 0 );
  
end
