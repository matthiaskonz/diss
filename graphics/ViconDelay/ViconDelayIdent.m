clear;

Ts = 5e-3;

load ViconDelayIdent.mat
AgentNr = 1;

t = meas{AgentNr}.vSerial(1,:);
kk = find( (t>0) & (t<10000) );
t = t(kk);

pENC = meas{AgentNr}.vSerial([2;4;6],kk);
pIMU = meas{AgentNr}.vSerial([3;5;7],kk);
qVIC = meas{AgentNr}.vSerial(8:11,kk);

tVIC = meas{AgentNr}.vSerial(15,kk);
tREC = meas{AgentNr}.vSerial(17,kk);
TPC = meas{AgentNr}.vSerial(18,kk) + 30e-3; % correct old estimate

pVIC = Quaternion2RollPitchYaw(qVIC);

% remove stationary difference between euler angles
k0 = find(t >= 0, 1, 'first');
for i=1:3
  pENC(i,:) = pENC(i,:) - pENC(i,k0) + pVIC(i,k0);
  pIMU(i,:) = pIMU(i,:) - pIMU(i,k0) + pVIC(i,k0);
end

% remove redundant vicon frames
[tREC, kk, ~] = unique(tREC);
tVICp = tVIC(kk);
pVICp = pVIC(:,kk);

TPC0 = mean(tREC - tVICp);

%% display measurement
fig1 = figure(1); clf;
sp(1) = subplot(3,1,1); hold on;
line(t, pENC(1,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
stairs(tREC, pVICp(1,:), 'Color', 'k');
line(t, pIMU(1,:), 'Color', 'g');
legend('encoder', 'Vicon', 'IMU');
grid on;
sp(2) = subplot(3,1,2); hold on;
line(t, pENC(2,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
stairs(tREC, pVICp(2,:), 'Color', 'k');
line(t, pIMU(2,:), 'Color', 'g');
grid on;
sp(3) = subplot(3,1,3); hold on;
line(t, pENC(3,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
stairs(tREC, pVICp(3,:), 'Color', 'k');
line(t, pIMU(3,:), 'Color', 'g');
grid on;
linkaxes(sp, 'x');

%%

figure(2); clf; hold on;
line(tREC, tREC - tVICp, 'Color', 'r', 'Linestyle', 'none', 'Marker', '.');
stairs(t, TPC, 'Color', 'b');
line(xlim, TPC0*[1 1], 'Linestyle', '-.', 'Color', 'k');
line(xlim, (TPC0+Ts)*[1 1], 'Linestyle', '--', 'Color', 'k');
line(xlim, (TPC0-Ts)*[1 1], 'Linestyle', '--', 'Color', 'k');
legend('tREC - tVICp', 'TPCh', 'mean(tREC - tVICp)');

tVICp = tVICp + TPC0 + 0.015;
tREC = tREC + 0.015;

%% oversample data for convolution
TsS = 1e-3;
tS = t(100):TsS:t(end-100);
pENCS = interp1(t', pENC', tS, 'pchip')';
pVICS = interp1(tVICp', pVICp', tS, 'pchip')';
pVIC2S = interp1(tREC', pVICp', tS, 'pchip')';
pIMUS = interp1(t', pIMU', tS)';

taT = -(length(tS)-1)*TsS : TsS : (length(tS)-1)*TsS;
autoCorrENC = xcorr(pENCS(3,:), pENCS(3,:), 'coeff');
autoCorrVIC = xcorr(pVICS(3,:), pVICS(3,:), 'coeff');
xCorrVICENC = xcorr(pVICS(3,:), pENCS(3,:), 'coeff');
%aT2 = xcorr(pVIC2S(3,:), pENCS(3,:), 'coeff');
[~, k] = max(xCorrVICENC);
TVICD = taT(k)

figure(3); clf;
line(taT, autoCorrENC, 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
line(taT, autoCorrVIC, 'Color', 'k', 'Linestyle', '--');
line(taT, xCorrVICENC, 'Color', 'r');
%line(taT, aT2, 'Color', 'b');
legend('autocorr. enc', 'autocorr. Vicon', 'xcorr', 'Location', 'Southwest');

grid on;
xlabel('time in s');
ylabel('xcorr.');
xlim([-5*TVICD 5*TVICD]);
ylim([0.9 1.001]);
line(TVICD*[1 1], ylim, 'Color', 'k', 'Linestyle', '--');

title(sprintf('identified delay T_D = %1.4f s', TVICD), 'interpreter', 'none');

%%
axes(sp(1));
%line(tREC-TVICD, pVICp(1,:), 'Color', 'm', 'Linestyle', '-', 'Marker', '.');
line(tVICp-TVICD, pVICp(1,:), 'Color', 'r', 'Linestyle', '-', 'Marker', '.');
axes(sp(2));
%line(tREC-TVICD, pVICp(2,:), 'Color', 'm', 'Linestyle', '-', 'Marker', '.');
line(tVICp-TVICD, pVICp(2,:), 'Color', 'r', 'Linestyle', '-', 'Marker', '.');
axes(sp(3));
%line(tREC-TVICD, pVICp(3,:), 'Color', 'm', 'Linestyle', '-', 'Marker', '.');
line(tVICp-TVICD, pVICp(3,:), 'Color', 'r', 'Linestyle', 'none', 'Marker', '.');

%% figure for paper
fig5 = figure(5); clf;
subplot(1,2,1); grid on; hold on;
line(t, 180/pi*pENC(3,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
stairs(tREC, 180/pi*pVICp(3,:), 'Color', 'b');
line(tVICp-TVICD, 180/pi*pVICp(3,:), 'Color', 'r', 'Linestyle', 'none', 'Marker', '.');
%line(tS, 180/pi*pVICS(3,:), 'Color', 'g');
legend('Encoder', 'Vicon', 'corrected', 'location', 'southeast');
%xlim([249 254]);
xlim([34 36]);
ylabel('yaw angle in DEG');
xlabel('t in s');

subplot(1,2,2); grid on;
line(taT, autoCorrENC, 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', '--');
line(taT, autoCorrVIC, 'Color', 'b', 'Linestyle', '-');
line(taT, xCorrVICENC, 'Color', 'r');
legend('autocorr. Encoder', 'autocorr. Vicon', 'crosscorr.', 'Location', 'Southwest');
xlabel('s');
%xlim([-5*TVICD 5*TVICD]);
xlim([-0.1 0.1]);
ylim([0.9 1.00]);
set(gca, 'ytick', [0.9 0.95 1]);
line(TVICD*[1 1], ylim, 'Color', 'k', 'Linestyle', '-.');
line([0 0], ylim, 'Color', 'k', 'Linestyle', '-.');

%% export figure
figSize = [16 8];
fig5.Units = 'centimeters';
fig5.Position = [0 0 figSize(1) figSize(2)];
tightfig(fig5);
fig5.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig5, 'ViconDelayIdent', 'pdf', 0 );

