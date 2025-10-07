clear;

validate = 0;
saveFigure = 1;

Ts = 0.005;   % sampling time

% bessel filter param
N = 4;        % filter order
w0 = 2*pi*7;  % cut-off frequency

% gauss filter param
FG = 49;
SG = 0.123;
TG = 0.07;

tS = 36.5;
tF = 38.5;
t0 = 31.5;

%%
load export4.mat

t = meas{1}.vSerial(1,:);
kk = find((t > 0) & (t < 1000));
t = t(kk) - t0;

u = meas{1}.vSerial(2,kk);
x = meas{1}.vSerial(3:6,kk);

%% IIR Bessel filter

[sysB, Ad, Bd] = JokerBesselFilterDiffEst(N, w0, Ts);

% impulse response
[BF, tBF] = impulse(sysB, 0.3);
BF = BF';
tBF = tBF';

%% FIR Gauss filter
[ tG, xG, tGF, GF ] = GaussDiffEstim( t, u, FG, SG, Ts, 0);

tG = tG+TG;
tGF = tGF+TG;

%%
fig = figure(1); clf;
spImpulse(1) = subplot(5,2,[1,3]);
line(tGF, GF(1,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(tBF, BF(1,:), 'Color', 'b');
line(TG*[1 1], ylim, 'color', 'k', 'Linestyle', '-.');
ylim([-1 15]);
xlim([min([tGF tBF]), max([tGF tBF])]);
ylabel('$y$', 'interpreter', 'LaTeX');
grid on;
title('impulse response', 'interpreter', 'LaTeX');

spImpulse(2) = subplot(5,2,5);
line(tGF, GF(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(tBF, BF(2,:), 'Color', 'b');
line(TG*[1 1], ylim, 'color', 'k', 'Linestyle', '-.');
ylabel('$\dot{y}$', 'interpreter', 'LaTeX');
grid on;

spImpulse(3) = subplot(5,2,7);
line(tGF, GF(3,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(tBF, BF(3,:), 'Color', 'b');
line(TG*[1 1], ylim, 'color', 'k', 'Linestyle', '-.');
ylabel('$\ddot{y}$', 'interpreter', 'LaTeX');
grid on;
spImpulse(4) = subplot(5,2,9);
line(tGF, GF(4,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(tBF, BF(4,:), 'Color', 'b');
line(TG*[1 1], ylim, 'color', 'k', 'Linestyle', '-.');
ylabel('$y^{(3)}$', 'interpreter', 'LaTeX');
grid on;
ylim(1e6*[-2 2]);
xlabel('$t$ in s', 'interpreter', 'LaTeX');
linkaxes(spImpulse, 'x');

spMeas(1) = subplot(5,2,[2,4]);
l(1) = stairs(t, u, 'color', 'k');
box off
l(2) = line(tG, xG(1,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
l(3) = line(t, x(1,:), 'Color', 'b');
ylim([-1.1 1.1]);
legend([l(1) l(3) l(2)], {'input', 'Bessel', 'Gauss'}, 'location', 'southwest');
grid on;
title('application', 'interpreter', 'LaTeX');

spMeas(2) = subplot(5,2,6);
line(tG, xG(2,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(t, x(2,:), 'Color', 'b');
ylim(20*[-1 1]);
grid on;

spMeas(3) = subplot(5,2,8);
line(tG, xG(3,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(t, x(3,:), 'Color', 'b');
ylim(300*[-1 1]);
grid on;

spMeas(4) = subplot(5,2,10);
line(tG, xG(4,:), 'Color', 0.7*[1 1 1], 'Linewidth', 3, 'Linestyle', ':');
line(t, x(4,:), 'Color', 'b');
ylim(1.2e4*[-1 1]);
grid on;
xlabel('$t$ in s', 'interpreter', 'LaTeX')

linkaxes(spMeas, 'x');
xlim([tS-t0 tF-t0]);

%% validate
if (validate > 0)
  xS = zeros(N+1, length(t));
  for k=2:length(t);
    xS(:,k) = Ad*xS(:,k-1) + Bd*u(k);
  end
  
  for i = 1:N
    axes(spMeas(i));
    line(t, xS(i,:), 'Color', 'r', 'Linestyle', '--');
  end
end

%% save figure
if (saveFigure > 0)
  figSize = [16 16];
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'RCFilter', 'pdf', 0 );
end