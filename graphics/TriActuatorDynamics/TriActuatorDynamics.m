clear;

%% system param
Ts = 5e-3;

aA = 0.240;   % [m] distance geo center of main body to propeller axis
hA = 0.0;   % [m] distance body fixed frame to tilt axis
phi1 =  pi/3;   % [RAD] angle arm 1
phi2 =  pi;
phi3 = -pi/3;
epsP1 = -1;   % propeller 1 spinning direction
epsP2 = 1;
epsP3 = 1;

JP  = 3.57e-5;
pP1 = 6.855e-3;
pP2 = 5.154;
pP3 = 1.876e2;
kappaF = 1.402e-05;
kappaT = JP*pP1;
kappaTR = kappaT;

B = [0 0 0 sqrt(0.3e1) / 0.2e1 0 -sqrt(0.3e1) / 0.2e1; 0 0 0 -0.1e1 / 0.2e1 1 -0.1e1 / 0.2e1; 1 1 1 0 0 0; sqrt(0.3e1) * aA / 0.2e1 0 -sqrt(0.3e1) * aA / 0.2e1 (sqrt(0.3e1) * kappaTR + kappaF * hA) / kappaF / 0.2e1 -hA (sqrt(0.3e1) * kappaT + kappaF * hA) / kappaF / 0.2e1; -aA / 0.2e1 aA -aA / 0.2e1 -(-kappaF * sqrt(0.3e1) * hA + kappaTR) / kappaF / 0.2e1 -kappaT / kappaF -(kappaF * sqrt(0.3e1) * hA - kappaT) / kappaF / 0.2e1; kappaTR / kappaF -kappaT / kappaF -kappaT / kappaF -aA -aA -aA;];

w0S   = 40.11;     % identified [ServoIdent.m]
zetaS = 0.9154;    % identified [ServoIdent.m]
pS0 = w0S^2;
pS1 = 2*zetaS * w0S;

FxMax = 5;     % max horizontal force
FzMin = 5;      % min vertical force
FzMax = 20;     % max vertical force
tauxMax = 0.75; % min/max roll/pitch torque
tauzMax = 0.5;  % min/max yaw torque
fMax = [  FxMax;  FxMax; FzMax;  tauxMax;  tauxMax;  tauzMax];
fMin = [ -FxMax; -FxMax; FzMin; -tauxMax; -tauxMax; -tauzMax];

NormMat = diag(.5*(fMax-fMin));

%% force filter
lambdaF = 75;
kF = exp(-lambdaF*Ts);
GF = tf(1-kF, [1 -kF], Ts);

%% servo dynamics
GS = tf(pS0*Ts^2, [1 Ts*pS1-2 1-Ts*pS1+Ts^2*pS0], Ts);

%% bode plots
fig = figure(1); clf;
sp(1) = subplot(1,2,1);
sp(2) = subplot(1,2,2);

ww = {1e0, 1e3};

% transfer function for thrust and servo
[GS_mag, GS_phase, GS_w] = bode(GS, ww);
[GF_mag, GF_phase, GF_w] = bode(GF, ww);

axes(sp(1));
linehandle(1) = semilogx(GF_w, 20*log10(squeeze(GF_mag)), 'Color', 0.75*[1 1 1], 'linestyle', '--', 'linewidth', 3);
hold on;
linehandle(2) = line(GS_w, 20*log10(squeeze(GS_mag)), 'Color', 0.75*[1 1 1], 'linestyle', '-.',  'linewidth', 3);
grid on;
sp(1).YLim = [-45 5];
sp(1).YTick = -50:10:0;
line(pi/Ts*[1 1], ylim, 'Color', 'k');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');

%%% straigt hover case
theta0 = [0; 0; 0];
J = inv(NormMat) * B * [ diag(cos(theta0)) -diag(sin(theta0)); diag(sin(theta0)) diag(cos(theta0)) ];
G = J * ss(blkdiag(GF, GF, GF, GS, GS, GS)) * inv(J);
[G_mag, G_phase, G_w] = bode(G, ww);

axes(sp(1));
for i=1:6
  for j=1:6
    if (i==j)
      linehandle(3) = line(G_w, 20*log10(squeeze(G_mag(i,j,:))), 'Color', 'r');
    else
      linehandle(4) = line(G_w, 20*log10(squeeze(G_mag(i,j,:))), 'Color', 'b');
    end
  end
end

axes(sp(2));
linehandle(1) = semilogx(GF_w, 20*log10(squeeze(GF_mag)), 'Color', 0.75*[1 1 1], 'linestyle', '--', 'linewidth', 3);
hold on;
linehandle(2) = line(GS_w, 20*log10(squeeze(GS_mag)), 'Color', 0.75*[1 1 1], 'linestyle', '-.',  'linewidth', 3);
grid on;
sp(2).YLim = [-45 5];
sp(2).YTick = -60:10:0;
line(pi/Ts*[1 1], ylim, 'Color', 'k');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');

%%% straigt forward movement
theta0 = pi/4*[1; 0; -1];
J = inv(NormMat) * B * [ diag(cos(theta0)) -diag(sin(theta0)); diag(sin(theta0)) diag(cos(theta0)) ];
G = J * ss(blkdiag(GF, GF, GF, GS, GS, GS)) * inv(J);
[G_mag, G_phase, G_w] = bode(G, ww);

axes(sp(2));
for i=1:6
  for j=1:6
    if (i==j)
      linehandle(3) = line(G_w, 20*log10(squeeze(G_mag(i,j,:))), 'Color', 'r');
    else
      linehandle(4) = line(G_w, 20*log10(squeeze(G_mag(i,j,:))), 'Color', 'b');
    end
%    fprintf('(%i,%i) max magnitude %f\n', i, j, max(20*log10(squeeze(G_mag(i,j,:)))));
  end
end
legend(linehandle, 'thrust', 'servo', 'diagonal entries', 'off-diagonal entries');
drawnow;

%% export figure
figSize = [16 8];
fig.Units = 'centimeters';
fig.Position = [0 0 figSize(1) figSize(2)];
tightfig(fig);
fig.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig, 'BodeTriActorDyn', 'pdf', 0 );
