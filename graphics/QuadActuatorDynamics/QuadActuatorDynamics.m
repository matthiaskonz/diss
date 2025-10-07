

% model param
Ts = 5e-3;

aA = 0.24;  % [m] distance geo center of main body to propeller axis
JP = 3.7e-5;        % [kg m^2] CAD model 2018-03-29
pP1 = 6.855e-3;
pP2 = 5.154;
pP3 = 1.876e2;
kappaF = 1.455e-05;
kappaT = JP*pP1;
kM     = JP*pP2;
tau0   = JP*pP3;

% continuous design param
wHover = sqrt(9.81/4/kappaF);
% filter design param
w0_Z    = 100;  zeta_Z    = sqrt(2)/2;
w0_Tilt = 125;  zeta_Tilt = sqrt(2)/2;
w0_Head = 66;

P0 = eye(4);
P1 = diag([ 2*zeta_Z/w0_Z; 2*zeta_Tilt/w0_Tilt; 2*zeta_Tilt/w0_Tilt; 1/w0_Head ]);
P2 = diag([ 1/w0_Z^2; 1/w0_Tilt^2; 1/w0_Tilt^2; 0 ]);

s = tf('s');
P = P0 + s*P1 + s^2*P2;
G = inv(P);

% corresponding time discrete filter
Gd = c2d(G, Ts, 'zoh');
[~, a_Z] = tfdata(Gd(1,1), 'v');
[~, a_Tilt] = tfdata(Gd(2,2), 'v');
[~, a_Head] = tfdata(Gd(4,4), 'v');

% resulting filter gains
p4 = JP/kappaT/wHover/2;
k0_Z = (a_Z(3) + a_Z(2) + 1) / Ts ^ 2;
k1_Z = (a_Z(2) + 2) / Ts;
k0_Tilt = (a_Tilt(3) + a_Tilt(2) + 1) / Ts ^ 2;
k1_Tilt = (a_Tilt(2) + 2) / Ts;
k0_Head = 1 / p4 * (a_Head(2) + 1) / Ts;
k1_Head = 1 / p4;

% return values
k0 = [k0_Z; k0_Tilt; k0_Tilt; k0_Head];
k1 = [k1_Z; k1_Tilt; k1_Tilt; k1_Head];

% min max body fixed forces and torques
FmMin = 4;      % min thrust magnitude
FmMax = 16;     % max thrust magnitude
tauxMax = 1.2;  % min/max roll/pitch torque
tauzMax = 1.0; % min/max yaw torque
fMax = [FmMax;  tauxMax;  tauxMax;  tauzMax];
fMin = [FmMin; -tauxMax; -tauxMax; -tauzMax];
NormMat = diag(.5*(fMax-fMin));

%% figure
fig = figure(1); clf;
sp(1) = subplot(1,2,1);
sp(1).XScale = 'log';
sp(1).YLim = [-30 5];
sp(1).XLabel.String = 'Frequency (rad/s)';
sp(1).YLabel.String = 'Magnitude (dB)';
line(pi/Ts*[1 1], ylim, 'Color', 'k', 'Linewidth', 1);
grid on;
sp(2) = subplot(1,2,2);
sp(2).XScale = 'log';
sp(2).YLim = [-30 5];
sp(2).XLabel.String = 'Frequency (rad/s)';
line(pi/Ts*[1 1], ylim, 'Color', 'k', 'Linewidth', 1);
grid on;

ww = {1e0, 1e3};

% continus design
%[G_mag, G_phase, G_w] = bode(G, ww);
%axes(sp(1));
%line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-.');
%line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-.');
%line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-.');

% discrete at hover
wP10 = wHover;
wP20 = wHover;
wP30 = wHover;
wP40 = wHover;

z = tf('z', Ts);
s = (z-1)/Ts;

p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;
G_Z = k0_Z / (s * k1_Z + s ^ 2 + k0_Z);
G_Tilt = k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt);
G_Head = (p4 * k0_Head * s + k0_Head) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head);
G_41 = p1 * k0_Z / (s * k1_Z + s ^ 2 + k0_Z) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);
G_42 = p2 * k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);
G_43 = p3 * k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);

G = [G_Z 0 0 0; 0 G_Tilt 0 0; 0 0 G_Tilt 0; G_41 G_42 G_43 G_Head];
[G_mag, G_phase, G_w] = bode(G, ww);

axes(sp(1));
line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-', 'Linewidth', 2);
line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-', 'Linewidth', 2);
line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-', 'Linewidth', 2);
axes(sp(2));
l(1) = line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-', 'Linewidth', 2);
l(2) = line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-', 'Linewidth', 2);
l(3) = line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-', 'Linewidth', 2);
drawnow;

ww0 = sqrt([1 6]/kappaF);
for wP10 = ww0
  for wP20 = ww0
    for wP30 = ww0
      for wP40 = ww0
        p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
        p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
        p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
        p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;
        G_Z = k0_Z / (s * k1_Z + s ^ 2 + k0_Z);
        G_Tilt = k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt);
        G_Head = (p4 * k0_Head * s + k0_Head) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head);
        G_41 = p1 * k0_Z / (s * k1_Z + s ^ 2 + k0_Z) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);
        G_42 = p2 * k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);
        G_43 = p3 * k0_Tilt / (s * k1_Tilt + s ^ 2 + k0_Tilt) / (s ^ 2 + (p4 * k0_Head + k1_Head) * s + k0_Head) * s ^ 2 * (k1_Head + s);
        
        G = [G_Z 0 0 0; 0 G_Tilt 0 0; 0 0 G_Tilt 0; G_41 G_42 G_43 G_Head];
        G = inv(NormMat) * G * NormMat;
        [G_mag, G_phase, G_w] = bode(G, ww);
        
        axes(sp(2));
        l(4) = line(G_w, 20*log10(squeeze(G_mag(4,1,:))), 'Color', 'm', 'Linestyle', '-');
        l(5) = line(G_w, 20*log10(squeeze(G_mag(4,2,:))), 'Color', 'c', 'Linestyle', '-');
        %line(G_w, 20*log10(squeeze(G_mag(4,3,:))), 'Color', 'k', 'Linestyle', '--');
        line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-');
        
        drawnow;
      end
    end
  end
end

legend(l, '$G_{11}$', '$G_{22}$', '$G_{44}$', '$G_{41}$', '$G_{42}$', 'Location', 'southwest');

%%
figSize = [16 8];
fig.Units = 'centimeters';
fig.Position = [0 0 figSize(1) figSize(2)];
tightfig(fig);
fig.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig, 'BodeQuadActorDyn', 'pdf', 0 );

%%
fig = figure(2); clf;
sp(1) = subplot(1,2,1);
sp(1).XScale = 'log';
sp(1).YLim = [-30 25];
sp(1).XLabel.String = 'Frequency (rad/s)';
sp(1).YLabel.String = 'Magnitude (dB)';
line(pi/Ts*[1 1], ylim, 'Color', 'k', 'Linewidth', 1);
grid on;
sp(2) = subplot(1,2,2);
sp(2).XScale = 'log';
sp(2).YLim = [-30 25];
sp(2).XLabel.String = 'Frequency (rad/s)';
line(pi/Ts*[1 1], ylim, 'Color', 'k', 'Linewidth', 1);
grid on;

% discrete at hover
wP10 = wHover;
wP20 = wHover;
wP30 = wHover;
wP40 = wHover;

p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;

zeta_Head = sqrt(2)/2;
tmp = c2d(tf(1, [1/w0_Head^2 2*zeta_Head/w0_Head 1]), Ts, 'zoh');
[~, blubb] = tfdata(tmp, 'v');
k0_Head = (blubb(3) + blubb(2) + 1) / Ts ^ 2;
k1_Head = (blubb(2) + 2) / Ts;

G = [k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) 0 0 0; 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0 0; 0 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0; p1 * s * k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) p2 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) p3 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) k0_Head * (p4 * s + 1) / (s * k1_Head + s ^ 2 + k0_Head);];
G = inv(NormMat) * G * NormMat;
[G_mag, G_phase, G_w] = bode(G, ww);

axes(sp(1));
l(1) = line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-', 'Linewidth', 2);
l(2) = line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-', 'Linewidth', 2);
l(3) = line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-', 'Linewidth', 2);

for wP10 = ww0
  for wP20 = ww0
    for wP30 = ww0
      for wP40 = ww0
        p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
        p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
        p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
        p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;
        
        G = [k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) 0 0 0; 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0 0; 0 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0; p1 * s * k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) p2 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) p3 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) k0_Head * (p4 * s + 1) / (s * k1_Head + s ^ 2 + k0_Head);];
        G = inv(NormMat) * G * NormMat;
        [G_mag, G_phase, G_w] = bode(G, ww);
        
        axes(sp(1));
        line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-');
        line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-');
        l(4) = line(G_w, 20*log10(squeeze(G_mag(4,1,:))), 'Color', 'm', 'Linestyle', '-');
        l(5) = line(G_w, 20*log10(squeeze(G_mag(4,2,:))), 'Color', 'c', 'Linestyle', '-');
        line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-');
                
        drawnow;
      end
    end
  end
end
%legend(l, '$G_{11}$', '$G_{22}$', '$G_{44}$', '$G_{41}$', '$G_{42}$', 'Location', 'southwest');

% discrete at hover
wP10 = wHover;
wP20 = wHover;
wP30 = wHover;
wP40 = wHover;

p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;

k0_Head = 1 / p4 * (a_Head(2) + 1) / Ts;
k1_Head = (a_Head(2) * p4 + Ts + p4) / Ts / p4;

G = [k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) 0 0 0; 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0 0; 0 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0; p1 * s * k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) p2 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) p3 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) k0_Head * (p4 * s + 1) / (s * k1_Head + s ^ 2 + k0_Head);];
G = inv(NormMat) * G * NormMat;
[G_mag, G_phase, G_w] = bode(G, ww);

axes(sp(2));
l(1) = line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-', 'Linewidth', 2);
l(2) = line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-', 'Linewidth', 2);
l(3) = line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-', 'Linewidth', 2);

for wP10 = ww0
  for wP20 = ww0
    for wP30 = ww0
      for wP40 = ww0
        p1 = JP / kappaF * (0.1e1 / wP40 - 0.1e1 / wP30 + 0.1e1 / wP20 - 0.1e1 / wP10) / 0.8e1;
        p2 = JP / aA / kappaF * (0.1e1 / wP20 - 0.1e1 / wP40) / 0.4e1;
        p3 = JP / aA / kappaF * (0.1e1 / wP10 - 0.1e1 / wP30) / 0.4e1;
        p4 = JP / kappaT * (0.1e1 / wP10 + 0.1e1 / wP20 + 0.1e1 / wP30 + 0.1e1 / wP40) / 0.8e1;
        
        G = [k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) 0 0 0; 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0 0; 0 0 k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) 0; p1 * s * k0_Z / (s * k1_Z + (s ^ 2) + k0_Z) p2 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) p3 * s * k0_Tilt / (s * k1_Tilt + (s ^ 2) + k0_Tilt) k0_Head * (p4 * s + 1) / (s * k1_Head + s ^ 2 + k0_Head);];
        G = inv(NormMat) * G * NormMat;
        [G_mag, G_phase, G_w] = bode(G, ww);
        
        axes(sp(2));
        line(G_w, 20*log10(squeeze(G_mag(1,1,:))), 'Color', 'r', 'Linestyle', '-');
        line(G_w, 20*log10(squeeze(G_mag(2,2,:))), 'Color', 'g', 'Linestyle', '-');
        l(4) = line(G_w, 20*log10(squeeze(G_mag(4,1,:))), 'Color', 'm', 'Linestyle', '-');
        l(5) = line(G_w, 20*log10(squeeze(G_mag(4,2,:))), 'Color', 'c', 'Linestyle', '-');
        line(G_w, 20*log10(squeeze(G_mag(4,4,:))), 'Color', 'b', 'Linestyle', '-');
                
        drawnow;
      end
    end
  end
end
legend(l, '$G_{11}$', '$G_{22}$', '$G_{44}$', '$G_{41}$', '$G_{42}$', 'Location', 'southwest');

%%
figSize = [16 8];
fig.Units = 'centimeters';
fig.Position = [0 0 figSize(1) figSize(2)];
tightfig(fig);
fig.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig, 'BodeQuadActorDynTraditional', 'pdf', 0 );
