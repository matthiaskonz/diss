clear;

filename = '2018-08-30_TriV42_2.mat';
AgentNr = 2;

tS = 54.5;
tF = tS + 20;

saveFigure = 0;
figSize = [15.9 23.3];

% t0 = 1;
% tEnd = 1000;

%% param
aA = 0.24;
hA = 0;
g = 9.81;

m   = 1.251;        % [kg] measured 2018-03-29
Jxx = 19.2e-3;      % [kg m^2] measured with pendulum
Jyy = 19.2e-3;      % [kg m^2] measured with pendulum
Jzz = 30.7e-3;      % [kg m^2] measured with pendulum

% identified propeller param
pP1 = 6.75e-3;
pP2 = 5.2;
pP3 = 100;

Ts = 5e-3;

%%
load(filename);
meas = meas{AgentNr};

kk = find((meas.vSerial(1,:) >= tS) & (meas.vSerial(1,:) <= tF));
tM = meas.vSerial(1,kk);
t = [tM(5):Ts:tM(end-5)];
q = interp1(tM', meas.vSerial([21; 2; 4; 6],kk)', t')';
w = interp1(tM', meas.vSerial([3 5 7],kk)', (t-Ts)')';
rd = interp1(tM', meas.vSerial([12 14 16],kk)', t')';
aIMU = interp1(tM', meas.vSerial([13 15 17],kk)', (t+Ts)')';
wP = interp1(tM', meas.vSerial([8:10],kk)', t'+1.5*Ts)';
pP = interp1(tM', meas.vSerial([18:20],kk)', t')';

iM1 = interp1(tM', meas.vSerial(11,kk)', t')';

t = t-t(1);

%% SG filter and derivative estimation
SG_ORDER = 5;
SG_WINDOW = 31;
[tSG, w1SG] = sgolayDiffEstim(t, w(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, w2SG] = sgolayDiffEstim(t, w(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, w3SG] = sgolayDiffEstim(t, w(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q0SG] = sgolayDiffEstim(t, q(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q1SG] = sgolayDiffEstim(t, q(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q2SG] = sgolayDiffEstim(t, q(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q3SG] = sgolayDiffEstim(t, q(4,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r1dSG] = sgolayDiffEstim(t, rd(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r2dSG] = sgolayDiffEstim(t, rd(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r3dSG] = sgolayDiffEstim(t, rd(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP1SG] = sgolayDiffEstim(t, wP(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP2SG] = sgolayDiffEstim(t, wP(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP3SG] = sgolayDiffEstim(t, wP(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, pP1SG] = sgolayDiffEstim(t, pP(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, pP2SG] = sgolayDiffEstim(t, pP(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, pP3SG] = sgolayDiffEstim(t, pP(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, iM1SG] = sgolayDiffEstim(t, iM1, SG_ORDER, SG_WINDOW, Ts, 2, 0);

qSG = [q0SG(1,:); q1SG(1,:); q2SG(1,:); q3SG(1,:)];

wSG   = [w1SG(1,:); w2SG(1,:); w3SG(1,:)];
wdSG  = [w1SG(2,:); w2SG(2,:); w3SG(2,:)];

rdSG  = [r1dSG(1,:); r2dSG(1,:); r3dSG(1,:)];
rddSG = [r1dSG(2,:); r2dSG(2,:); r3dSG(2,:)];

N = length(tSG); % number of measurement points
xiSG = zeros(6, N);
xidSG = zeros(6, N);
for k = 1:N
  qSG(:,k) = qSG(:,k)/norm(qSG(:,k));
  R = Quaternion2RotMat(qSG(:,k));
  xiSG(:,k) = [R'*rdSG(:,k); wSG(:,k)];
  xidSG(:,k) = [R'*rddSG(:,k) - wed_so3(xiSG(4:6,k))*xiSG(1:3,k); wdSG(:,k)];
end

%%
% figure(2); clf;
% sp_xid(1) = subplot(3,2,1); grid on;
% line(tSG, xidSG(1,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('vxd in m/s^2');
% sp_xid(2) = subplot(3,2,3); grid on;
% line(tSG, xidSG(2,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('vyd in m/s^2');
% sp_xid(3) = subplot(3,2,5); grid on;
% line(tSG, xidSG(3,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('vzd in m/s^2');
% xlabel('t in s');
% sp_xid(4) = subplot(3,2,2); grid on;
% line(tSG, xidSG(4,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('wxd in RAD/s^2');
% sp_xid(5) = subplot(3,2,4); grid on;
% line(tSG, xidSG(5,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('wyd in RAD/s^2');
% sp_xid(6) = subplot(3,2,6); grid on;
% line(tSG, xidSG(6,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3, 'Linestyle', '-');
% ylabel('wzd in RAD/s^2');
% xlabel('t in s');
% 
% linkaxes(sp_xid, 'x');

%% LQ param identification
fprintf('Setting up identification equations...'); tic;

nE = 6;
nP = 19;

AA = zeros(nE*N, nP);
bb = zeros(nE*N, 1);
for k = 1:N
  R = Quaternion2RotMat(qSG(:,k));
  Rxx = R(1,1);
  Rxy = R(1,2);
  Rxz = R(1,3);
  Ryx = R(2,1);
  Ryy = R(2,2);
  Rzx = R(3,1);
  Rzy = R(3,2);
  Rzz = R(3,3);
  wx = wSG(1,k);
  wy = wSG(2,k);
  wz = wSG(3,k);
  wxd = wdSG(1,k);
  wyd = wdSG(2,k);
  wzd = wdSG(3,k);
  
  vx = xiSG(1,k);
  vy = xiSG(2,k);
  vz = xiSG(3,k);
  vxd = xidSG(1,k);
  vyd = xidSG(2,k);
  vzd = xidSG(3,k);
  
  wP1   = wP1SG(1,k);
  wP2   = wP2SG(1,k);
  wP3   = wP3SG(1,k);
  wP1d  = wP1SG(2,k);
  wP2d  = wP2SG(2,k);
  wP3d  = wP3SG(2,k);
  
  theta1    = pP1SG(1,k);
  theta2    = pP2SG(1,k);
  theta3    = pP3SG(1,k);
  theta1d   = pP1SG(2,k);
  theta2d   = pP2SG(2,k);
  theta3d   = pP3SG(2,k);
  theta1dd  = pP1SG(3,k);
  theta2dd  = pP2SG(3,k);
  theta3dd  = pP3SG(3,k);
  
  A = [0 0 0 -m * wy ^ 2 - m * wz ^ 2 m * wx * wy - m * wzd m * wx * wz + m * wyd vx 0 0 0 -sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 / 0.2e1 + sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 / 0.2e1 0 0 -1 0 0 0 0 0; 0 0 0 m * wx * wy + m * wzd -m * wx ^ 2 - m * wz ^ 2 m * wy * wz - m * wxd vy 0 0 0 sin(theta1) * wP1 ^ 2 / 0.2e1 - sin(theta2) * wP2 ^ 2 + sin(theta3) * wP3 ^ 2 / 0.2e1 0 0 0 -1 0 0 0 0; 0 0 0 m * wx * wz - m * wyd m * wy * wz + m * wxd -m * wx ^ 2 - m * wy ^ 2 0 vz 0 0 -cos(theta1) * wP1 ^ 2 - cos(theta2) * wP2 ^ 2 - cos(theta3) * wP3 ^ 2 0 0 0 0 -1 0 0 0; -wx * wz + wyd wx * wy + wzd wy ^ 2 - wz ^ 2 0 g * m * Rzz - m * vx * wy + m * vy * wx + m * vzd -g * m * Rzy - m * vx * wz + m * vz * wx - m * vyd 0 0 wx 0 -wP1 ^ 2 * (sqrt(0.3e1) * cos(theta1) * aA + sin(theta1) * hA) / 0.2e1 - wP3 ^ 2 * (-sqrt(0.3e1) * cos(theta3) * aA + sin(theta3) * hA) / 0.2e1 + sin(theta2) * hA * wP2 ^ 2 theta1dd / 0.2e1 - theta2dd + theta3dd / 0.2e1 + wz * sqrt(0.3e1) * theta3d / 0.2e1 - wz * sqrt(0.3e1) * theta1d / 0.2e1 -theta1d * cos(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 - theta3d * sqrt(0.3e1) * cos(theta3) * wP3 / 0.2e1 - sqrt(0.3e1) * sin(theta1) * wP1d / 0.2e1 - sqrt(0.3e1) * sin(theta3) * wP3d / 0.2e1 - wP1 ^ 2 * sqrt(0.3e1) * sin(theta1) * pP1 / 0.2e1 - wP3 ^ 2 * sqrt(0.3e1) * sin(theta3) * pP1 / 0.2e1 + wy * cos(theta3) * wP3 - wy * cos(theta1) * wP1 + wy * cos(theta2) * wP2 - wz * sin(theta1) * wP1 / 0.2e1 + wz * sin(theta3) * wP3 / 0.2e1 - wz * sin(theta2) * wP2 0 0 0 -1 0 0; wy * wz + wxd -wx ^ 2 + wz ^ 2 -wx * wy + wzd -g * m * Rzz + m * vx * wy - m * vy * wx - m * vzd 0 g * m * Rzx - m * vy * wz + m * vz * wy + m * vxd 0 0 wy 0 -wP1 ^ 2 * (sqrt(0.3e1) * sin(theta1) * hA - aA * cos(theta1)) / 0.2e1 - wP2 ^ 2 * aA * cos(theta2) + wP3 ^ 2 * (sin(theta3) * hA * sqrt(0.3e1) + aA * cos(theta3)) / 0.2e1 -sqrt(0.3e1) * theta3dd / 0.2e1 + sqrt(0.3e1) * theta1dd / 0.2e1 + wz * theta1d / 0.2e1 - wz * theta2d + wz * theta3d / 0.2e1 -wz * sin(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 - wz * sqrt(0.3e1) * sin(theta3) * wP3 / 0.2e1 - sin(theta3) * wP3d / 0.2e1 + sin(theta1) * wP1d / 0.2e1 + sin(theta2) * wP2d + wP1 ^ 2 * sin(theta1) * pP1 / 0.2e1 + wP2 ^ 2 * sin(theta2) * pP1 - wP3 ^ 2 * sin(theta3) * pP1 / 0.2e1 + wx * cos(theta1) * wP1 - wx * cos(theta2) * wP2 - wx * cos(theta3) * wP3 + theta2d * cos(theta2) * wP2 - theta3d * cos(theta3) * wP3 / 0.2e1 + theta1d * cos(theta1) * wP1 / 0.2e1 0 0 0 0 -1 0; wx ^ 2 - wy ^ 2 -wy * wz + wxd wx * wz + wyd g * m * Rzy + m * vx * wz - m * vz * wx + m * vyd -g * m * Rzx + m * vy * wz - m * vz * wy - m * vxd 0 0 0 0 wz wP1 ^ 2 * aA * sin(theta1) + wP2 ^ 2 * aA * sin(theta2) + wP3 ^ 2 * aA * sin(theta3) wx * sqrt(0.3e1) * theta1d / 0.2e1 - wx * sqrt(0.3e1) * theta3d / 0.2e1 - wy * theta1d / 0.2e1 + wy * theta2d - wy * theta3d / 0.2e1 wx * sin(theta1) * wP1 / 0.2e1 + wx * sin(theta2) * wP2 - wx * sin(theta3) * wP3 / 0.2e1 + theta1d * sin(theta1) * wP1 - theta2d * sin(theta2) * wP2 - theta3d * sin(theta3) * wP3 + wy * sin(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 + wy * sqrt(0.3e1) * sin(theta3) * wP3 / 0.2e1 + cos(theta2) * wP2d - cos(theta1) * wP1d + cos(theta3) * wP3d - wP1 ^ 2 * cos(theta1) * pP1 + wP2 ^ 2 * cos(theta2) * pP1 + wP3 ^ 2 * cos(theta3) * pP1 0 0 0 0 0 -1;];
  b = [-m * (Rzx * g - vy * wz + vz * wy + vxd); -m * (Rzy * g + vx * wz - vz * wx + vyd); -m * (Rzz * g - vx * wy + vy * wx + vzd); wy * wz * Jyy - wy * wz * Jzz - Jxx * wxd; -wx * wz * Jxx + wx * wz * Jzz - Jyy * wyd; wx * wy * Jxx - wx * wy * Jyy - Jzz * wzd;];

  AA(1+nE*(k-1):nE*k, :) = A;
  bb(1+nE*(k-1):nE*k, :) = b;
  
end
fprintf('%f\n', toc);

%%
fprintf('Solving... '); tic;
p = AA \ bb;
fprintf('%f\n', toc);

Jxy = p(1);
Jxz = p(2);
Jyz = p(3);
sx = p(4);
sy = p(5);
sz = p(6);
dvxy = p(7);
dvz = p(8);
dwxy = p(9);
dwz = p(10);
kappaF = p(11);
JAC = p(12);
JP = p(13);
FBx = p(14);
FBy = p(15);
FBz = p(16);
tauBx = p(17);
tauBy = p(18);
tauBz = p(19);

kappaT = pP1*JP;

fprintf('center of mass:        [ sx, sy, sz ] = [ %1.2f, %1.2f, %1.2f ] mm \n', sx*1e3, sy*1e3, sz*1e3);
fprintf('deviation inertia:     [ Jxy, Jxz, Jyz ] = [ %1.2e, %1.2e, %1.2e ] kg*m^2 \n', Jxy, Jxz, Jyz);
fprintf('translational damping: dvxy = %1.2e, dvz = %1.2e kg/s \n', dvxy, dvz);
fprintf('rotational damping:    dwxy = %1.2e, dwz = %1.2e kg*m^2/s \n', dwxy, dwz);
fprintf('bias force:            FB = [ %1.2e, %1.2e, %1.2e ] N, tauB = [ %1.2e, %1.2e, %1.2e ] Nm \n', FBx, FBy, FBz, tauBx, tauBy, tauBz);
fprintf('propeller param:       JAC = %1.3e, JP = %1.3e, kappaF = %1.3e, kappaT = %1.3e \n', JAC, JP, kappaF, kappaT);

%%
fprintf('Computing model accelerations... '); tic;

xidH1 = zeros(7, N);
for k = 1:N;
  R = Quaternion2RotMat(qSG(:,k));
  Rxx = R(1,1);
  Rxy = R(1,2);
  Rxz = R(1,3);
  Ryx = R(2,1);
  Ryy = R(2,2);
  Rzx = R(3,1);
  Rzy = R(3,2);
  Rzz = R(3,3);
  wx = wSG(1,k);
  wy = wSG(2,k);
  wz = wSG(3,k);
  wxd = wdSG(1,k);
  wyd = wdSG(2,k);
  wzd = wdSG(3,k);
  
  vx = xiSG(1,k);
  vy = xiSG(2,k);
  vz = xiSG(3,k);
  vxd = xidSG(1,k);
  vyd = xidSG(2,k);
  vzd = xidSG(3,k);
  
  wP1   = wP1SG(1,k);
  wP2   = wP2SG(1,k);
  wP3   = wP3SG(1,k);
  wP1d  = wP1SG(2,k);
  wP2d  = wP2SG(2,k);
  wP3d  = wP3SG(2,k);
  
  theta1    = pP1SG(1,k);
  theta2    = pP2SG(1,k);
  theta3    = pP3SG(1,k);
  theta1d   = pP1SG(2,k);
  theta2d   = pP2SG(2,k);
  theta3d   = pP3SG(2,k);
  theta1dd  = pP1SG(3,k);
  theta2dd  = pP2SG(3,k);
  theta3dd  = pP3SG(3,k);
  
  Mb = [m 0 0 0 m * sz -m * sy; 0 m 0 -m * sz 0 m * sx; 0 0 m m * sy -m * sx 0; 0 -m * sz m * sy Jxx Jxy Jxz; m * sz 0 -m * sx Jxy Jyy Jyz; -m * sy m * sx 0 Jxz Jyz Jzz;];
  fb = [m * sx * wy ^ 2 + m * sx * wz ^ 2 - m * sy * wx * wy - m * sz * wx * wz - g * m * Rzx + m * vy * wz - m * vz * wy - dvxy * vx + sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 * kappaF / 0.2e1 - sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 * kappaF / 0.2e1 + FBx; -m * sx * wx * wy + m * sy * wx ^ 2 + m * sy * wz ^ 2 - m * sz * wy * wz - g * m * Rzy - m * vx * wz + m * vz * wx - dvxy * vy - sin(theta1) * wP1 ^ 2 * kappaF / 0.2e1 + sin(theta2) * wP2 ^ 2 * kappaF - sin(theta3) * wP3 ^ 2 * kappaF / 0.2e1 + FBy; -m * sx * wx * wz - m * sy * wy * wz + m * sz * wx ^ 2 + m * sz * wy ^ 2 - g * m * Rzz + m * vx * wy - m * vy * wx - dvz * vz + cos(theta1) * wP1 ^ 2 * kappaF + cos(theta2) * wP2 ^ 2 * kappaF + cos(theta3) * wP3 ^ 2 * kappaF + FBz; theta1d * JP * cos(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 + theta3d * JP * sqrt(0.3e1) * cos(theta3) * wP3 / 0.2e1 + kappaF * cos(theta1) * sqrt(0.3e1) * aA * wP1 ^ 2 / 0.2e1 - kappaF * cos(theta3) * sqrt(0.3e1) * aA * wP3 ^ 2 / 0.2e1 + JP * wP1 ^ 2 * sqrt(0.3e1) * sin(theta1) * pP1 / 0.2e1 + JP * wP3 ^ 2 * sqrt(0.3e1) * sin(theta3) * pP1 / 0.2e1 + kappaF * sin(theta1) * hA * wP1 ^ 2 / 0.2e1 + kappaF * sin(theta3) * hA * wP3 ^ 2 / 0.2e1 - JAC * theta1dd / 0.2e1 + JAC * theta2dd - JAC * theta3dd / 0.2e1 - dwxy * wx + sqrt(0.3e1) * sin(theta1) * JP * wP1d / 0.2e1 + sqrt(0.3e1) * sin(theta3) * JP * wP3d / 0.2e1 - sin(theta2) * hA * wP2 ^ 2 * kappaF - wy ^ 2 * Jyz + wz ^ 2 * Jyz + vx * wy * m * sy - vy * wx * m * sy + vx * wz * m * sz - vz * m * sz * wx + Rzy * g * m * sz - Rzz * g * m * sy - wy * JP * cos(theta3) * wP3 + wy * JP * cos(theta1) * wP1 - wy * JP * cos(theta2) * wP2 + wz * JP * sin(theta1) * wP1 / 0.2e1 - wz * JP * sin(theta3) * wP3 / 0.2e1 - wz * JAC * sqrt(0.3e1) * theta3d / 0.2e1 + wz * JP * sin(theta2) * wP2 + wz * JAC * sqrt(0.3e1) * theta1d / 0.2e1 + wy * wz * Jyy + wz * wx * Jxy - wy * wx * Jxz - wy * wz * Jzz + tauBx; -kappaF * cos(theta3) * aA * wP3 ^ 2 / 0.2e1 - JP * wP1 ^ 2 * sin(theta1) * pP1 / 0.2e1 - JP * wP2 ^ 2 * sin(theta2) * pP1 + JP * wP3 ^ 2 * sin(theta3) * pP1 / 0.2e1 - kappaF * cos(theta1) * aA * wP1 ^ 2 / 0.2e1 + kappaF * wP2 ^ 2 * aA * cos(theta2) + wz * JP * sin(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 + wz * JP * sqrt(0.3e1) * sin(theta3) * wP3 / 0.2e1 + vy * wx * m * sx + vy * wz * m * sz - vx * wy * m * sx - vz * wy * m * sz - dwxy * wy + kappaF * sqrt(0.3e1) * sin(theta1) * hA * wP1 ^ 2 / 0.2e1 - kappaF * sin(theta3) * sqrt(0.3e1) * hA * wP3 ^ 2 / 0.2e1 + wx ^ 2 * Jxz - wz ^ 2 * Jxz - wx * wz * Jxx + wx * wy * Jyz + wx * wz * Jzz - wz * wy * Jxy - wz * JAC * theta1d / 0.2e1 + wz * JAC * theta2d - wz * JAC * theta3d / 0.2e1 - wx * JP * cos(theta1) * wP1 + wx * JP * cos(theta2) * wP2 + wx * JP * cos(theta3) * wP3 - theta2d * JP * cos(theta2) * wP2 + theta3d * JP * cos(theta3) * wP3 / 0.2e1 - theta1d * JP * cos(theta1) * wP1 / 0.2e1 + sqrt(0.3e1) * JAC * theta3dd / 0.2e1 - sqrt(0.3e1) * JAC * theta1dd / 0.2e1 - sin(theta1) * JP * wP1d / 0.2e1 - sin(theta2) * JP * wP2d + sin(theta3) * JP * wP3d / 0.2e1 - Rzx * g * m * sz + Rzz * g * m * sx + tauBy; -vx * m * sx * wz + vz * wx * m * sx + vz * wy * m * sy - vy * wz * m * sy - wx * JAC * sqrt(0.3e1) * theta1d / 0.2e1 + wx * JAC * sqrt(0.3e1) * theta3d / 0.2e1 - wx * JP * sin(theta1) * wP1 / 0.2e1 - wx * JP * sin(theta2) * wP2 + wx * JP * sin(theta3) * wP3 / 0.2e1 - theta1d * JP * sin(theta1) * wP1 + theta2d * JP * sin(theta2) * wP2 + theta3d * JP * sin(theta3) * wP3 + wx * wy * Jxx - wx * wy * Jyy - wx * wz * Jyz + wy * wz * Jxz + wy * JAC * theta1d / 0.2e1 - wy * JAC * theta2d + wy * JAC * theta3d / 0.2e1 - wy * JP * sin(theta1) * sqrt(0.3e1) * wP1 / 0.2e1 - wy * JP * sqrt(0.3e1) * sin(theta3) * wP3 / 0.2e1 + JP * wP1 ^ 2 * cos(theta1) * pP1 - JP * wP2 ^ 2 * cos(theta2) * pP1 - JP * wP3 ^ 2 * cos(theta3) * pP1 - sin(theta3) * aA * kappaF * wP3 ^ 2 - sin(theta1) * aA * kappaF * wP1 ^ 2 - sin(theta2) * aA * kappaF * wP2 ^ 2 - wx ^ 2 * Jxy + wy ^ 2 * Jxy - dwz * wz + Rzx * g * m * sy - Rzy * g * m * sx + cos(theta1) * JP * wP1d - cos(theta2) * JP * wP2d - cos(theta3) * JP * wP3d + tauBz;];

  xidH1(1:6,k) = Mb \ fb;
  
  xidH1(7,k) = pP2*iM1SG(1,k) - pP3 - pP1*wP1^2;
end
fprintf('%f\n', toc);

%% LQ param identification
fprintf('Setting up identification equations...'); tic;

nE = 6;
nP = 9;

AA = zeros(nE*N, nP);
bb = zeros(nE*N, 1);
for k = 1:N
  R = Quaternion2RotMat(qSG(:,k));
  Rxx = R(1,1);
  Rxy = R(1,2);
  Rxz = R(1,3);
  Ryx = R(2,1);
  Ryy = R(2,2);
  Rzx = R(3,1);
  Rzy = R(3,2);
  Rzz = R(3,3);
  wx = wSG(1,k);
  wy = wSG(2,k);
  wz = wSG(3,k);
  wxd = wdSG(1,k);
  wyd = wdSG(2,k);
  wzd = wdSG(3,k);
  
  vx = xiSG(1,k);
  vy = xiSG(2,k);
  vz = xiSG(3,k);
  vxd = xidSG(1,k);
  vyd = xidSG(2,k);
  vzd = xidSG(3,k);
  
  wP1   = wP1SG(1,k);
  wP2   = wP2SG(1,k);
  wP3   = wP3SG(1,k);
  wP1d  = wP1SG(2,k);
  wP2d  = wP2SG(2,k);
  wP3d  = wP3SG(2,k);
  
  theta1    = pP1SG(1,k);
  theta2    = pP2SG(1,k);
  theta3    = pP3SG(1,k);
  theta1d   = pP1SG(2,k);
  theta2d   = pP2SG(2,k);
  theta3d   = pP3SG(2,k);
  theta1dd  = pP1SG(3,k);
  theta2dd  = pP2SG(3,k);
  theta3dd  = pP3SG(3,k);
  
  A = [vx -sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 / 0.2e1 + sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 / 0.2e1 0 -1 0 0 0 0 0; vy sin(theta1) * wP1 ^ 2 / 0.2e1 - sin(theta2) * wP2 ^ 2 + sin(theta3) * wP3 ^ 2 / 0.2e1 0 0 -1 0 0 0 0; vz -cos(theta1) * wP1 ^ 2 - cos(theta2) * wP2 ^ 2 - cos(theta3) * wP3 ^ 2 0 0 0 -1 0 0 0; 0 -wP1 ^ 2 * (sqrt(0.3e1) * cos(theta1) * aA + sin(theta1) * hA) / 0.2e1 - wP3 ^ 2 * (-sqrt(0.3e1) * cos(theta3) * aA + sin(theta3) * hA) / 0.2e1 + sin(theta2) * hA * wP2 ^ 2 -sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 / 0.2e1 - sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 / 0.2e1 0 0 0 -1 0 0; 0 -wP1 ^ 2 * (sqrt(0.3e1) * sin(theta1) * hA - aA * cos(theta1)) / 0.2e1 - wP2 ^ 2 * aA * cos(theta2) + wP3 ^ 2 * (sin(theta3) * hA * sqrt(0.3e1) + aA * cos(theta3)) / 0.2e1 sin(theta2) * wP2 ^ 2 + sin(theta1) * wP1 ^ 2 / 0.2e1 - sin(theta3) * wP3 ^ 2 / 0.2e1 0 0 0 0 -1 0; 0 wP1 ^ 2 * aA * sin(theta1) + wP2 ^ 2 * aA * sin(theta2) + wP3 ^ 2 * aA * sin(theta3) cos(theta2) * wP2 ^ 2 - cos(theta1) * wP1 ^ 2 + cos(theta3) * wP3 ^ 2 0 0 0 0 0 -1;];
  b = [-m * (Rzx * g - vy * wz + vz * wy + vxd); -m * (Rzy * g + vx * wz - vz * wx + vyd); -m * (Rzz * g - vx * wy + vy * wx + vzd); wy * wz * Jyy - wy * wz * Jzz - Jxx * wxd; -wx * wz * Jxx + wx * wz * Jzz - Jyy * wyd; wx * wy * Jxx - wx * wy * Jyy - Jzz * wzd;];

  AA(1+nE*(k-1):nE*k, :) = A;
  bb(1+nE*(k-1):nE*k, :) = b;
  
end
fprintf('%f\n', toc);

%%
fprintf('Solving... '); tic;
p = AA \ bb;
fprintf('%f\n', toc);

dv = p(1);
kappaF = p(2);
kappaT = p(3);
FBx = p(4);
FBy = p(5);
FBz = p(6);
tauBx = p(7);
tauBy = p(8);
tauBz = p(9);
Jxy = 0;
Jxz = 0;
Jyz = 0;
sx = 0;
sy = 0;
sz = 0;
JP = 0;
JAC = 0;
dvxy = dv;
dvz = dv;
dwxy = 0;
dwz = 0;

fprintf('center of mass:        [ sx, sy, sz ] = [ %1.2f, %1.2f, %1.2f ] mm \n', sx*1e3, sy*1e3, sz*1e3);
fprintf('deviation inertia:     [ Jxy, Jxz, Jyz ] = [ %1.2e, %1.2e, %1.2e ] kg*m^2 \n', Jxy, Jxz, Jyz);
fprintf('translational damping: dvxy = %1.2e, dvz = %1.2e kg/s \n', dvxy, dvz);
fprintf('rotational damping:    dwxy = %1.2e, dwz = %1.2e kg*m^2/s \n', dwxy, dwz);
fprintf('bias force:            FB = [ %1.2e, %1.2e, %1.2e ] N, tauB = [ %1.2e, %1.2e, %1.2e ] Nm \n', FBx, FBy, FBz, tauBx, tauBy, tauBz);
fprintf('propeller param:       JAC = %1.3e, JP = %1.3e, kappaF = %1.3e, kappaT = %1.3e \n', JAC, JP, kappaF, kappaT);

%%
fprintf('Computing model accelerations... '); tic;

xidH2 = zeros(6, N);
for k = 1:N;
  R = Quaternion2RotMat(qSG(:,k));
  Rxx = R(1,1);
  Rxy = R(1,2);
  Rxz = R(1,3);
  Ryx = R(2,1);
  Ryy = R(2,2);
  Rzx = R(3,1);
  Rzy = R(3,2);
  Rzz = R(3,3);
  wx = wSG(1,k);
  wy = wSG(2,k);
  wz = wSG(3,k);
  wxd = wdSG(1,k);
  wyd = wdSG(2,k);
  wzd = wdSG(3,k);
  
  vx = xiSG(1,k);
  vy = xiSG(2,k);
  vz = xiSG(3,k);
  vxd = xidSG(1,k);
  vyd = xidSG(2,k);
  vzd = xidSG(3,k);
  
  wP1   = wP1SG(1,k);
  wP2   = wP2SG(1,k);
  wP3   = wP3SG(1,k);
  wP1d  = wP1SG(2,k);
  wP2d  = wP2SG(2,k);
  wP3d  = wP3SG(2,k);
  
  theta1    = pP1SG(1,k);
  theta2    = pP2SG(1,k);
  theta3    = pP3SG(1,k);
  theta1d   = pP1SG(2,k);
  theta2d   = pP2SG(2,k);
  theta3d   = pP3SG(2,k);
  theta1dd  = pP1SG(3,k);
  theta2dd  = pP2SG(3,k);
  theta3dd  = pP3SG(3,k);
  
  Mb = [m 0 0 0 m * sz -m * sy; 0 m 0 -m * sz 0 m * sx; 0 0 m m * sy -m * sx 0; 0 -m * sz m * sy Jxx Jxy Jxz; m * sz 0 -m * sx Jxy Jyy Jyz; -m * sy m * sx 0 Jxz Jyz Jzz;];
  fb = [-m * sy * wzd + m * sz * wyd - g * m * Rzx + m * vy * wz - m * vz * wy - dv * vx + sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 * kappaF / 0.2e1 - sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 * kappaF / 0.2e1 + FBx; m * sx * wzd - m * sz * wxd - g * m * Rzy - m * vx * wz + m * vz * wx - dv * vy - sin(theta1) * wP1 ^ 2 * kappaF / 0.2e1 + sin(theta2) * wP2 ^ 2 * kappaF - sin(theta3) * wP3 ^ 2 * kappaF / 0.2e1 + FBy; -m * sx * wyd + m * sy * wxd - g * m * Rzz + m * vx * wy - m * vy * wx - dv * vz + cos(theta1) * wP1 ^ 2 * kappaF + cos(theta2) * wP2 ^ 2 * kappaF + cos(theta3) * wP3 ^ 2 * kappaF + FBz; m * sy * vzd - m * sz * vyd + Jxy * wyd + Jxz * wzd - kappaF * cos(theta3) * sqrt(0.3e1) * aA * wP3 ^ 2 / 0.2e1 + kappaT * sqrt(0.3e1) * sin(theta3) * wP3 ^ 2 / 0.2e1 + kappaF * sin(theta3) * hA * wP3 ^ 2 / 0.2e1 + kappaF * cos(theta1) * sqrt(0.3e1) * aA * wP1 ^ 2 / 0.2e1 + kappaT * sqrt(0.3e1) * sin(theta1) * wP1 ^ 2 / 0.2e1 + kappaF * sin(theta1) * hA * wP1 ^ 2 / 0.2e1 - sin(theta2) * hA * wP2 ^ 2 * kappaF + wy * wz * Jyy - wy * wz * Jzz + tauBx; -m * sx * vzd + m * sz * vxd + Jxy * wxd + Jyz * wzd + kappaF * wP2 ^ 2 * aA * cos(theta2) - kappaT * sin(theta2) * wP2 ^ 2 + kappaF * sqrt(0.3e1) * sin(theta1) * hA * wP1 ^ 2 / 0.2e1 - kappaF * cos(theta1) * aA * wP1 ^ 2 / 0.2e1 - kappaT * sin(theta1) * wP1 ^ 2 / 0.2e1 - kappaF * sin(theta3) * sqrt(0.3e1) * hA * wP3 ^ 2 / 0.2e1 - kappaF * cos(theta3) * aA * wP3 ^ 2 / 0.2e1 + kappaT * sin(theta3) * wP3 ^ 2 / 0.2e1 - wx * wz * Jxx + wx * wz * Jzz + tauBy; m * sx * vyd - m * sy * vxd + Jxz * wxd + Jyz * wyd - sin(theta2) * aA * kappaF * wP2 ^ 2 - cos(theta2) * kappaT * wP2 ^ 2 - sin(theta1) * aA * kappaF * wP1 ^ 2 + cos(theta1) * kappaT * wP1 ^ 2 + wx * wy * Jxx - wx * wy * Jyy + tauBz - sin(theta3) * aA * kappaF * wP3 ^ 2 - cos(theta3) * kappaT * wP3 ^ 2;];

  xidH2(:,k) = Mb \ fb;
end
fprintf('%f\n', toc);

%%
ColorMeas = 0.7*[1 1 1];
LineWidthMeas = 3;
LineStyleMeas = '-';

LineWidthIdent = 1;

ColorIdent1 = 'b';
ColorIdent2 = 'r';
ColorIdent3 = 'm';
LineStyleIdent1 = '-';
LineStyleIdent2 = '--';
LineStyleIdent3 = '--';

fig = figure(1); clf;
sp_xi(1) = subplot(5,1,1); grid on; hold on;
line(tSG, xidSG(1,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(1,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH2(1,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
ylabel('vxd in m/s^2', 'Interpreter', 'none');

sp_xi(2) = subplot(5,1,2); grid on; hold on;
line(tSG, xidSG(3,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(3,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH2(3,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
ylabel('vzd in m/s^2', 'Interpreter', 'none');

sp_xi(3) = subplot(5,1,3); grid on; hold on;
line(tSG, xidSG(5,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(5,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH2(5,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
ylim(12*[-1 1]);
ylabel('wyd in RAD/s^2', 'Interpreter', 'none');

sp_xi(4) = subplot(5,1,4); grid on; hold on;
line(tSG, xidSG(6,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(6,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH2(6,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
%line(tSG, xidH3(6,:), 'Color', ColorIdent3, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent3);
ylim(3.5*[-1 1]);
ylabel('wzd in RAD/s2', 'Interpreter', 'none');

sp_xi(5) = subplot(5,1,5); grid on; hold on;
line(tSG, wP1SG(2,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(7,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH1(7,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
ylim(500*[-1 1]);
ylabel('wP1d in RAD/s^2', 'Interpreter', 'none');
xlabel('time in s', 'Interpreter', 'none');

linkaxes(sp_xi, 'x');

axes(sp_xi(1));
%l_xi = legend('measured', 'full model', 'simp. model', 'oversimp. model', 'location', 'northeast');
l_xi = legend('measured', 'full model', 'simp. model', 'location', 'northeast');

%%
if (saveFigure > 0)
  
  sp_xi(1).UserData.LatexYLabel = '$\vxd$ in $\sfrac{\unit{m}}{\unit{s}^2}$';
  sp_xi(2).UserData.LatexYLabel = '$\vzd$ in $\sfrac{\unit{m}}{\unit{s}^2}$';
  sp_xi(3).UserData.LatexYLabel = '$\wyd$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
  sp_xi(4).UserData.LatexYLabel = '$\wzd$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
  sp_xi(5).UserData.LatexYLabel = '$\PropVeld[1]$ in $\sfrac{\unit{RAD}}{\unit{s}^2}$';
  sp_xi(5).UserData.LatexXLabel = '$t$ in $\unit{s}$';
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'ModelValidationTri', 'pdf', 0 );
  
end