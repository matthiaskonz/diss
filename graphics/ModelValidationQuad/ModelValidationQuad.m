clear;

filename = '2018-08-30_QuadV4.mat';
tS = 70;
tF = 90;

saveFigure = 1;
figSize = [15.9 23.3];

%% param
Ts = 5e-3;

m = 1.001;           % [kg] QuadV4 mass
Jxx = 17.9e-3;
Jyy = 18.0e-3;
Jzz = 30.7e-3;
aA = 0.24;
g = 9.81;

% identified propeller param
pP1 = 6.855e-3;
pP2 = 5.154;
pP3 = 1.876e2;

%%
load(filename);
meas = meas{1};

kk = find((meas.vSerial(1,:) >= tS) & (meas.vSerial(1,:) <= tF));
tM = meas.vSerial(1,kk);
t = [tM(5) : Ts : tM(end-5)];
q = interp1(tM', meas.vSerial([21; 2; 4; 6],kk)', t')';
w = interp1(tM', meas.vSerial([3 5 7],kk)', t'+Ts)';
rd = interp1(tM', meas.vSerial([12 14 16],kk)', t')';
a = interp1(tM', meas.vSerial([13 15 17],kk)', t'+Ts)';
wP = interp1(tM', meas.vSerial([8:11],kk)', t'+1.5*Ts)';
iM = interp1(tM', meas.vSerial([18:20],kk)', t')';
iM = [iM; zeros(size(t))];

t = t-t(1);

%% SG filter and derivative estimation
SG_ORDER = 5;
SG_WINDOW = 51;
[tSG, q0SG] = sgolayDiffEstim(t, q(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q1SG] = sgolayDiffEstim(t, q(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q2SG] = sgolayDiffEstim(t, q(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, q3SG] = sgolayDiffEstim(t, q(4,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, w1SG] = sgolayDiffEstim(t, w(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, w2SG] = sgolayDiffEstim(t, w(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, w3SG] = sgolayDiffEstim(t, w(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r1dSG] = sgolayDiffEstim(t, rd(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r2dSG] = sgolayDiffEstim(t, rd(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, r3dSG] = sgolayDiffEstim(t, rd(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP1SG] = sgolayDiffEstim(t, wP(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP2SG] = sgolayDiffEstim(t, wP(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP3SG] = sgolayDiffEstim(t, wP(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, wP4SG] = sgolayDiffEstim(t, wP(4,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, iM1SG] = sgolayDiffEstim(t, iM(1,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, iM2SG] = sgolayDiffEstim(t, iM(2,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, iM3SG] = sgolayDiffEstim(t, iM(3,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);
[~, iM4SG] = sgolayDiffEstim(t, iM(4,:), SG_ORDER, SG_WINDOW, Ts, 2, 0);

N = length(tSG); % number of measurement points

qSG = [q0SG(1,:); q1SG(1,:); q2SG(1,:); q3SG(1,:)];
for k = 1:N
  qSG(:,k) = qSG(:,k) / norm(qSG(:,k));
end

wSG = [w1SG(1,:); w2SG(1,:); w3SG(1,:)];
wdSG = [w1SG(2,:); w2SG(2,:); w3SG(2,:)];

rdSG = [r1dSG(1,:); r2dSG(1,:); r3dSG(1,:)];
rddSG = [r1dSG(2,:); r2dSG(2,:); r3dSG(2,:)];

wPSG = [wP1SG(1,:); wP2SG(1,:); wP3SG(1,:); wP4SG(1,:)];
wPdSG = [wP1SG(2,:); wP2SG(2,:); wP3SG(2,:); wP4SG(2,:)];

%iMSG = interp1(t', iM', tSG')';
iMSG = [iM1SG(1,:); iM2SG(1,:); iM3SG(1,:); iM4SG(1,:)];

% aSG = [a1SG(1,:); a2SG(1,:); a3SG(1,:)];
% 
xiSG = zeros(10, N);
xidSG = zeros(10, N);
%vdIMU = zeros(3, N);
for k = 1:N
  R = Quaternion2RotMat(qSG(:,k));
  xiSG(:,k) = [R'*rdSG(:,k); wSG(:,k); wPSG(:,k)];
  xidSG(:,k) = [R'*rddSG(:,k) - wed_so3(xiSG(4:6,k))*xiSG(1:3,k); wdSG(:,k); wPdSG(:,k)];
%  vdIMU(:,k) = aSG(:,k) + R'*[0;0;-g] - wed_so3(xiSG(4:6,k))*xiSG(1:3,k);
end

xidN = zeros(10,N-1);
for k = 1:N-1
  xidN(:,k) = 1/Ts*(xiSG(:,k+1)-xiSG(:,k));
end

%% LQ param identification
fprintf('Setting up identification equations...'); tic;

nE = 6;
nP = 18;

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
  
  wP1 = wPSG(1,k);
  wP2 = wPSG(2,k);
  wP3 = wPSG(3,k);
  wP4 = wPSG(4,k);
  wP1d = wPdSG(1,k);
  wP2d = wPdSG(2,k);
  wP3d = wPdSG(3,k);
  wP4d = wPdSG(4,k);
  
  A = [-m * wy ^ 2 - m * wz ^ 2 m * wx * wy - m * wzd m * wx * wz + m * wyd 0 0 0 vx 0 0 0 0 0 -1 0 0 0 0 0; m * wx * wy + m * wzd -m * wx ^ 2 - m * wz ^ 2 m * wy * wz - m * wxd 0 0 0 vy 0 0 0 0 0 0 -1 0 0 0 0; m * wx * wz - m * wyd m * wy * wz + m * wxd -m * wx ^ 2 - m * wy ^ 2 0 0 0 0 vz 0 0 -wP1 ^ 2 - wP2 ^ 2 - wP3 ^ 2 - wP4 ^ 2 0 0 0 -1 0 0 0; 0 g * m * Rzz - m * vx * wy + m * vy * wx + m * vzd -g * m * Rzy - m * vx * wz + m * vz * wx - m * vyd -wx * wz + wyd wx * wy + wzd wy ^ 2 - wz ^ 2 0 0 wx 0 -aA * wP2 ^ 2 + aA * wP4 ^ 2 wP1 * wy - wP2 * wy + wP3 * wy - wP4 * wy 0 0 0 -1 0 0; -g * m * Rzz + m * vx * wy - m * vy * wx - m * vzd 0 g * m * Rzx - m * vy * wz + m * vz * wy + m * vxd wy * wz + wxd -wx ^ 2 + wz ^ 2 -wx * wy + wzd 0 0 wy 0 aA * wP1 ^ 2 - aA * wP3 ^ 2 -wP1 * wx + wP2 * wx - wP3 * wx + wP4 * wx 0 0 0 0 -1 0; g * m * Rzy + m * vx * wz - m * vz * wx + m * vyd -g * m * Rzx + m * vy * wz - m * vz * wy - m * vxd 0 wx ^ 2 - wy ^ 2 -wy * wz + wxd wx * wz + wyd 0 0 0 wz 0 wP1 ^ 2 * pP1 - wP2 ^ 2 * pP1 + wP3 ^ 2 * pP1 - wP4 ^ 2 * pP1 + wP1d - wP2d + wP3d - wP4d 0 0 0 0 0 -1;];
  b = [-g * m * Rzx + m * vy * wz - m * vz * wy - m * vxd; -g * m * Rzy - m * vx * wz + m * vz * wx - m * vyd; -g * m * Rzz + m * vx * wy - m * vy * wx - m * vzd; Jyy * wy * wz - Jzz * wy * wz - Jxx * wxd; -Jxx * wx * wz + Jzz * wx * wz - Jyy * wyd; Jxx * wx * wy - Jyy * wx * wy - Jzz * wzd;];

  AA(1+nE*(k-1):nE*k, :) = A;
  bb(1+nE*(k-1):nE*k, :) = b;
  
end
fprintf('%f\n', toc);

%%
fprintf('Solving... '); tic;
p = AA \ bb;
fprintf('%f\n', toc);

sx = p(1);
sy = p(2);
sz = p(3);
Jxy = p(4);
Jxz = p(5);
Jyz = p(6);
dvxy = p(7);
dvz = p(8);
dwxy = p(9);
dwz = p(10);
kappaF = p(11);
JP = p(12);
FBx = p(13);
FBy = p(14);
FBz = p(15);
tauBx = p(16);
tauBy = p(17);
tauBz = p(18);

kappaT = pP1*JP;

fprintf('center of mass:        [ sx, sy, sz ] = [ %1.2f, %1.2f, %1.2f ] mm \n', sx*1e3, sy*1e3, sz*1e3);
fprintf('deviation inertia:     [ Jxy, Jxz, Jyz ] = [ %1.2e, %1.2e, %1.2e ] kg*m^2 \n', Jxy, Jxz, Jyz);
fprintf('translational damping: dvxy = %1.2e, dvz = %1.2e kg/s \n', dvxy, dvz);
fprintf('rotational damping:    dwxy = %1.2e, dwz = %1.2e kg*m^2/s \n', dwxy, dwz);
fprintf('bias force:            FB = [ %1.2e, %1.2e, %1.2e ] N, tauB = [ %1.2e, %1.2e, %1.2e ] Nm \n', FBx, FBy, FBz, tauBx, tauBy, tauBz);
fprintf('propeller param:       JP = %1.3e , kappaF = %1.3e, kappaT = %1.3e \n', JP, kappaF, kappaT);
%fprintf('motor param:           kM = [%2.1f, %2.1f, %2.1f, %2.1f] * 1e-6, tau0 = [%2.1f, %2.1f, %2.1f %2.1f] * 1e-3 \n', kM1*1e6, kM2*1e6, kM3*1e6, kM4*1e6, tauM01*1e3, tauM02*1e3, tauM03*1e3, tauM04*1e3);

%%
fprintf('Computing model accelerations... '); tic;

xidH1 = zeros(10, N);
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
  
  wP1 = wPSG(1,k);
  wP2 = wPSG(2,k);
  wP3 = wPSG(3,k);
  wP4 = wPSG(4,k);
  wP1d = wPdSG(1,k);
  wP2d = wPdSG(2,k);
  wP3d = wPdSG(3,k);
  wP4d = wPdSG(4,k);
  iM1 = iMSG(1,k);
  iM2 = iMSG(2,k);
  iM3 = iMSG(3,k);
  iM4 = iMSG(4,k);
  
  Mb = [m 0 0 0 m * sz -m * sy; 0 m 0 -m * sz 0 m * sx; 0 0 m m * sy -m * sx 0; 0 -m * sz m * sy Jxx Jxy Jxz; m * sz 0 -m * sx Jxy Jyy Jyz; -m * sy m * sx 0 Jxz Jyz Jzz;];
  fb = [m * sx * wy ^ 2 + m * sx * wz ^ 2 - m * sy * wx * wy - m * sz * wx * wz - g * m * Rzx + m * vy * wz - m * vz * wy - dvxy * vx + FBx; -m * sx * wx * wy + m * sy * wx ^ 2 + m * sy * wz ^ 2 - m * sz * wy * wz - g * m * Rzy - m * vx * wz + m * vz * wx - dvxy * vy + FBy; -m * sx * wx * wz - m * sy * wy * wz + m * sz * wx ^ 2 + m * sz * wy ^ 2 - g * m * Rzz + kappaF * wP1 ^ 2 + kappaF * wP2 ^ 2 + kappaF * wP3 ^ 2 + kappaF * wP4 ^ 2 + m * vx * wy - m * vy * wx - dvz * vz + FBz; -dwxy * wx - JP * wP1 * wy - Jxz * wx * wy + Jyy * wy * wz + JP * wP2 * wy + JP * wP4 * wy - Jzz * wy * wz + Jxy * wx * wz - JP * wP3 * wy + m * sy * vx * wy - m * sy * vy * wx + m * sz * vx * wz - m * sz * vz * wx + Rzy * g * m * sz - Rzz * g * m * sy - Jyz * wy ^ 2 + Jyz * wz ^ 2 + tauBx + aA * kappaF * wP2 ^ 2 - aA * kappaF * wP4 ^ 2; -dwxy * wy + JP * wP3 * wx + JP * wP1 * wx - JP * wP4 * wx - JP * wP2 * wx - Jxx * wx * wz - Jxy * wy * wz + Jyz * wx * wy + Jzz * wx * wz - m * sx * vx * wy + m * sx * vy * wx + m * sz * vy * wz - m * sz * vz * wy - Rzx * g * m * sz + Rzz * g * m * sx + Jxz * wx ^ 2 - Jxz * wz ^ 2 + tauBy - aA * kappaF * wP1 ^ 2 + aA * kappaF * wP3 ^ 2; -dwz * wz + wP4 ^ 2 * pP1 * JP + wP2 ^ 2 * pP1 * JP - wP3 ^ 2 * pP1 * JP - wP1 ^ 2 * pP1 * JP + Jxx * wx * wy + Jxz * wy * wz - m * sx * vx * wz + m * sx * vz * wx - m * sy * vy * wz + m * sy * vz * wy + Rzx * g * m * sy - Rzy * g * m * sx - Jxy * wx ^ 2 + Jxy * wy ^ 2 - JP * wP1d - JP * wP3d + JP * wP2d + JP * wP4d + tauBz - Jyy * wx * wy - Jyz * wx * wz;];

  xidH1(1:6,k) = Mb \ fb;
  
  xidH1(7,k)  = pP2*iM1 - pP3 - pP1*wP1^2;
  xidH1(8,k)  = pP2*iM2 - pP3 - pP1*wP2^2;
  xidH1(9,k)  = pP2*iM3 - pP3 - pP1*wP3^2;
  xidH1(10,k) = pP2*iM4 - pP3 - pP1*wP4^2;
end
fprintf('%f\n', toc);


%% LQ param identification
fprintf('Setting up identification equations...'); tic;

nE = 6;
nP = 10;

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
  
  wP1 = wPSG(1,k);
  wP2 = wPSG(2,k);
  wP3 = wPSG(3,k);
  wP4 = wPSG(4,k);
  wP1d = wPdSG(1,k);
  wP2d = wPdSG(2,k);
  wP3d = wPdSG(3,k);
  wP4d = wPdSG(4,k);
  
  A = [vx 0 0 0 -1 0 0 0 0 0; vy 0 0 0 0 -1 0 0 0 0; vz 0 -wP1 ^ 2 - wP2 ^ 2 - wP3 ^ 2 - wP4 ^ 2 0 0 0 -1 0 0 0; 0 0 -aA * wP2 ^ 2 + aA * wP4 ^ 2 0 0 0 0 -1 0 0; 0 0 aA * wP1 ^ 2 - aA * wP3 ^ 2 0 0 0 0 0 -1 0; 0 wz 0 pP1 * wP1 ^ 2 - pP1 * wP2 ^ 2 + pP1 * wP3 ^ 2 - pP1 * wP4 ^ 2 + wP1d - wP2d + wP3d - wP4d 0 0 0 0 0 -1;];
  b = [-g * m * Rzx + m * vy * wz - m * vz * wy - m * vxd; -g * m * Rzy - m * vx * wz + m * vz * wx - m * vyd; -g * m * Rzz + m * vx * wy - m * vy * wx - m * vzd; Jyy * wy * wz - Jzz * wy * wz - Jxx * wxd; -Jxx * wx * wz + Jzz * wx * wz - Jyy * wyd; Jxx * wx * wy - Jyy * wx * wy - Jzz * wzd;];

  AA(1+nE*(k-1):nE*k, :) = A;
  bb(1+nE*(k-1):nE*k, :) = b;
  
end
fprintf('%f s\n', toc);

%%
fprintf('Solving... '); tic;
p = AA \ bb;
fprintf('%f s\n', toc);
% 
dv = p(1);
dwz = p(2);
kappaF = p(3);
JP = p(4);
FBx = p(5);
FBy = p(6);
FBz = p(7);
tauBx = p(8);
tauBy = p(9);
tauBz = p(10);

% dwz = 0;
% JP = 4.2e-5;
kappaT = pP1*JP;

sx = 0;
sy = 0;
sz = 0;
Jxy = 0;
Jxz = 0;
Jyz = 0;
dvxy = dv;
dvz = dv;
dwxy = 0;

fprintf('center of mass:        [ sx, sy, sz ] = [ %1.2f, %1.2f, %1.2f ] mm \n', sx*1e3, sy*1e3, sz*1e3);
fprintf('deviation inertia:     [ Jxy, Jxz, Jyz ] = [ %1.2e, %1.2e, %1.2e ] kg*m^2 \n', Jxy, Jxz, Jyz);
fprintf('translational damping: dvxy = %1.2e, dvz = %1.2e kg/s \n', dvxy, dvz);
fprintf('rotational damping:    dwxy = %1.2e, dwz = %1.2e kg*m^2/s \n', dwxy, dwz);
fprintf('bias force:            FB = [ %1.2e, %1.2e, %1.2e ] N, tauB = [ %1.2e, %1.2e, %1.2e ] Nm \n', FBx, FBy, FBz, tauBx, tauBy, tauBz);
fprintf('propeller param:       JP = %1.3e , kappaF = %1.3e, kappaT = %1.3e \n', JP, kappaF, kappaT);


%%
fprintf('Computing model accelerations... '); tic;

nE = 6;
xidH2 = zeros(nE, N);
xidH3 = zeros(nE, N);
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
  
  wP1 = wPSG(1,k);
  wP2 = wPSG(2,k);
  wP3 = wPSG(3,k);
  wP4 = wPSG(4,k);
  wP1d = wPdSG(1,k);
  wP2d = wPdSG(2,k);
  wP3d = wPdSG(3,k);
  wP4d = wPdSG(4,k);
  
  Mb = [m 0 0 0 0 0; 0 m 0 0 0 0; 0 0 m 0 0 0; 0 0 0 Jxx 0 0; 0 0 0 0 Jyy 0; 0 0 0 0 0 Jzz;];
  fb = [-g * m * Rzx + m * vy * wz - m * vz * wy - dv * vx + FBx; -g * m * Rzy - m * vx * wz + m * vz * wx - dv * vy + FBy; -g * m * Rzz + kappaF * wP1 ^ 2 + kappaF * wP2 ^ 2 + kappaF * wP3 ^ 2 + kappaF * wP4 ^ 2 + m * vx * wy - m * vy * wx - dv * vz + FBz; aA * kappaF * wP2 ^ 2 - aA * kappaF * wP4 ^ 2 + Jyy * wy * wz - Jzz * wy * wz + tauBx; -aA * kappaF * wP1 ^ 2 + aA * kappaF * wP3 ^ 2 - Jxx * wx * wz + Jzz * wx * wz + tauBy; -wP1 ^ 2 * pP1 * JP + wP2 ^ 2 * pP1 * JP - wP3 ^ 2 * pP1 * JP + wP4 ^ 2 * pP1 * JP + Jxx * wx * wy - Jyy * wx * wy - JP * wP1d + JP * wP2d - JP * wP3d + JP * wP4d - dwz * wz + tauBz;];

  xidH2(:,k) = Mb \ fb;
  
  Mb = [m 0 0 0 0 0; 0 m 0 0 0 0; 0 0 m 0 0 0; 0 0 0 Jxx 0 0; 0 0 0 0 Jyy 0; 0 0 0 0 0 Jzz;];
  fb = [-g * m * Rzx + m * vy * wz - m * vz * wy - dv * vx + FBx; -g * m * Rzy - m * vx * wz + m * vz * wx - dv * vy + FBy; -g * m * Rzz + kappaF * wP1 ^ 2 + kappaF * wP2 ^ 2 + kappaF * wP3 ^ 2 + kappaF * wP4 ^ 2 + m * vx * wy - m * vy * wx - dv * vz + FBz; aA * kappaF * wP2 ^ 2 - aA * kappaF * wP4 ^ 2 + Jyy * wy * wz - Jzz * wy * wz + tauBx; -aA * kappaF * wP1 ^ 2 + aA * kappaF * wP3 ^ 2 - Jxx * wx * wz + Jzz * wx * wz + tauBy; -wP1 ^ 2 * pP1 * JP + wP2 ^ 2 * pP1 * JP - wP3 ^ 2 * pP1 * JP + wP4 ^ 2 * pP1 * JP + Jxx * wx * wy - Jyy * wx * wy - dwz * wz + tauBz;];
  
  xidH3(:,k) = Mb \ fb;
end
fprintf('%f\n', toc);

%%
ColorMeas = 0.7*[1 1 1];
LineWidthMeas = 3;
LineStyleMeas = '-';

LineWidthIdent = 1;

ColorIdent1 = 'b';
ColorIdent2 = 'c';
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
ylim(40*[-1 1]);
ylabel('wyd in RAD/s^2', 'Interpreter', 'none');

sp_xi(4) = subplot(5,1,4); grid on; hold on;
line(tSG, xidSG(6,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(6,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH2(6,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
%line(tSG, xidH3(6,:), 'Color', ColorIdent3, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent3);
ylim(3*[-1 1]);
ylabel('wzd in RAD/s2', 'Interpreter', 'none');

sp_xi(5) = subplot(5,1,5); grid on; hold on;
line(tSG, xidSG(7,:), 'Color', ColorMeas, 'Linewidth', LineWidthMeas, 'LineStyle', LineStyleMeas);
line(tSG, xidH1(7,:), 'Color', ColorIdent1, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent1);
line(tSG, xidH1(7,:), 'Color', ColorIdent2, 'Linewidth', LineWidthIdent, 'LineStyle', LineStyleIdent2);
ylim(2200*[-1 1]);
ylabel('wP1d in RAD/s^2', 'Interpreter', 'none');
xlabel('time in s', 'Interpreter', 'none');

linkaxes(sp_xi, 'x');

axes(sp_xi(4));
%l_xi = legend('measured', 'full model', 'simp. model', 'oversimp. model', 'location', 'northeast');
l_xi = legend('measured', 'full model', 'simp. model', 'location', 'northeast');


%set(groot, 'DefaultTextInterpreter', 'LaTeX');
%set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
%set(groot, 'DefaultAxesFontName', 'LaTeX');
%set(groot, 'DefaultLegendInterpreter', 'LaTeX');

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
  
  JokerPrintFig( fig, 'ModelValidationQuad', 'pdf', 0 );
  
end



