function param = QuadCtrl_Param(param)
% compute constant control parameters 

% mode param
m   = param.m;
J1  = param.J1;
J2  = param.J2;
J3  = param.J3;
aG  = param.aG;

param.M = [m 0 0 0 0 0; 0 m 0 0 0 0; 0 0 m 0 0 0; 0 0 0 J1 0 0; 0 0 0 0 J2 0; 0 0 0 0 0 J3;];
param.D = zeros(6,6);
param.B = [0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;];
param.iB = [0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0;];
param.Bt = [1 0 0 0 0 0; 0 1 0 0 0 0;];

% desired closed loop roots
w01   = param.w01; % X,Y,Z translation
zeta1 = param.zeta1;
w02   = param.w02; % roll pitch
zeta2 = param.zeta2;
w03   = param.w03; % yaw
zeta3 = param.zeta3;

mc = m; % arbitrary
% matched parameters [Maple code export "MatchingQuadcopter.mw"]
kc = w01 ^ 2 * w02 ^ 2 * (w01 * zeta2 + w02 * zeta1) / (4 * w01 ^ 2 * w02 * zeta1 * zeta2 ^ 2 + 4 * w01 * w02 ^ 2 * zeta1 ^ 2 * zeta2 + w01 ^ 3 * zeta2 + w02 ^ 3 * zeta1) * m;
dc = 2 * m * (w01 * zeta2 + w02 * zeta1) ^ 2 / (4 * w01 ^ 2 * w02 * zeta1 * zeta2 ^ 2 + 4 * w01 * w02 ^ 2 * zeta1 ^ 2 * zeta2 + w01 ^ 3 * zeta2 + w02 ^ 3 * zeta1) * w01 * w02;
sc = aG * (w01 * zeta2 + w02 * zeta1) / (4 * w01 ^ 2 * w02 * zeta1 * zeta2 ^ 2 + 4 * w01 * w02 ^ 2 * zeta1 ^ 2 * zeta2 + w01 ^ 3 * zeta2 + w02 ^ 3 * zeta1);
hc = aG * (w01 * zeta1 + w02 * zeta2) / w01 / w02 / (w01 * zeta2 + w02 * zeta1);

J1c = mc * sc * (hc - sc);
J2c = mc * sc * (hc - sc);
kap1 = mc * (hc - sc) * aG;
kap2 = mc * (hc - sc) * aG;
lc = hc;
sig1 = 0;
sig2 = 0;

kap3 = J3*w03^2;
sig3 = J3*2*zeta3*w03;

param.mc = mc;
param.kc = kc;
param.dc = dc;
param.sc = sc;
param.hc = hc;
param.kap3 = kap3;
param.sig3 = sig3;

param.Mc = [mc 0 0 0 sc * mc 0; 0 mc 0 -sc * mc 0 0; 0 0 mc 0 0 0; 0 -sc * mc 0 mc * sc ^ 2 + J1c 0 0; sc * mc 0 0 0 mc * sc ^ 2 + J2c 0; 0 0 0 0 0 J3;];
param.Dc = [dc 0 0 0 lc * dc 0; 0 dc 0 -lc * dc 0 0; 0 0 dc 0 0 0; 0 -lc * dc 0 dc * lc ^ 2 + sig1 0 0; lc * dc 0 0 0 dc * lc ^ 2 + sig2 0; 0 0 0 0 0 sig3;];
param.Kc = [kc 0 0 0 kc * hc 0; 0 kc 0 -kc * hc 0 0; 0 0 kc 0 0 0; 0 -kc * hc 0 kc * hc ^ 2 + kap1 0 0; kc * hc 0 0 0 kc * hc ^ 2 + kap2 0; 0 0 0 0 0 kap3;];

% for flatness based approach
polyCtrl = conv([1 2*zeta1*w01 w01^2], [1 2*zeta2*w02 w02^2]);
param.px3 = polyCtrl(2); 
param.px2 = polyCtrl(3);
param.px1 = polyCtrl(4);
param.px0 = polyCtrl(5);
param.pz0 = kc/mc;
param.pz1 = dc/mc;
param.pp0 = w03^2;
param.pp1 = 2*zeta3*w03;

% attitude gains for cascade approach
param.Kc0 = w01^2 * eye(3);
param.Kc1 = 2*zeta1*w01 * eye(3);
param.Kc_R = diag([w02^2; w02^2; w03^2]);
param.Dc_R = 2*diag([zeta2*w02; zeta2*w02; zeta3*w03]);

