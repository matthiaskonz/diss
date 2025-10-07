function param = BicopterCtrlParam(param, figNr)

%% model parameters 
m = param.m;
Jx = param.Jx;
Jy = param.Jy;
Jz = param.Jz;
g = param.g;
lY = param.lY;
lZ = param.lZ;
betaF = param.betaF;
bF = param.bF;

epsx = -Jy / m / lZ;
epsy = Jx / m * sin(betaF) * lY / (cos(betaF) * lY ^ 2 - sin(betaF) * lY * lZ);

param.epsx = epsx;
param.epsy = epsy;

%% XY control
w1   = param.w1;
chi1 = param.chi1;
w2   = param.w2;
chi2 = param.chi2;

polyCtrl = conv([1 2*chi1*w1 w1^2], [1 2*chi2*w2 w2^2]);
p3 = polyCtrl(2);
p2 = polyCtrl(3);
p1 = polyCtrl(4);
p0 = polyCtrl(5);

mc = m; % arbitrary choice
kc = -p0 * p1 * mc / (p0 * p3 - p1 * p2);
dc = -p1 ^ 2 * mc / (p0 * p3 - p1 * p2);
hcz = (p1 * epsy + g * p3) / p1;
scz = (epsy * p0 * p3 - epsy * p1 * p2 - g * p1) / (p0 * p3 - p1 * p2);
Jcy = mc * (hcz - scz) * (scz - epsx);
Jcx = mc * (hcz - scz) * (scz - epsy);
lcz = hcz;
sigcx = 0;
sigcy = 0;
kapcx = mc * g * (hcz - scz);
kapcy = mc * g * (hcz - scz);

%% YAW control
wZ   = param.wZ;
chiZ = param.chiZ;

Jcz = Jz;
sigcz = Jz*2*chiZ*wZ;
kapcz = Jz*wZ^2;

%% store
param.mc   = mc;
param.scz  = scz;
param.Jcx  = Jcx;
param.Jcy  = Jcy;
param.Jcz  = Jcz;

param.dc   = dc;
param.lcz  = lcz;
param.sigcx = sigcx;
param.sigcy = sigcy;
param.sigcz = sigcz;

param.kc   = kc;
param.hcz  = hcz;
param.kapcx = kapcx;
param.kapcy = kapcy;
param.kapcz = kapcz;

%% input constraints
FMin = param.FMin;
FMax = param.FMax;
TiltMin = param.FTiltMin;
TiltMax = param.FTiltMax;

% polygonal approximation of constraint
tilt = linspace(param.FTiltMax, param.FTiltMin, 4);
P = [ ...
 FMin*tan(TiltMin) FMin*tan(TiltMax) FMax*sin(tilt);
 FMin FMin FMax*cos(tilt)
];
P = [P P(:,1)];

% compute corresponding inequalities
C = size(P,2)-1;
Z = zeros(C, 2);
Zb = zeros(C, 1);
for k=1:C
  p1 = P(:, k);
  p2 = P(:, k+1);
  d = p2-p1;
  d = d/norm(d);
  m = [d(2) -d(1)];
  Z(k,:) = m;
  Zb(k,:) = m*p1;
end
param.Z = blkdiag(Z,Z);
param.Zb = [Zb; Zb];

% visualize
if (figNr > 0)
  fig = figure(figNr); clf;
  
  % original constraints
  tilt1 = linspace(param.FTiltMin, param.FTiltMax, 10);
  tilt2 = linspace(param.FTiltMin, param.FTiltMax, 20);
  A = [ ...
    param.FMin*sin(tilt1) param.FMax*sin(fliplr(tilt2));
    param.FMin*cos(tilt1) param.FMax*cos(fliplr(tilt2));
  ];
  A = [A A(:,1)];
  fill(A(1,:), A(2,:), 0.7*[1 1 1], 'LineStyle', '-', 'EdgeColor', 0.75*[1 1 1]);
  
  % approximated constraints
  line(P(1,:), P(2,:), 'Color', 'b', 'LineWidth', 0.75);
  
  xlabel('$u_1 = F_1 \sin\theta_1$', 'Interpreter', 'LaTeX');
  ylabel('$u_2 = F_1 \cos\theta_1$', 'Interpreter', 'LaTeX');
  axis equal;
  grid on;
  xlim(round([FMax*sin(TiltMin) FMax*sin(TiltMax)]) + [-1 1]);
  ylim([0 round(FMax)+1]);
  drawnow;
  
  figSize = [7 7];
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  JokerPrintFig( fig, 'BicopterInputConstraintApprox', 'pdf', 0 );
end
