clear;

saveFigure = 1;

%%% Control modes:
%  1.. approach 1 (Bullo & Murray)
% 11.. approach 1 with alternative transport map
%  2.. approach 2 (Gauss)
%  3.. approach 3 

CtrlModes = [2];
LineStyles = {'-', '--', '-.'};

%% param
param = struct;
% model param
param.m = 1;
param.Jx = 0.017;
param.Jy = 0.010;
param.Jz = 0.025;
param.g = 9.81;
param.lY = 0.24;
param.lZ = 0.05;
param.betaF = 20*pi/180;
param.bF = 0.0;
param.FMin = 2;
param.FMax = 14;
param.FTiltMin = -30*pi/180;
param.FTiltMax = 30*pi/180;
% desired closed loop roots
param.w1   = 3;
param.chi1 = 0.8;
param.w2   = 12;
param.chi2 = 0.8;
param.wZ   = 9;
param.chiZ = 1;
% compute resulting control param
param = BicopterCtrlParam(param, 0);

%% initial conditions
r0  = [0.5; 1; -0.5];
v0 = [0; 0; 2];
R0  = RollPitchYaw2RotMat([0;0;-90]*pi/180);
w0 = [0;0;0];

res = struct;
res.x0 = [r0; R0(:,1); R0(:,2); R0(:,3)];
res.xi0 = [v0; w0];
clear r0 R0 v0 w0;

%% reference traj
Ts = 2e-2;

tic;
traj = BicopterTrajCircle( param, Ts, 5, 1.5, 2*pi*0.5 );
toc;

TSim = traj.t(traj.N);

%% setup plots
fig = figure(1); clf;
sp_r   = subplot(7,1,1:2); grid on;
posRefLine = line(traj.t, traj.x(1,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3);
line(traj.t, traj.x(2,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3);
line(traj.t, traj.x(3,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3);
ylabel('position in m', 'Interpreter', 'LaTeX');
sp_RPY = subplot(7,1,3:4); grid on;
RPYRefLine = line(traj.t, traj.RPY(1,:)*180/pi, 'Color', 0.8*[1 1 1], 'Linewidth', 3);
line(traj.t, traj.RPY(2,:)*180/pi, 'Color', 0.8*[1 1 1], 'Linewidth', 3);
line(traj.t, traj.RPY(3,:)*180/pi, 'Color', 0.8*[1 1 1], 'Linewidth', 3);
sp_RPY.YLim = 180*[-1 1];
sp_RPY.YTick = 180*[-1 -.5 0 .5 1];
ylabel('euler angles in DEG', 'Interpreter', 'LaTeX');
sp_F = subplot(7,1,5); grid on;
sp_F.YLim = [param.FMin param.FMax];
ylabel('thrust in N', 'Interpreter', 'LaTeX');
sp_tilt = subplot(7,1,6); grid on;
sp_tilt.YLim = 180/pi*[param.FTiltMin param.FTiltMax];
sp_tilt.YLim = 30*[-1 1];
sp_tilt.YTick = 15*[-6:1:6];
ylabel('tilt in DEG', 'Interpreter', 'LaTeX');
sp_W = subplot(7,1,7); grid on;
sp_W.YLim = [-20 10];
ylabel('Energies', 'Interpreter', 'LaTeX');
xlabel('time in s', 'Interpreter', 'LaTeX');
linkaxes([sp_r sp_RPY sp_F sp_tilt sp_W], 'x');
sp_r.XLim = [0 TSim];
drawnow;

%TSim = 3;

%% sim controlled system
for i = 1:length(CtrlModes)
  param.CtrlMode = CtrlModes(i);
  
  tic;
  z0 = [res.x0; res.xi0];
  simOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-3*ones(size(z0)), 'MaxStep', 1e-2);
  [t,z] = ode45(@BicopterModel, [0 TSim], z0, simOpt, param, traj);
  
  res.t = traj.t;
  res.N = length(res.t);
  res.x = interp1(t, z(:,1:12), res.t')';
  res.xi = interp1(t, z(:,13:18), res.t')';
  clear t z z0;

  res.u   = zeros(4, res.N);
  res.Wc  = zeros(1, res.N);
  res.Pc  = zeros(1, res.N);
  res.RPY = zeros(3, res.N);
  for k=1:res.N
    [res.u(:,k), res.Wc(k), res.Pc(k)] = BicopterCtrl(res.t(k), [res.x(:,k); res.xi(:,k)], param, traj);
    res.RPY(:,k) = RotMat2RollPitchYaw([res.x(4:6,k) res.x(7:9,k) res.x(10:12,k)]);
  end
  res.F = [sqrt(res.u(1,:).^2 + res.u(2,:).^2); sqrt(res.u(3,:).^2 + res.u(4,:).^2)];
  res.tilt = [atan2(res.u(1,:), res.u(2,:)); atan2(res.u(3,:), res.u(4,:))];
  toc
  
  % plot
  axes(sp_r);
  hLine(1) = line(res.t, res.x(1,:), 'Color', 'r', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(2) = line(res.t, res.x(2,:), 'Color', 'g', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(3) = line(res.t, res.x(3,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  if (i==1) leg_Pos = legend([posRefLine hLine(1:3)], {'ref' '$r_{\rm{x}}$', '$r_{\rm{y}}$', '$r_{\rm{z}}$'}, 'Interpreter', 'LaTeX'); end
  axes(sp_RPY);
  hLine(1) = line(res.t, 180/pi*res.RPY(1,:), 'Color', 'r', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(2) = line(res.t, 180/pi*res.RPY(2,:), 'Color', 'g', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(3) = line(res.t, 180/pi*res.RPY(3,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  if (i==1) leg_RPY = legend([RPYRefLine hLine(1:3)], {'ref' 'roll', 'pitch', 'yaw'}, 'Interpreter', 'LaTeX'); end
  axes(sp_F);
  hLine(1) = line(res.t, res.F(1,:), 'Color', 'r', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(2) = line(res.t, res.F(2,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  if (i==1) leg_F = legend(hLine(1:2), {'$F_1$', '$F_2$'}, 'Interpreter', 'LaTeX'); end
  axes(sp_tilt);
  hLine(1) = line(res.t, 180/pi*res.tilt(1,:), 'Color', 'r', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(2) = line(res.t, 180/pi*res.tilt(2,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  if (i==1) leg_Tilt = legend(hLine(1:2), {'$\theta_1$', '$\theta_2$'}, 'Interpreter', 'LaTeX'); end
  axes(sp_W);
  hLine(1) = line(res.t, res.Wc, 'Color', 'b', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(2) = line(res.t(1:end-1)+.5*diff(res.t), diff(res.Wc)./diff(res.t), 'Color', 'c', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  hLine(3) = line(res.t, -2*res.Pc, 'Color', 'r', 'LineStyle', LineStyles{i}, 'Linewidth', 1);
  if (i==1) leg_W = legend(hLine(1:3), {'$\bar{W}$', '$\dot{\bar{W}}$', '$\bar{P}$'}, 'Interpreter', 'LaTeX'); end
  
  drawnow;

end

if (saveFigure > 0)
  figSize = [16 23];
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  leg_Pos.Location  = 'NorthEast';
  leg_RPY.Location  = 'NorthEast';
  leg_F.Location    = 'NorthEast';
  leg_Tilt.Location = 'NorthEast';
  leg_W.Location    = 'NorthEast';
  JokerPrintFig( fig, 'BicopterCircleSimRes', 'pdf', 0 );
end



break

%% animation
fig = figure(2); clf;
view(135, 25);

param.ForceScale = 0.05;
param.BasisScale = 0.5;
%BicopterAnimate(fig, param, res.t, res.x, traj.x, res.F, res.tilt, 0.25);

if (saveFigure > 0)
  fig = figure(2); clf;
  view(135, 10);
  BicopterSnapshots(fig, param, res.t, res.x, traj.x, res.F, res.tilt, [0 .3 .65 1.1 1.6 2.1 2.5 3]);
  
  figSize = [14 14];
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  JokerPrintFig( fig, 'BicopterCircleSimSnapshots', 'pdf', 0 );
end