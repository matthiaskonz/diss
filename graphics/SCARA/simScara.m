clear;

TSim = 5;
tSnap = [0 1 2.2 2.8 4];
saveFigure = 1;

CtrlModes = [1 5];
LineStyles = {'-', '--', '-.'};

%% model param
param.l1 = 0.3;
param.l2 = 0.2;
param.m1 = 1;
param.m2 = 2;
param.s2x = param.l2/2;
param.s2y = 0;
param.J1z = param.m1*param.l1^2/12;
param.J2z = param.m2*param.l2^2/12;

%% ctrl param
param.w0c = 4;
param.zetac = 1;

%% initials
z_0 = [pi/180*-30; pi/180*90; 0; 0];

%% reference traj
Ts = 0.005;

T = [0 1 4 TSim];
X = [
  0.4 0.4 -0.4 -0.4 ;
  0 0 0 0;
  0 0 0 0;
  ];
[x, t] = TransitionTrajectory(X, T, Ts);
Y = [
  0.2 0.2 0.2 0.2;
  0 0 0 0;
  0 0 0 0;
  ];
[y, t] = TransitionTrajectory(Y, T, Ts);

rx = x(1,:);
ry = y(1,:);
rxd = x(2,:);
ryd = y(2,:);
rxdd = x(3,:);
rydd = y(3,:);
clear T X Y x y

% rx = param.l1*ones(size(t));
% rxd = zeros(size(t));
% rxdd = zeros(size(t));
% A = 0.35;
% w = 2*pi*0.5;
% ry = A*cos(w*t);
% ryd = -A*w*sin(w*t);
% rydd = -A*w^2*cos(w*t);
% clear A w

a2 = -acos( (rx.^2+ry.^2-param.l1^2-param.l2^2) / (2*param.l1*param.l2) );
a1 = atan2( ry.*(param.l1+param.l2*cos(a2)) - rx * param.l2.*sin(a2), rx.*(param.l1+param.l2*cos(a2)) + ry * param.l2.*sin(a2) );
a1d = (rxd.*cos(a1+a2) + ryd.*sin(a1+a2)) ./ (param.l1*sin(a2));
a2d = -(rxd.*(param.l1*cos(a1) + param.l2*cos(a1+a2)) + ryd.*(param.l1*sin(a1) + param.l2*sin(a1+a2))) ./ (param.l1*param.l2*sin(a2));
a1dd = ((rxdd.*cos(a1+a2) + rydd.*sin(a1+a2)) + (-rxd.*sin(a1+a2)+ryd.*cos(a1+a2)).*(a1d+a2d) - param.l1*cos(a2).*a1d.*a2d) ./ (param.l1*sin(a2));
a2dd = -((rydd.*sin(a1) + rxdd.*cos(a1))*param.l1 + (rydd.*sin(a1+a2) + rxdd.*cos(a1+a2))*param.l2 + (ryd.*cos(a1)-rxd.*sin(a1))*param.l1.*a1d + (ryd.*cos(a1+a2) - rxd.*sin(a1+a2))*param.l2.*(a1d+a2d) + param.l1*param.l2*cos(a2).*a2d.^2 ) ./ (param.l1*param.l2*sin(a2));

traj.t = t;
traj.x = [a1; a2];
traj.xi = [a1d; a2d];
traj.xid = [a1dd; a2dd];
traj.r = param.l1 * [cos(a1); sin(a1)] + param.l2 * [cos(a1+a2); sin(a1+a2)];
clear t rx ry rxd ryd rxdd rydd a1 a2 a1d a2d a1dd a2dd

% figure;
% % line(tR, a2d);
% % line(tR(1:end-1)+.5*diff(a2), diff(a2)./diff(tR), 'Linestyle', '--', 'color', 'r');
% line(tR, a2dd);
% line(tR(1:end-1)+.5*diff(tR), diff(a2d)./diff(tR), 'Linestyle', '--', 'color', 'r');

%% setup plots
fig1 = figure(1); clf;
sp_x = subplot(4,1,1); grid on;
ylabel('joint angle in DEG', 'Interpreter', 'LaTex');
l_aR = line(traj.t, 180/pi*traj.x(1,:), 'Color', 0.8*[1 1 1], 'linewidth', 3);
line(traj.t, 180/pi*traj.x(2,:), 'Color', 0.8*[1 1 1], 'linewidth', 3);
sp_r = subplot(4,1,2); grid on;
ylabel('tool pos.\ in m', 'Interpreter', 'LaTex');
l_rR = line(traj.t, traj.r(1,:), 'Color', 0.8*[1 1 1], 'linewidth', 3);
line(traj.t, traj.r(2,:), 'Color', 0.8*[1 1 1], 'linewidth', 3);
sp_u = subplot(4,1,3); grid on;
ylabel('joint torques in Nm', 'Interpreter', 'LaTex');
sp_W = subplot(4,1,4); grid on;
ylabel('Total energy', 'Interpreter', 'LaTex');
xlabel('t in s', 'Interpreter', 'LaTex');
linkaxes([sp_x sp_r sp_u sp_W], 'x');
xlim([0 TSim]);

%% sim controlled system
odeOpt = odeset('RelTol', 1e-4, 'AbsTol', 1e-4*ones(size(z_0)), 'MaxStep', 1e-2);
for i = 1:length(CtrlModes)
  tic;
  param.CtrlMode = CtrlModes(i);
  [t,z] = ode45(@ModelScara, [0 TSim], z_0, odeOpt, param, traj);
  
  res.t = [t(1):0.05:t(end)];
  res.x = interp1(t, z(:,1:2), res.t')';
  res.xi = interp1(t, z(:,3:4), res.t')';
  clear t z;
  
  res.xR  = interp1(traj.t', traj.x', res.t')';
  res.xiR = interp1(traj.t', traj.xi', res.t')';
  res.xiRd = interp1(traj.t', traj.xid', res.t')';
  
  res.u = zeros(2, length(res.t));
  res.Wc = zeros(1, length(res.t));
  res.Pc = zeros(1, length(res.t));
  for k=1:length(res.t)
    [res.u(:,k), res.Wc(k), res.Pc(k)] = CtrlScara(res.t(k), [res.x(:,k); res.xi(:,k)], param, traj);
  end
  res.r = param.l1 * [cos(res.x(1,:)); sin(res.x(1,:))] + param.l2 * [cos(res.x(1,:)+res.x(2,:)); sin(res.x(1,:)+res.x(2,:))];  
  toc
  
  %% plot
  axes(sp_x);
  l_a1 = line(res.t, 180/pi*res.x(1,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  l_a2 = line(res.t, 180/pi*res.x(2,:), 'Color', 'g', 'LineStyle', LineStyles{i});
  if i==1; leg_a = legend([l_aR, l_a1, l_a2], {'ref', '$\theta_1$', '$\theta_2$'}, 'Interpreter', 'LaTeX'); end
  axes(sp_r);
  l_rx = line(res.t, res.r(1,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  l_ry = line(res.t, res.r(2,:), 'Color', 'g', 'LineStyle', LineStyles{i});
  if i==1; leg_r = legend([l_rR, l_rx, l_ry], {'ref', '$r_\mathrm{Tx}$', '$r_\mathrm{Ty}$'}, 'Interpreter', 'LaTeX'); end
  axes(sp_u);
  l_u1 = line(res.t, res.u(1,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  l_u2 = line(res.t, res.u(2,:), 'Color', 'g', 'LineStyle', LineStyles{i});
  if i==1; leg_u = legend([l_u1, l_u2], {'$\tau_1$', '$\tau_2$'}, 'Interpreter', 'LaTeX'); end
 % ylim(10*[-1 1]);
  axes(sp_W);
  l_W = line(res.t, res.Wc, 'Color', 'r', 'LineStyle', LineStyles{i});
  l_P = line(res.t, -2*res.Pc, 'Color', 'b', 'LineStyle', LineStyles{i});
  if i==1; leg_W = legend([l_W, l_P], {'$\bar{W}$', '$-2\bar{P}$'}, 'Interpreter', 'LaTeX'); end
  %line(res.t(1:end-1)+diff(res.t)/2, diff(res.Wc)./diff(res.t), 'Color', 'm', 'LineStyle', LineStyles{i});

  figSnap = figure(10+i); clf;
  axSnap{i} = subplot(1,1,1);
  axis equal; grid on;
  xlabel('$r_\mathrm{Tx}$', 'Interpreter', 'LaTeX');
  ylabel('$r_\mathrm{Ty}$', 'Interpreter', 'LaTeX');
  snapshotsScara(axSnap{i}, param, res.t, res.x, res.xR, tSnap)
  
  drawnow;
end

%%
if (saveFigure)
  fig = fig1;
  figSize = [16 13];
  fig.Units = 'centimeters';
  fig.Position = [fig.Position(1) fig.Position(2) figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [fig.Position(1) fig.Position(2) figSize(1) figSize(2)];
  
  sp_x.YLim = 180*[-1 1];
  sp_x.YTick = 90*[-2:2];
  sp_u.YLim = [-1 3];
  sp_W.YLim = [-95 50];
  leg_a.Location  = 'SouthEast';
  leg_r.Location  = 'NorthEast';
  leg_u.Location  = 'NorthEast';
  leg_W.Location  = 'SouthEast';
  
  JokerPrintFig(fig, 'ScaraSimRes', 'pdf', 0 );
  
  for i = 1:length(CtrlModes);
  fig = figure(10+i);
  figSize = [7.9 6];
  fig.Units = 'centimeters';
  fig.Position = [fig.Position(1) fig.Position(2) figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [fig.Position(1) fig.Position(2) figSize(1) figSize(2)];
  
  axSnap{i}.XLim = 0.51*[-1 1];
  axSnap{i}.YLim = [-0.2 0.51];
  axSnap{i}.XTick = -0.4:0.2:0.4;
  axSnap{i}.YTick = -0.4:0.2:0.4;
  JokerPrintFig(fig, sprintf('ScaraSnapApproach%i', i), 'pdf', 0 );
  end
  
end

