clear;

TSim = 10;
saveFigure = 1;

CtrlModes = [1 2 3];
LineStyles = {'-', '--', ':'};

%% param
% model param
param.m = 1;
param.sx = 0;
param.sy = 0;
param.Jzz = 0.30;

% closed loop param
param.w0_v = 2;
param.zeta_v = 1;
param.w0_w = 6;
param.zeta_w = 0.95;

param.mc = param.m;
param.dc = param.mc * 2*param.zeta_v*param.w0_v;
param.kc = param.mc * param.w0_v^2;
param.Jczz = param.Jzz;
param.sigczz = param.Jczz * 2*param.zeta_w*param.w0_w;
param.kapczz = param.Jczz * param.w0_w^2;
param.scx = 0;
param.scy = 0;
param.lcx = 0;
param.lcy = 0;
param.hcx = 0;
param.hcy = 0;

%% initials
ini.r = [0; 1];
ini.a = 100 * pi/180;
ini.v = [0; 0];
ini.w = 0;
ini.z = [ini.r; sin(ini.a); cos(ini.a); ini.v; ini.w];

%% reference traj
traj.Ts = 0.005;
traj.t = 0:traj.Ts:TSim;

traj.r = [0; 0] * ones(size(traj.t));

% traj.A = 20*pi/180;
% traj.w0 = 2*pi*1;
% traj.a = traj.A * sin(traj.w0 * traj.t);
% traj.w = traj.A * traj.w0 * cos(traj.w0 * traj.t);
% traj.wd = -traj.w0^2 * traj.a;
traj.w0 = 2*pi*0.5;
traj.a = traj.w0 * traj.t;
traj.w = traj.w0 * ones(size(traj.t));
traj.wd = 0 * ones(size(traj.t));

traj.x = [ traj.r; sin(traj.a); cos(traj.a) ];
traj.xi = [ 0*traj.r; traj.w ];
traj.xid = [ 0*traj.r; traj.wd ];

traj.a = atan2(traj.x(3,:), traj.x(4,:));

%%
Mc = diag([param.mc param.mc]);
Dc = [param.dc -param.mc*traj.w0; param.mc*traj.w0 param.dc ];
Kc = diag([param.kc param.kc]);

e = polyeig(Kc, Dc, Mc)

e(1) = -(param.dc + i*param.mc*traj.w0)/2/param.mc + sqrt((param.dc + i*param.mc*traj.w0)^2/4/param.mc^2 - param.kc);
e(2) = -(param.dc - i*param.mc*traj.w0)/2/param.mc + sqrt((param.dc - i*param.mc*traj.w0)^2/4/param.mc^2 - param.kc);
e(3) = -(param.dc + i*param.mc*traj.w0)/2/param.mc - sqrt((param.dc + i*param.mc*traj.w0)^2/4/param.mc^2 - param.kc);
e(4) = -(param.dc - i*param.mc*traj.w0)/2/param.mc - sqrt((param.dc - i*param.mc*traj.w0)^2/4/param.mc^2 - param.kc);
e
break

%% setup plots
fig1 = figure(1); clf;
sp_r = subplot(5,1,1); grid on;
line_rR = line(traj.t, traj.x(1,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3);
line(traj.t, traj.x(2,:), 'Color', 0.8*[1 1 1], 'Linewidth', 3);
ylabel('position in m');
sp_a = subplot(5,1,2); grid on;
line_aR = line(traj.t, traj.a(1,:)*180/pi, 'Color', 0.8*[1 1 1], 'Linewidth', 3);
ylabel('orientation angle in DEG');
sp_F = subplot(5,1,3); grid on;
ylabel('Force in N');
sp_tau = subplot(5,1,4); grid on;
ylabel('torque in Nm');
sp_W = subplot(5,1,5); grid on;
ylabel('Total energy');
xlabel('t in s');
linkaxes([sp_r sp_a sp_F sp_tau sp_W], 'x');
xlim([0 TSim]);

fig2 = figure(2); clf;
sp_e = subplot(4,2,[1 3]); grid on; axis equal;
xlabel('$e_{\rm{x}}$ in m', 'interpreter', 'LaTeX');
ylabel('$e_{\rm{y}}$ in m', 'interpreter', 'LaTeX');
sp_e.XLim = 1.05*[-1 1];
sp_e.YLim = 1.05*[-1 1];
sp_e.XTick = [-1:0.5:1];
sp_e.YTick = [-1:0.5:1];
sp_rE = subplot(4,2,[2 4]); grid on; axis equal;
xlabel('$r_{\rm{xE}}$ in m', 'interpreter', 'LaTeX');
ylabel('$r_{\rm{yE}}$ in m', 'interpreter', 'LaTeX');
sp_rE.XLim = 1.05*[-1 1];
sp_rE.YLim = 1.05*[-1 1];
sp_rE.XTick = [-1:0.5:1];
sp_rE.YTick = [-1:0.5:1];
sp_Wr = subplot(4,2,[5 6]); grid on;
%xlabel('$t$ in s', 'interpreter', 'LaTeX');
ylabel('$\bar{W}_r$ in J', 'interpreter', 'LaTeX');
sp_Wr.YLim = [0 7];
sp_d = subplot(4,2,[7 8]); grid on;
xlabel('$t$ in s', 'interpreter', 'LaTeX');
ylabel('$||\mathbf{e}||$ in m', 'interpreter', 'LaTeX');
sp_d.YLim = 1.05*[0 1];

sp_Wr.XLim = [0 5];
sp_d.XLim = [0 5];

%% sim controlled system
for i = 1:length(CtrlModes)
  
  param.CtrlMode = CtrlModes(i);
  
  tic;
  options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3*ones(size(ini.z)), 'MaxStep', 1e-2);
  [t,z] = ode45(@ModelPlanarRigidBody, [0 TSim], ini.z, options, param, traj);
  
  res.t = t';
  res.x = interp1(t, z(:,1:4), res.t')';
  res.xi = interp1(t, z(:,5:7), res.t')';
  
  res.a = atan2(res.x(3,:), res.x(4,:));
  res.u = zeros(3, length(res.t));
  res.Wc = zeros(1, length(res.t));
  res.Pc = zeros(1, length(res.t));
  res.Wcr = zeros(1, length(res.t));
  for k=1:length(res.t)
    [res.u(:,k), res.Wc(k), res.Pc(k), res.Wcr(k)] = CtrlPlanarRigidBody(res.t(k), [res.x(:,k); res.xi(:,k)], param, traj);
  end
  
  res.xR = interp1(traj.t', traj.x', res.t')';
  res.e = res.x(1:2,:) - res.xR(1:2,:);
  res.rE = [res.xR(4,:).*res.e(1,:) + res.xR(3,:).*res.e(2,:); -res.xR(3,:).*res.e(1,:) + res.xR(4,:).*res.e(2,:)];
    
  %% plot
  axes(sp_r);
  line_rx = line(res.t, res.x(1,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  line_ry = line(res.t, res.x(2,:), 'Color', 'g', 'LineStyle', LineStyles{i});
  if (i == 1)
    legend([line_rR, line_rx, line_ry], {'ref', 'rx', 'ry'}, 'location', 'northeast');
  end
  axes(sp_a);
  line_a = line(res.t, res.a*180/pi, 'Color', 'r', 'LineStyle', LineStyles{i});
  if (i == 1)
    legend([line_aR, line_a], {'ref', '\alpha'}, 'location', 'northeast');
  end
  axes(sp_F);
  line_Fx = line(res.t, res.u(1,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  line_Fy = line(res.t, res.u(2,:), 'Color', 'g', 'LineStyle', LineStyles{i});
  if (i == 1)
    legend([line_Fx, line_Fy], {'F_x', 'F_y'}, 'location', 'northeast');
  end
  axes(sp_tau);
  line_tau = line(res.t, res.u(3,:), 'Color', 'r', 'LineStyle', LineStyles{i});
  if (i == 1)
    legend([line_tau], {'\tau_z'}, 'location', 'northeast');
  end
  axes(sp_W);
  line_Wc = line(res.t, res.Wc, 'Color', 'b', 'LineStyle', LineStyles{i});
  line_Pc = line(res.t, -2*res.Pc, 'Color', 'r', 'LineStyle', LineStyles{i});
  line(res.t(1:end-1)+diff(res.t)/2, diff(res.Wc)./diff(res.t), 'Color', 'm', 'Linestyle', '--', 'Linewidth', 2);
  if (i == 1)
    legend([line_Wc, line_Pc], {'W', '-2P'}, 'location', 'northeast');
  end
  
  axes(sp_e);
  line(res.e(1,:), res.e(2,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'LineWidth', 1);
  axes(sp_rE);
  line(res.rE(1,:), res.rE(2,:), 'Color', 'b', 'LineStyle', LineStyles{i}, 'LineWidth', 1);
  axes(sp_Wr);
  line(res.t, res.Wcr, 'Color', 'b', 'LineStyle', LineStyles{i}, 'LineWidth', 1);
  axes(sp_d);
  line(res.t, sqrt(res.e(1,:).^2 + res.e(2,:).^2), 'Color', 'b', 'LineStyle', LineStyles{i}, 'LineWidth', 1);
  
  toc
  
  drawnow;
end

%% save figure
if (saveFigure > 0)
  figSize = [16 16];
  fig2.Units = 'centimeters';
  fig2.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig2);
  fig2.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig2, 'RBCtrlPosError', 'pdf', 0 );
end
