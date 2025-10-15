clear;

Linestyle1 = '-';
Linestyle2 = '--';

saveFigures = 1;

%% param
% model parameters of QuadV4
param.m = 1;
param.J1 = 0.015;
param.J2 = param.J1;
param.J3 = 0.03;
param.aG = 9.81;
param.lF = 0.24;
% desired closed loop roots
param.w01 = 5; % X,Y,Z translation
param.zeta1 = 0.8;
param.w02 = 20; % roll pitch
param.zeta2 = 0.95;
param.w03 = 3; % yaw
param.zeta3 = 0.8;

% compute other constant control parameters
param = QuadCtrl_Param(param);

A = [1.5; 1.5; 0];
Ts = 2e-2;
T_0 = 0.5;
T_Transit = 2.1;
T_Circle = 3.8;

TSnap = [0.1 1.6 2:0.25:3.5];

%% reference traj
tic;
%traj = QuadTraj_Transition( param, [0;0;0], [1;0;0], [0 0.2 1.2 3], Ts );
traj = QuadTraj_Circle( param, A, T_0, T_Transit, T_Circle, Ts, 1);
TSim = traj.t(end);

traj.r(3,:) = traj.r(3,:) + 1;
traj.x(3,:) = traj.r(3,:);

alpha = acos(traj.x(12,:));
phi = atan2(traj.x(5,:)-traj.x(7,:), traj.x(4,:)+traj.x(8,:));

toc;

%% setup plots
fig1 = figure(1); clf;
sp_r   = subplot(5,1,1); grid on;
l(1) = line(traj.t, traj.x(1,:), 'Color', 'r', 'Linewidth', 1);
l(2) = line(traj.t, traj.x(2,:), 'Color', 'g', 'Linewidth', 1);
l(3) = line(traj.t, traj.x(3,:), 'Color', 'b', 'Linewidth', 1);
ylabel('m', 'Interpreter', 'LaTeX');
legend(l, {'$r_{\rm{x}}$', '$r_{\rm{y}}$', '$r_{\rm{z}}$'}, 'Interpreter', 'LaTeX');
ylim(max(A)*[-1.1 1.1]);
sp_p   = subplot(5,1,2); grid on;
l(1) = line(traj.t, traj.RPY(1,:)*180/pi, 'Color', 'r', 'Linewidth', 1);
l(2) = line(traj.t, traj.RPY(2,:)*180/pi, 'Color', 'g', 'Linewidth', 1);
l(3) = line(traj.t, traj.RPY(3,:)*180/pi, 'Color', 'b', 'Linewidth', 1);
set(sp_p, 'ytick', 90*[-2 -1 0 1 2]);
ylim(180*[-1 1]); 
ylabel('deg', 'Interpreter', 'LaTeX');
legend(l, {'roll', 'pitch', 'yaw'}, 'Interpreter', 'LaTeX');
sp_a   = subplot(5,1,3); grid on;
l(1) = line(traj.t, alpha*180/pi, 'Color', 'r', 'Linewidth', 1);
l(2) = line(traj.t, phi*180/pi, 'Color', 'b', 'Linewidth', 1);
set(sp_a, 'ytick', 90*[-2 -1 0 1 2]);
ylim(180*[-1 1]); 
ylabel('deg', 'Interpreter', 'LaTeX');
legend(l(1:2), {'tilt', 'heading'}, 'Interpreter', 'LaTeX');
sp_F = subplot(5,1,4); grid on;
l(1) = line(traj.t, traj.u(1,:), 'Color', 'k', 'Linewidth', 1);
sp_F.YLim = [8 21];
legend(l(1), {'$F_{\rm{z}}$'}, 'Interpreter', 'LaTeX');
ylabel('N', 'Interpreter', 'LaTeX');
sp_tau = subplot(5,1,5); grid on;
l(1) = line(traj.t, traj.u(2,:), 'Color', 'r', 'Linewidth', 1);
l(2) = line(traj.t, traj.u(3,:), 'Color', 'g', 'Linewidth', 1);
l(3) = line(traj.t, traj.u(4,:), 'Color', 'b', 'Linewidth', 1);
%line(traj.t, sqrt([1 1 1] * traj.u(2:4,:).^2), 'Color', 'k', 'Linewidth', 1);
ylim(0.35*[-1 1]);
ylabel('Nm', 'Interpreter', 'LaTeX');
legend(l, {'$\tau_{\rm{x}}$', '$\tau_{\rm{y}}$', '$\tau_{\rm{z}}$'}, 'Interpreter', 'LaTeX');

xlabel('time in s');
linkaxes([sp_r sp_p sp_a sp_F sp_tau], 'x');
xlim([0 TSim]);

drawnow;

%%
fig2 = figure(2); clf;
ax2 = axes;
QuadSnapshots(ax2, param, traj.t, traj.x, 0*traj.x, zeros(4,length(traj.t)), TSnap);
view(ax2, 120, 35);

ax2.XLim = 1.8*[-1 1];
ax2.YLim = 1.8*[-1 1];
ax2.ZLim = [0 1.4];
ax2.XTick = [-2:2];
ax2.YTick = [-2:2];
ax2.ZTick = [0 1];

%% 
tic;
traj = QuadTraj_Circle( param, A, T_0, T_Transit, T_Circle, Ts, 0);
TSim = traj.t(end);

alpha = acos(traj.x(12,:));
phi = atan2(traj.x(5,:)-traj.x(7,:), traj.x(4,:)+traj.x(8,:));

toc;

axes(sp_r);
% line(traj.t, traj.x(1,:), 'Color', 'r', 'Linewidth', 1, 'LineStyle', Linestyle2);
% line(traj.t, traj.x(2,:), 'Color', 'g', 'Linewidth', 1, 'LineStyle', Linestyle2);
% line(traj.t, traj.x(3,:), 'Color', 'b', 'Linewidth', 1, 'LineStyle', Linestyle2);
axes(sp_p);
line(traj.t, traj.RPY(1,:)*180/pi, 'Color', 'r', 'Linewidth', 1, 'LineStyle', Linestyle2);
line(traj.t, traj.RPY(2,:)*180/pi, 'Color', 'g', 'Linewidth', 1, 'LineStyle', Linestyle2);
line(traj.t, traj.RPY(3,:)*180/pi, 'Color', 'b', 'Linewidth', 1, 'LineStyle', Linestyle2);
axes(sp_a);
%line(traj.t, alpha*180/pi, 'Color', 'r', 'Linewidth', 1, 'LineStyle', Linestyle2);
line(traj.t, phi*180/pi, 'Color', 'b', 'Linewidth', 1, 'LineStyle', Linestyle2);
axes(sp_F);
%line(traj.t, traj.u(1,:), 'Color', 'k', 'Linewidth', 1, 'LineStyle', Linestyle2);
axes(sp_tau);
line(traj.t, traj.u(2,:), 'Color', 'r', 'Linewidth', 1, 'LineStyle', Linestyle2);
line(traj.t, traj.u(3,:), 'Color', 'g', 'Linewidth', 1, 'LineStyle', Linestyle2);
line(traj.t, traj.u(4,:), 'Color', 'b', 'Linewidth', 1, 'LineStyle', Linestyle2);
%line(traj.t, sqrt([1 1 1] * traj.u(2:4,:).^2), 'Color', 'k', 'Linewidth', 1, 'LineStyle', Linestyle2);

drawnow;

%%
% fig3 = figure(3); clf;
% ax3 = axes;
% QuadSnapshots(ax3, param, traj.t, traj.x, 0*traj.x, zeros(4,length(traj.t)), TSnap);
% view(ax3, 40, 20);

%%
if (saveFigures > 0)
  figSize = [16 14];
  fig1.Units = 'centimeters';
  fig1.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig1);
  fig1.Position = [0 0 figSize(1) figSize(2)];

  JokerPrintFig( fig1, 'QuadCircleTrajGraph', 'pdf', 0 );
  
  figure(fig2);
  figSize = [16 9];
  fig2.Units = 'centimeters';
  fig2.Position = [0 0 figSize(1) figSize(2)];
  zoom(1.15); drawnow;

  JokerPrintFig( fig2, 'QuadCircleTrajSnapshots', 'pdf', 0 );
end

