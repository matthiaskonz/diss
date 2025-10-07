clear;

saveFigure = 1;
figSize = [15.9 10];

CtrlModes = [1 2 0];
LineStyles = {'-', '--', ':'};

%% param
param = struct;

param.Jz = 0.2;

param.w0 = 3;
param.zeta = 0.7;

param.Jcz = 1;
param.kapcz = param.w0^2;
param.sigcz = 2 * param.zeta * param.w0;


%% initials
a_0 = -45*pi/180;
w_0 = 0;
z_0 = [a_0; w_0];

%% reference traj
Ts = 0.05;
X = [
  -pi/2 -pi/2 pi/2 0 0;
  0 0 0 0 0;
  0 0 0 0 0;
  0 0 0 0 0;
  ];
T = [0 0.5 2.5 2.6 3];
[xR, tR] = TransitionTrajectory(X, T, Ts);

A = 45 * pi/180;
w0 = 2*pi * 0.4;
xR = [ A*sin(w0*tR); A*w0*cos(w0*tR); -A*w0^2*sin(w0*tR) ];

param.tR   = tR;
param.xR   = xR(1,:);
param.xiR  = xR(2,:);
param.xiRd = xR(3,:);
param.uR = param.Jz*param.xiRd;

TSim = tR(end);

%% setup plots
fig = figure(1); clf;
sp_x = subplot(3,1,1); grid on;
line(param.tR, 180/pi*param.xR, 'Color', 0.8*[1 1 1], 'Linewidth', 3);
ylabel('a in DEG');
sp_u = subplot(3,1,2); grid on;
line(param.tR, param.uR, 'Color', 0.8*[1 1 1], 'Linewidth', 3); 
ylabel('tau in Nm');
sp_W = subplot(3,1,3); grid on;
ylabel('Wc in J');
xlabel('t in s');
linkaxes([sp_x sp_u sp_W], 'x');
xlim([0 TSim]);

drawnow;


%% sim controlled system
for i = 1:length(CtrlModes)
  param.CtrlMode = CtrlModes(i);
  
  tic;
  options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3*ones(size(z_0)), 'MaxStep', 5e-2);
  [t,z] = ode45(@ModelRevoluteJoint, [0 TSim], z_0, options, param);
  
  t = t';
  z = z';
  a = z(1,:);
  
  u  = zeros(1, length(t));
  Wc = zeros(1, length(t));
  Pc = zeros(1, length(t));
  for k=1:length(t)
    [u(:,k), Wc(k), Pc(k)] = CtrlRevoluteJoint(t(k), z(:,k), param);
  end
  toc
  
  % plot
  axes(sp_x);
  line(t, 180/pi*(mod(a+pi, 2*pi)-pi), 'Color', 'r', 'LineStyle', LineStyles{i});
  axes(sp_u);
  line(t, u, 'Color', 'r', 'LineStyle', LineStyles{i});
  axes(sp_W);
  line(t, Wc, 'Color', 'b', 'LineStyle', LineStyles{i});
  %line(t, -2*Pc, 'Color', 'r', 'LineStyle', LineStyles{i});
  %legend('Hc', 'Hcd');
  %line(t(1:end-1)+diff(t)/2, diff(Wc)./diff(t), 'Color', 'm', 'LineStyle', '--');
  
  drawnow;
  
end

%%
if (saveFigure)
  sp_u.YLim = [-1.5 3.5];
  sp_W.YLim = [0 6];
  
  sp_x.UserData.LatexYLabel = '$\jointAngle$ in $\unit{DEG}$';
  sp_u.UserData.LatexYLabel = '$\tau$ in $\unit{Nm}$';
	sp_W.UserData.LatexYLabel = '$\totalEnergyC$ in $\unit{J}$';
  
  %l2.UserData.LatexString = {'$e_{\rm{Obs}}$', '$e_{\rm{Ctrl}}$'};
  
  fig.Units = 'centimeters';
  fig.Position = [0 0 figSize(1) figSize(2)];
  tightfig(fig);
  fig.Position = [0 0 figSize(1) figSize(2)];
  
  JokerPrintFig( fig, 'RevoluteJointSimRes', 'pdf', 0 );
end
