function QuadSnapshots(ax, param, t, x, xR, F, Tsnap)

a = param.lF;
rProp = a/2;
ForceScale = .05;
BasisScale = 0.5;

rLand = .13;
dLandx = 0.6*rLand;
dLandz = .1;

%% setup figure
axes(ax);
%set(gcf,'color','w');
%sp1 = subplot(121);
axis equal
grid on
%view(30, 15);
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;

% sp2 = subplot(122);
% axis equal
axis([...
 min(x(1,:))-2*a max(x(1,:))+2*a ...
 min(x(2,:))-2*a max(x(2,:))+2*a ...
 min(x(3,:))-2*a max(x(3,:))+2*a ...
 ]) ;
% grid on
% view(30, 20);
% xlabel('x') ;
% ylabel('y') ;
% zlabel('z') ;

%% plot object templates

k = 1;
r = x(1:3,k);
q = x(4:7,k);
R = Quaternion2RotMat(q);
g = [R r; 0 0 0 1];

rR = xR(1:3,k);
qR = xR(4:7,k);

line(xR(1,:), xR(2,:), xR(3,:), 'Color', 'k', 'Linewidth', 1);
line(x(1,:), x(2,:), x(3,:), 'Color', 'g', 'Linewidth', 1);

% Frame MockUp
center = [0;0;0;1];
arm1 = [a;0;0;1];  % rotor 1
arm2 = [0;a;0;1];  % rotor 2
arm3 = [-a;0;0;1]; % rotor 3
arm4 = [0;-a;0;1]; % rotor 4
land11 = [dLandx;0;0;1];
land21 = [0;dLandx;0;1];
land31 = [-dLandx;0;0;1];
land41 = [0;-dLandx;0;1];
land12 = [rLand;0;-dLandz;1];
land22 = [0;rLand;-dLandz;1];
land32 = [-rLand;0;-dLandz;1];
land42 = [0;-rLand;-dLandz;1];

tmp = linspace(0, 2*pi, 20);
propCircle = [rProp*cos(tmp); rProp*sin(tmp); zeros(size(tmp)); zeros(size(tmp))];
LandCircle = [rLand*cos(tmp); rLand*sin(tmp); -dLandz*ones(size(tmp)); zeros(size(tmp))];

breakLine = [NaN; NaN; NaN; 1];
MockUp0 = [
 center arm2 breakLine ...
 center arm3 breakLine ...
 center arm4 breakLine ...
 arm2*ones(size(tmp))+propCircle breakLine ...
 arm3*ones(size(tmp))+propCircle breakLine ...
 arm4*ones(size(tmp))+propCircle breakLine ...
 land11 land12 breakLine ...
 land21 land22 breakLine ...
 land31 land32 breakLine ...
 land41 land42 breakLine ...
 center*ones(size(tmp))+LandCircle ...
];
MockUpProp10 = [
 center arm1 breakLine ...
 arm1*ones(size(tmp))+propCircle
];

% MockUpObj = line(MockUp(1,:), MockUp(2,:), MockUp(3,:), 'Color', 'b', 'LineWidth', 3);
% 
% ForceLines0 = [
%  arm1 arm1+[0;0;ForceScale*F(1,k); 0] breakLine ...
%  arm2 arm2+[0;0;ForceScale*F(2,k); 0] breakLine ...
%  arm3 arm3+[0;0;ForceScale*F(3,k); 0] breakLine ...
%  arm4 arm4+[0;0;ForceScale*F(4,k); 0]
% ];
% ForceLines = g*ForceLines0;
% 
% timeTextObj = title(sprintf('t = %2.2fs', t(1)));
% ForceLinesObj = line(ForceLines(1,:), ForceLines(2,:), ForceLines(3,:), 'Color', 'r', 'LineWidth', 3);

line(BasisScale*[0 1 0 0 0 0], BasisScale*[0 0 0 1 0 0], BasisScale*[0 0 0 0 0 1], 'Color', 'k', 'LineWidth', 3)
  
for T = Tsnap
  k = find(t >= T, 1, 'first');

  r = x(1:3,k);
  R = [x(4:6,k) x(7:9,k) x(10:12,k)];
  g = [R r; 0 0 0 1];

  MockUp = g*MockUp0;
	line(MockUp(1,:), MockUp(2,:), MockUp(3,:), 'Color', 'b', 'LineWidth', 2);
  MockUpProp1 = g*MockUpProp10;
	line(MockUpProp1(1,:), MockUpProp1(2,:), MockUpProp1(3,:), 'Color', 'c', 'LineWidth', 2);

%   ForceLines0(3,2) = ForceScale*F(1,k);
%   ForceLines0(3,5) = ForceScale*F(2,k);
%   ForceLines0(3,8) = ForceScale*F(3,k);
%   ForceLines0(3,11) = ForceScale*F(4,k);
%   ForceLines = g*ForceLines0;
%   set(ForceLinesObj, 'XData', ForceLines(1,:), 'YData', ForceLines(2,:), 'ZData', ForceLines(3,:));
  
end

drawnow;

