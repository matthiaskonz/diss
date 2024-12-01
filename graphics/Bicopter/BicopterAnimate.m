function BicopterAnimate(fig, param, t, x, xR, F, p, timeScale)

lY = param.lY;
lZ = param.lZ;
betaF = param.betaF;
rProp = .5*lY;

ForceScale = param.ForceScale;
BasisScale = param.BasisScale;

%% setup figure
figure(fig);
%set(gcf,'color','w');
axis equal
grid on
%view(30, 15);
xlabel('x') ;
ylabel('y') ;
zlabel('z') ;
axis([...
 min(x(1,:))-2*lY max(x(1,:))+2*lY ...
 min(x(2,:))-2*lY max(x(2,:))+2*lY ...
 min(x(3,:))-2*lY max(x(3,:))+2*lY ...
 ]);

%% frame models
Center = [0; 0; 0; 1];
Ex = [BasisScale; 0; 0; 1];
Ey = [0; BasisScale; 0; 1];
Ez = [0; 0; BasisScale; 1];
Prop1 = [0; lY; lZ; 1];
Prop2 = [0; -lY; lZ; 1];
Top = [0; 0; lZ+lY*tan(betaF); 1];
Body1 = [0; .4*lY; lZ+.6*lY*tan(betaF); 1];
Body2 = [0; .5*lY; 0; 1];
Body3 = [0; -.5*lY; 0; 1];
Body4 = [0; -.4*lY; lZ+.6*lY*tan(betaF); 1];
BreakLine = [NaN; NaN; NaN; 1];

RRFrame0 = [Ex Center Ey BreakLine Center Ez];

MockUp0 = [ Prop1 Top Prop2 BreakLine Body1 Body2 Body3 Body4 ];

tmp = linspace(0, 2*pi, 20);
PropCircle0 = [rProp*cos(tmp); rProp*sin(tmp); zeros(size(tmp)); ones(size(tmp))];

%% plot object templates
k = 1;
% ref traj
rR = xR(1:3,k);
RR = [ xR(4:6,k) xR(7:9,k) xR(10:12,k) ];
gR = [RR rR; 0 0 0 1];

rRTrajLineObj = line(rR(1,1:k), rR(2,1:k), rR(3,1:k), 'Color', 0.7*[1 1 1], 'Linewidth', 1.5);
RRFrame = gR*RRFrame0;
RRLineObj = line(RRFrame(1,:), RRFrame(2,:), RRFrame(3,:), 'Color', 0.7*[1 1 1], 'LineWidth', 3);

% quad pos line and mockup
r = x(1:3,k);
R = [ x(4:6,k) x(7:9,k) x(10:12,k) ];
p1 = p(1,k);
p2 = p(2,k);
F1 = F(1,k);
F2 = F(2,k);

G = [R r; 0 0 0 1];
GP1 = [cos(p1) 0 sin(p1) 0; -sin(betaF) * sin(p1) cos(betaF) sin(betaF) * cos(p1) lY; -cos(betaF) * sin(p1) -sin(betaF) cos(betaF) * cos(p1) lZ; 0 0 0 1;];
GP2 = [cos(p2) 0 sin(p2) 0; sin(betaF) * sin(p2) cos(betaF) -sin(betaF) * cos(p2) -lY; -cos(betaF) * sin(p2) sin(betaF) cos(betaF) * cos(p2) lZ; 0 0 0 1;];

rTrajLineObj = line(x(1,1:k), x(2,1:k), x(3,1:k), 'Color', 'g', 'Linewidth', 1);

MockUp = G*MockUp0;
MockUpObj = line(MockUp(1,:), MockUp(2,:), MockUp(3,:), 'Color', 'b', 'LineWidth', 3);

PropCircle1 = G*GP1*PropCircle0;
PropCircle1Obj = line(PropCircle1(1,:), PropCircle1(2,:), PropCircle1(3,:), 'Color', 'b', 'LineWidth', 3);
PropCircle2 = G*GP2*PropCircle0;
PropCircle2Obj = line(PropCircle2(1,:), PropCircle2(2,:), PropCircle2(3,:), 'Color', 'b', 'LineWidth', 3);

ForceLine1 = G*GP1*[0 0; 0 0; 0 ForceScale*F1; 1 1];
ForceLine1Obj = line(ForceLine1(1,:), ForceLine1(2,:), ForceLine1(3,:), 'Color', 'r', 'LineWidth', 3);
ForceLine2 = G*GP2*[0 0; 0 0; 0 ForceScale*F2; 1 1];
ForceLine2Obj = line(ForceLine2(1,:), ForceLine2(2,:), ForceLine2(3,:), 'Color', 'r', 'LineWidth', 3);

timeTextObj = title(sprintf('t = %2.2fs', t(1)));

tic;
while ~isempty(k)
  
  rR = xR(1:3,k);
  RR = [ xR(4:6,k) xR(7:9,k) xR(10:12,k) ];
  gR = [RR rR; 0 0 0 1];
  set(rRTrajLineObj, 'XData', xR(1,1:k), 'YData', xR(2,1:k), 'ZData', xR(3,1:k));
  RRFrame = gR*RRFrame0;
  set(RRLineObj, 'XData', RRFrame(1,:), 'YData', RRFrame(2,:), 'ZData', RRFrame(3,:));
  
  r = x(1:3,k);
  R = [ x(4:6,k) x(7:9,k) x(10:12,k) ];
  p1 = p(1,k);
  p2 = p(2,k);
  F1 = F(1,k);
  F2 = F(2,k);
  G = [R r; 0 0 0 1];
  GP1 = [cos(p1) 0 sin(p1) 0; -sin(betaF) * sin(p1) cos(betaF) sin(betaF) * cos(p1) lY; -cos(betaF) * sin(p1) -sin(betaF) cos(betaF) * cos(p1) lZ; 0 0 0 1;];
  GP2 = [cos(p2) 0 sin(p2) 0; sin(betaF) * sin(p2) cos(betaF) -sin(betaF) * cos(p2) -lY; -cos(betaF) * sin(p2) sin(betaF) cos(betaF) * cos(p2) lZ; 0 0 0 1;];

  set(rTrajLineObj, 'XData', x(1,1:k), 'YData', x(2,1:k), 'ZData', x(3,1:k));
  MockUp = G*MockUp0;
  set(MockUpObj, 'XData', MockUp(1,:), 'YData', MockUp(2,:), 'ZData', MockUp(3,:));
  PropCircle1 = G*GP1*PropCircle0;
  set(PropCircle1Obj, 'XData', PropCircle1(1,:), 'YData', PropCircle1(2,:), 'ZData', PropCircle1(3,:));
  PropCircle2 = G*GP2*PropCircle0;
  set(PropCircle2Obj, 'XData', PropCircle2(1,:), 'YData', PropCircle2(2,:), 'ZData', PropCircle2(3,:));

  ForceLine1 = G*GP1*[0 0; 0 0; 0 ForceScale*F1; 1 1];
  set(ForceLine1Obj, 'XData', ForceLine1(1,:), 'YData', ForceLine1(2,:), 'ZData', ForceLine1(3,:));
  ForceLine2 = G*GP2*[0 0; 0 0; 0 ForceScale*F2; 1 1];
  set(ForceLine2Obj, 'XData', ForceLine2(1,:), 'YData', ForceLine2(2,:), 'ZData', ForceLine2(3,:));
  
  % show time
  set(timeTextObj, 'string', sprintf('t = %2.2fs', t(k)));
  
  drawnow;
  
  k = find(t >= toc*timeScale, 1, 'first') ;
  
end

