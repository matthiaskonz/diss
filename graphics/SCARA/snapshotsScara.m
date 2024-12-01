function snapshotsScara(ax, param, t, x, xR, tSnap)

axes(ax);

% workspace limits
tmp = linspace(0,2*pi, 360);
line((param.l1-param.l2)*cos(tmp),(param.l1-param.l2)*sin(tmp), 'color', 'k');
line((param.l1+param.l2)*cos(tmp),(param.l1+param.l2)*sin(tmp), 'color', 'k');

% reference traj
rR = param.l1 * [cos(xR(1,:)); sin(xR(1,:))] + param.l2 * [cos(xR(1,:)+xR(2,:)); sin(xR(1,:)+xR(2,:))];
line(rR(1,:), rR(2,:), 'Color', 0.8*[1 1 1], 'linewidth', 3);

% actual traj
r = param.l1 * [cos(x(1,:)); sin(x(1,:))] + param.l2 * [cos(x(1,:)+x(2,:)); sin(x(1,:)+x(2,:))];
line(r(1,:), r(2,:), 'Color', 'b', 'linewidth', 1);

% SCARA snapshots
for i = 1:length(tSnap)
  k = find(t >= tSnap(i), 1, 'first');
  j1 = [0;0];
  j2 = j1 + param.l1 * [cos(x(1,k)); sin(x(1,k))];
  j3 = j2 + param.l2 * [cos(x(1,k)+x(2,k)); sin(x(1,k)+x(2,k))];
  jj = [j1 j2 j3];
  
  col = [0 0 1] + [1 1 0]*(1 - i/length(tSnap));
  line(jj(1,:) , jj(2,:), 'Color', col, 'Linewidth', 4, 'Marker', '.', 'MarkerSize', 30);
end

end

