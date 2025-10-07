clear;

r = [1; 2; 1.5];
R = RollPitchYaw2RotMat(pi/180*[25; 0; 55]);

ex = [1; 0; 0];
ey = [0; 1; 0];
ez = [0; 0; 1];

px = r + R*ex;
py = r + R*ey;
pz = r + R*ez;

pb = 2*[
 1.0  0.4 -0.8 0.3;
 1.0 -0.6  0.1 0.3;
 0.0  0.0  0.0 1.2;
];
pb = [
  2.0 -0.3 -0.3 0.2;
 -0.1  2.0 -0.1 0.2;
  0.0  0.0  0.0 1.5;
];

p = r*ones(1,size(pb,2)) + R*pb;

fig1 = figure(1); clf;
line([0 ex(1)], [0 ex(2)], [0 ex(3)], 'Linewidth', 2, 'Color', 'k');
line([0 ey(1)], [0 ey(2)], [0 ey(3)], 'Linewidth', 2, 'Color', 'k');
line([0 ez(1)], [0 ez(2)], [0 ez(3)], 'Linewidth', 2, 'Color', 'k');
line([r(1) px(1)], [r(2) px(2)], [r(3) px(3)], 'Linewidth', 2, 'Color', 'b');
line([r(1) py(1)], [r(2) py(2)], [r(3) py(3)], 'Linewidth', 2, 'Color', 'b');
line([r(1) pz(1)], [r(2) pz(2)], [r(3) pz(3)], 'Linewidth', 2, 'Color', 'b');

line(p(1,:), p(2,:), p(3,:), 'Linestyle', 'none', 'Color', 'c', 'Marker', '.', 'Markersize', 30);

for i=1:size(p, 2)
  for j=i+1:size(p, 2)
    line([p(1,i) p(1,j)], [p(2,i) p(2,j)], [p(3,i) p(3,j)], 'Linewidth', 1, 'Color', 'c');
  end
end

grid on;
axis equal;
axis off;
view(120, 20);


