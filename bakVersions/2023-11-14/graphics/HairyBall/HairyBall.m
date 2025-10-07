clear;

[X, Y, Z] = sphere(10);

alpha = linspace(0, pi, 31);
phi = linspace(0, 2*pi, 3*4*3+1);

Alpha = alpha' * ones(size(phi));
Phi = ones(size(alpha))' * phi;

Rxz = cos(Phi).*sin(Alpha);
Ryz = sin(Phi).*sin(Alpha);
Rzz = cos(Alpha);

Ones = ones(size(Rzz));
Rxx = Ones - Rxz.^2./(Ones+Rzz);
Ryx = -Rxz.*Ryz./(Ones+Rzz);
Rzx = -Rxz;
Rxy = -Rxz.*Ryz./(Ones+Rzz);
Ryy = Ones - Ryz.^2./(Ones+Rzz);
Rzy = -Ryz;


fig = figure(1);
mesh(Rxz,Ryz,Rzz, 'EdgeColor', 0.75*[1 1 1], 'FaceColor', 0.95*[1 1 1]);
axis equal;
axis off;
view(135, 45);
%view(135, 40);
zoom(2);

d = 1.005;
c = 0.15;
Linewidth = 1.5;

for i = 1:3:size(Rzz,1)-1
  for j = 1:3:size(Rzz,2)
    line(d*[Rxz(i,j) Rxz(i,j)+c*Rxx(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryx(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzx(i,j)], 'Color', 'r', 'Linewidth', Linewidth);
    line(d*[Rxz(i,j) Rxz(i,j)+c*Rxy(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryy(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzy(i,j)], 'Color', 'g', 'Linewidth', Linewidth);
    line(d*[Rxz(i,j) Rxz(i,j)+c*Rxz(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryz(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzz(i,j)], 'Color', 'b', 'Linewidth', Linewidth);
  end
end

%% export
figSize = [8 8];
fig.Units = 'centimeters';
fig.Position = [0 0 figSize(1) figSize(2)];
%tightfig(fig);
%fig.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig, 'HairyBallUpper', 'pdf', 0 );

view(135, -45);
JokerPrintFig( fig, 'HairyBallLower', 'pdf', 0 );

break;
%%
figure(fig); clf;
mesh(Rxz,Ryz,Rzz, 'EdgeColor', 0.75*[1 1 1], 'FaceColor', 0.95*[1 1 1]);
axis equal;
axis off;
view(135, 45);
%view(135, 40);
zoom(1.6);

for i = 1:3:size(Rzz,1)-1
  for j = 1:3:size(Rzz,2)
    line(d*[Rxz(i,j) Rxz(i,j)+c*Rxx(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryx(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzx(i,j)], 'Color', 'r', 'Linewidth', Linewidth);
    line(d*[Rxz(i,j) Rxz(i,j)+c*Rxy(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryy(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzy(i,j)], 'Color', 'g', 'Linewidth', Linewidth);
    %line(d*[Rxz(i,j) Rxz(i,j)+c*Rxz(i,j)], d*[Ryz(i,j) Ryz(i,j)+c*Ryz(i,j)], d*[Rzz(i,j) Rzz(i,j)+c*Rzz(i,j)], 'Color', 'b', 'Linewidth', Linewidth);
  end
end

view(135, 90);
JokerPrintFig( fig, 'HairyBallTop', 'pdf', 0 );

view(-135, -90);
JokerPrintFig( fig, 'HairyBallBottom', 'pdf', 0 );

