function traj = BicopterCircleTraj( param, Ts, TSim, R, w )

t = 0:Ts:TSim;
N = length(t);

g = param.g;
eps = param.epsy;

%R = 2;
%w = 2*pi*0.5;
y = R*[cos(w*t); sin(w*t)];
yd = R*w*[-sin(w*t); cos(w*t)];
ydd = -w^2*y;
yddd = -w^2*yd;
ydddd = -w^2*ydd;

a = atan2(-ydd(1,:), ydd(2,:)+g);
sa = sin(a);
ca = cos(a);
w = (yddd(1,:).*ca + yddd(2,:).*sa)./(ydd(1,:).*sa - (ydd(2,:)+g).*ca);
wd = (ydddd(1,:).*ca + ydddd(2,:).*sa - 2*(yddd(1,:).*sa - yddd(2,:).*ca).*w + (ydd(1,:).*ca + (ydd(2,:)+g).*sa).*w.^2)./(ydd(1,:).*sa - (ydd(2,:)+g).*ca);

r = y + eps*[sa; (1-ca)];
rd = yd + eps*[ca.*w; sa.*w];
rdd = ydd + eps*[ca.*wd; sa.*wd] + eps*[-sa.*w.^2; ca.*w.^2];

x = [r; a];
xi = [ca.*rd(1,:)+sa.*rd(2,:); -sa.*rd(1,:)+ca.*rd(2,:); w];
xid = [ca.*rdd(1,:)+sa.*rdd(2,:)-sa.*w.*rd(1,:)+ca.*w.*rd(2,:); -sa.*rdd(1,:)+ca.*rdd(2,:)-ca.*w.*rd(1,:)-sa.*w.*rd(2,:); wd];

rx = zeros(1, N);
ry = r(1,:);
rz = r(2,:);
Rxx = ones(1, N);
Ryx = zeros(1, N);
Rzx = zeros(1, N);
Rxy = zeros(1, N);
Ryy = cos(a);
Rzy = sin(a);
Rxz = zeros(1, N);
Ryz = -sin(a);
Rzz = cos(a);

vx = zeros(1, N);
vy = ca.*rd(1,:)+sa.*rd(2,:);
vz = -sa.*rd(1,:)+ca.*rd(2,:);
wx = w;
wy = zeros(1, N);
wz = zeros(1, N);

vxd = zeros(1, N);
vyd = ca.*rdd(1,:)+sa.*rdd(2,:)-sa.*w.*rd(1,:)+ca.*w.*rd(2,:);
vzd = -sa.*rdd(1,:)+ca.*rdd(2,:)-ca.*w.*rd(1,:)-sa.*w.*rd(2,:);
wxd = wd;
wyd = zeros(1, N);
wzd = zeros(1, N);

traj = struct;
traj.N = N;
traj.t = t;
traj.x = [ rx; ry; rz; Rxx; Ryx; Rzx; Rxy; Ryy; Rzy; Rxz; Ryz; Rzz ];
traj.xi = [ vx; vy; vz; wx; wy; wz ];
traj.xid = [ vxd; vyd; vzd; wxd; wyd; wzd ];

for k=1:traj.N
  traj.RPY(:,k) = RotMat2RollPitchYaw([traj.x(4:6,k) traj.x(7:9,k) traj.x(10:12,k)]);
end

end
