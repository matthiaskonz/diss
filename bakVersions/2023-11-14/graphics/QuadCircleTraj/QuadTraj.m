function traj = QuadTraj( param, traj, optHeading )

N = length(traj.t);

gVec = [0;0;-param.aG];

e3 = zeros(3,N);
e3d = zeros(3,N);
e3dd = zeros(3,N);
e3ddd = zeros(3,N);
e3dddd = zeros(3,N);
for k=1:N
  [ ~, ~, ~, ~, ~, e3(:,k), e3d(:,k), e3dd(:,k), e3ddd(:,k), e3dddd(:,k) ] = Vec2MagDir( param.m*(traj.rdd(:,k)-gVec), param.m * traj.rddd(:,k), param.m * traj.rdddd(:,k), param.m * traj.rddddd(:,k), param.m * traj.rdddddd(:,k) );
end

if (optHeading == 0)
  Rvec  = zeros(9,N);
  w     = zeros(3,N);
  wd    = zeros(3,N);
  wdd   = zeros(3,N);
  wddd  = zeros(3,N);
  for k=1:N
    [R, w(:,k), wd(:,k), wdd(:,k), wddd(:,k)] = ez2MinRotMat(e3(:,k), e3d(:,k), e3dd(:,k), e3ddd(:,k), e3dddd(:,k));
    Rvec(:,k) = [R(:,1); R(:,2); R(:,3)];
  end
elseif (optHeading == 1)
  [R0, ~, ~, ~, ~] = ez2MinRotMat(e3(:,1), zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1));
  [Rvec, w, wd, wdd, wddd] = VecRotVel2RotMat( traj.t, R0, e3d, e3dd, e3ddd, e3dddd, zeros(1,N), zeros(1,N), zeros(1,N), zeros(1,N) );
else
  error('Specify option for heading');
end

traj.x = [traj.r; Rvec];
traj.xi = zeros(6, N);
traj.xid = zeros(6, N);

traj.u = zeros(6, N);
traj.RPY = zeros(3, N);

for k=1:length(traj.t)
  R = reshape(Rvec(:,k), [3,3]);
  traj.RPY(:,k) = RotMat2RollPitchYaw(R);
  traj.xi(:,k) = [R' * traj.rd(:,k); w(:,k)];
  traj.xid(:,k) = [R' * traj.rdd(:,k) + (R*SkwMat(w(:,k)))' * traj.rd(:,k); wd(:,k)];
  
  fV = [param.aG * param.m * R(3,:)'; 0; 0; 0;];
  traj.u(:,k) = param.iB*(param.M*traj.xid(:,k) - aad_se3(traj.xi(:,k))' * param.M * traj.xi(:,k) + fV); 
end

traj.wdd = wdd;
traj.wddd = wddd;

end

