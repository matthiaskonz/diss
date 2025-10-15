function traj = QuadTraj_Circle( param, A, T_0, T_Transit, T_Circle, Ts, optHeading)

Ax = A(1);
Ay = A(2);
Az = A(3);

w = 2*2*pi/T_Circle;

t0 = 0;
t1 = T_0;
t2 = t1 + T_Transit;
t3 = t2 + T_Circle;
t4 = t3 + T_Transit;
t5 = t4 + T_0;

traj.t = 0:Ts:t5;
N = length(traj.t);
k0 = 1;
k1 = find(traj.t >= t1, 1, 'first');
k2 = find(traj.t >= t2, 1, 'first');
k3 = find(traj.t >= t3, 1, 'first');
k4 = find(traj.t >= t4, 1, 'first');
k5 = N;

% the circle
traj.r      =      [ Ax*cos(w*(traj.t-t2)); -Ay*sin(w*(traj.t-t2)); -Az*sin(w*(traj.t-t2))];
traj.rd     =  w  *[-Ax*sin(w*(traj.t-t2)); -Ay*cos(w*(traj.t-t2)); -Az*cos(w*(traj.t-t2))];
traj.rdd    =  w^2*[-Ax*cos(w*(traj.t-t2));  Ay*sin(w*(traj.t-t2));  Az*sin(w*(traj.t-t2))];
traj.rddd     = -w^2*traj.rd;
traj.rdddd    = -w^2*traj.rdd;
traj.rddddd   = -w^2*traj.rddd;
traj.rdddddd  = -w^2*traj.rdddd;

% transition from origin to circle
T = [t0 t1 t2];
X = [
  0 0 traj.r(1,k2);
  0 0 traj.rd(1,k2);
  0 0 traj.rdd(1,k2);
  0 0 traj.rddd(1,k2);
  0 0 traj.rdddd(1,k2);
  0 0 traj.rddddd(1,k2);
  0 0 traj.rdddddd(1,k2)
  ];
[rx, ~] = TransitionTrajectory(X, T, Ts);
Y = [
  0 0 traj.r(2,k2);
  0 0 traj.rd(2,k2);
  0 0 traj.rdd(2,k2);
  0 0 traj.rddd(2,k2);
  0 0 traj.rdddd(2,k2);
  0 0 traj.rddddd(2,k2);
  0 0 traj.rdddddd(2,k2)
  ];
[ry, ~] = TransitionTrajectory(Y, T, Ts);
Z = [
  0 0 traj.r(3,k2);
  0 0 traj.rd(3,k2);
  0 0 traj.rdd(3,k2);
  0 0 traj.rddd(3,k2);
  0 0 traj.rdddd(3,k2);
  0 0 traj.rddddd(3,k2);
  0 0 traj.rdddddd(3,k2)
  ];
[rz, ~] = TransitionTrajectory(Z, T, Ts);
traj.r(:,k0:k2)        = [ rx(1,:); ry(1,:); rz(1,:); ];
traj.rd(:,k0:k2)       = [ rx(2,:); ry(2,:); rz(2,:); ];
traj.rdd(:,k0:k2)      = [ rx(3,:); ry(3,:); rz(3,:); ];
traj.rddd(:,k0:k2)     = [ rx(4,:); ry(4,:); rz(4,:); ];
traj.rdddd(:,k0:k2)    = [ rx(5,:); ry(5,:); rz(5,:); ];
traj.rddddd(:,k0:k2)   = [ rx(6,:); ry(6,:); rz(6,:); ];
traj.rdddddd(:,k0:k2)  = [ rx(7,:); ry(7,:); rz(7,:); ];

% transition from circle to origin
T = [t3 t4 t5];
X = [
  traj.r(1,k3) 0 0;
  traj.rd(1,k3) 0 0;
  traj.rdd(1,k3) 0 0;
  traj.rddd(1,k3) 0 0;
  traj.rdddd(1,k3) 0 0;
  traj.rddddd(1,k3) 0 0;
  traj.rdddddd(1,k3) 0 0
  ];
[rx, ~] = TransitionTrajectory(X, T, Ts);
Y = [
  traj.r(2,k3) 0 0;
  traj.rd(2,k3) 0 0;
  traj.rdd(2,k3) 0 0;
  traj.rddd(2,k3) 0 0;
  traj.rdddd(2,k3) 0 0;
  traj.rddddd(2,k3) 0 0;
  traj.rdddddd(2,k3) 0 0
  ];
[ry, ~] = TransitionTrajectory(Y, T, Ts);
Z = [
  traj.r(3,k3) 0 0;
  traj.rd(3,k3) 0 0;
  traj.rdd(3,k3) 0 0;
  traj.rddd(3,k3) 0 0;
  traj.rdddd(3,k3) 0 0;
  traj.rddddd(3,k3) 0 0;
  traj.rdddddd(3,k3) 0 0
  ];
[rz, ~] = TransitionTrajectory(Z, T, Ts);
traj.r(:,k3:k5)        = [ rx(1,:); ry(1,:); rz(1,:); ];
traj.rd(:,k3:k5)       = [ rx(2,:); ry(2,:); rz(2,:); ];
traj.rdd(:,k3:k5)      = [ rx(3,:); ry(3,:); rz(3,:); ];
traj.rddd(:,k3:k5)     = [ rx(4,:); ry(4,:); rz(4,:); ];
traj.rdddd(:,k3:k5)    = [ rx(5,:); ry(5,:); rz(5,:); ];
traj.rddddd(:,k3:k5)   = [ rx(6,:); ry(6,:); rz(6,:); ];
traj.rdddddd(:,k3:k5)  = [ rx(7,:); ry(7,:); rz(7,:); ];

% compute corresponding orientation
traj = QuadTraj(param, traj, optHeading);

end

