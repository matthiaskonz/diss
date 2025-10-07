function dz = ModelPlanarRigidBody(t, z, param, traj)

% controller
u = CtrlPlanarRigidBody(t, z, param, traj);

% param
m = param.m;
sx = param.sx;
sy = param.sy;
Jzz = param.Jzz;

Lambda = 1;

% state
x = z(1:4);
xi = z(5:7);

A = [x(4) -x(3) 0; x(3) x(4) 0; 0 0 x(4); 0 0 -x(3);];
M = [m 0 -m * sy; 0 m m * sx; -m * sy m * sx Jzz;];
f = [-m * sx * xi(3) ^ 2 - m * xi(2) * xi(3); -m * sy * xi(3) ^ 2 + m * xi(1) * xi(3); m * sx * xi(1) * xi(3) + m * sy * xi(2) * xi(3);];
B = [1 0 0; 0 1 0; 0 0 1;];
phi = [x(3) ^ 2 + x(4) ^ 2 - 1;];
Phi = [0 0 2 * x(3) 2 * x(4);];

% kinematic equation
xd = A*xi - pinv(Phi)*Lambda*phi;
% kinetic equation
xid = -M \ (f - B*u);

dz = [xd; xid];

[maxErr, idx] = max(abs(phi));
if (maxErr > 1e-3)
  % should not occur in theory - only due to numeric errors or wrong initial conditions
  fprintf('WARNING: Geometric constraint %i not fulfilled at time=%f, error=%e \n', idx, t, maxErr);
end

end

