function dz = ModelScara(t, z, param, traj)

% controller
u = CtrlScara(t, z, param, traj);

% model param
l1  = param.l1;
l2  = param.l2;
J1z = param.J1z;
m2  = param.m2;
J2z = param.J2z;
s2x = param.s2x;
s2y = param.s2y;

% model state
x = z(1:2);
xi = z(3:4);

% model matrices
A = [1 0; 0 1;];
M = [l1 ^ 2 * m2 + (0.2e1 * cos(x(2)) * s2x - 0.2e1 * sin(x(2)) * s2y) * m2 * l1 + J1z + J2z (cos(x(2)) * s2x - sin(x(2)) * s2y) * m2 * l1 + J2z; (cos(x(2)) * s2x - sin(x(2)) * s2y) * m2 * l1 + J2z J2z;];
f = [(-xi(2) * (2 * xi(1) + xi(2)) * s2x * sin(x(2)) - xi(2) * (2 * xi(1) + xi(2)) * s2y * cos(x(2))) * m2 * l1; ((xi(1) ^ 2) * s2x * sin(x(2)) + (xi(1) ^ 2) * cos(x(2)) * s2y) * m2 * l1;];
B = [1 0; 0 1;];

% kinematic equation
xd = A*xi;
% kinetic equation
xid = -M \ (f - B*u);

dz = [xd; xid];

end

