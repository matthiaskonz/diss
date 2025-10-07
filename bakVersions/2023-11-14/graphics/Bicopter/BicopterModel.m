function zd = BicopterModel(t, z, param, traj)

% controller
u = BicopterCtrl(t, z, param, traj);

% model param
m = param.m;
Jx = param.Jx;
Jy = param.Jy;
Jz = param.Jz;
g = param.g;
lY = param.lY;
lZ = param.lZ;
betaF = param.betaF;
bF = param.bF;

L = 10*eye(7);

% model state
x = z(1:12);
xi = z(13:18);

% model matrices
A = [x(4) x(7) x(10) 0 0 0; x(5) x(8) x(11) 0 0 0; x(6) x(9) x(12) 0 0 0; 0 0 0 0 -x(10) x(7); 0 0 0 0 -x(11) x(8); 0 0 0 0 -x(12) x(9); 0 0 0 x(10) 0 -x(4); 0 0 0 x(11) 0 -x(5); 0 0 0 x(12) 0 -x(6); 0 0 0 -x(7) x(4) 0; 0 0 0 -x(8) x(5) 0; 0 0 0 -x(9) x(6) 0;];
M = [m 0 0 0 0 0; 0 m 0 0 0 0; 0 0 m 0 0 0; 0 0 0 Jx 0 0; 0 0 0 0 Jy 0; 0 0 0 0 0 Jz;];
f = [g * m * x(6) - m * xi(2) * xi(6) + m * xi(3) * xi(5); g * m * x(9) + m * xi(1) * xi(6) - m * xi(3) * xi(4); g * m * x(12) - m * xi(1) * xi(5) + m * xi(2) * xi(4); -Jy * xi(5) * xi(6) + Jz * xi(5) * xi(6); Jx * xi(4) * xi(6) - Jz * xi(4) * xi(6); -Jx * xi(4) * xi(5) + Jy * xi(4) * xi(5);];
B = [1 0 1 0; 0 sin(betaF) 0 -sin(betaF); 0 cos(betaF) 0 cos(betaF); bF -sin(betaF) * lZ + cos(betaF) * lY -bF sin(betaF) * lZ - cos(betaF) * lY; lZ sin(betaF) * bF lZ sin(betaF) * bF; -lY cos(betaF) * bF lY -cos(betaF) * bF;];
phi = [x(4) ^ 2 + x(5) ^ 2 + x(6) ^ 2 - 1; x(4) * x(7) + x(5) * x(8) + x(6) * x(9); x(4) * x(10) + x(5) * x(11) + x(6) * x(12); x(7) ^ 2 + x(8) ^ 2 + x(9) ^ 2 - 1; x(7) * x(10) + x(8) * x(11) + x(9) * x(12); x(10) ^ 2 + x(11) ^ 2 + x(12) ^ 2 - 1; x(4) * x(8) * x(12) - x(4) * x(9) * x(11) - x(5) * x(7) * x(12) + x(5) * x(9) * x(10) + x(6) * x(7) * x(11) - x(6) * x(8) * x(10) - 1;];
Phi = [0 0 0 2 * x(4) 2 * x(5) 2 * x(6) 0 0 0 0 0 0; 0 0 0 x(7) x(8) x(9) x(4) x(5) x(6) 0 0 0; 0 0 0 x(10) x(11) x(12) 0 0 0 x(4) x(5) x(6); 0 0 0 0 0 0 2 * x(7) 2 * x(8) 2 * x(9) 0 0 0; 0 0 0 0 0 0 x(10) x(11) x(12) x(7) x(8) x(9); 0 0 0 0 0 0 0 0 0 2 * x(10) 2 * x(11) 2 * x(12); 0 0 0 x(8) * x(12) - x(9) * x(11) -x(7) * x(12) + x(9) * x(10) x(7) * x(11) - x(8) * x(10) -x(5) * x(12) + x(6) * x(11) x(4) * x(12) - x(6) * x(10) -x(4) * x(11) + x(5) * x(10) x(5) * x(9) - x(6) * x(8) -x(4) * x(9) + x(6) * x(7) x(4) * x(8) - x(5) * x(7);];

% EoM
xd = A*xi - pinv(Phi)*L*phi;
xid = M \ (B*u - f);

zd = [xd; xid];

[v,i] = max(abs(phi));
if (v > 1e-3)
  fprintf('WARNING: Geometric constraint #%i has error=%e at time %f\n', i, v, t);
end

end


