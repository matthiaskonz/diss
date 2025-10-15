function [ x_t, w_t, wd_t, wdd_t, wddd_t ] = VecRotVel2RotMat( t, R0, b3d_t, b3dd_t, b3ddd_t, b3dddd_t, w3_t, w3d_t, w3dd_t, w3ddd_t )

l = 1/(t(2)-t(1)); % param for numeric stabilization

N = length(t);

x_t     = zeros(9,N);
w_t     = zeros(3,N);
wd_t    = zeros(3,N);
wdd_t   = zeros(3,N);
wddd_t  = zeros(3,N);

R = R0;
x_t(:,1) = [R(:,1); R(:,2); R(:,3)];

for k = 1:N-1

  Rd    = [zeros(3,2) b3d_t(:,k)];
  Rdd   = [zeros(3,2) b3dd_t(:,k)];
  Rddd  = [zeros(3,2) b3ddd_t(:,k)];
  Rdddd = [zeros(3,2) b3dddd_t(:,k)];
  
  w     = [0; 0; w3_t(k)];
  wd    = [0; 0; w3d_t(k)];
  wdd   = [0; 0; w3dd_t(k)];
  wddd  = [0; 0; w3ddd_t(k)];
  
  w(1) = -R(:,2)' * Rd(:,3);
  w(2) =  R(:,1)' * Rd(:,3);
  W = wed_so3(w);
  Rd(:,1:2) = R * W(:,1:2);

  wd(1) = -Rd(:,2)' * Rd(:,3) - R(:,2)' * Rdd(:,3);
  wd(2) =  Rd(:,1)' * Rd(:,3) + R(:,1)' * Rdd(:,3);
  Wd = wed_so3(wd);
  Rdd(:,1:2) = Rd * W(:,1:2) + R * Wd(:,1:2);

  wdd(1) = -Rdd(:,2)' * Rd(:,3) - 2*Rd(:,2)' * Rdd(:,3) - R(:,2)' * Rddd(:,3);
  wdd(2) =  Rdd(:,1)' * Rd(:,3) + 2*Rd(:,1)' * Rdd(:,3) + R(:,1)' * Rddd(:,3);
  Wdd = wed_so3(wdd);
  Rddd(:,1:2) = Rdd * W(:,1:2) + 2*Rd * Wd(:,1:2) + R * Wdd(:,1:2);
 
  wddd(1) = -Rddd(:,2)' * Rd(:,3) - 3*Rdd(:,2)' * Rdd(:,3) - 3*Rd(:,2)' * Rddd(:,3) - R(:,2)' * Rdddd(:,3);
  wddd(2) =  Rddd(:,1)' * Rd(:,3) + 3*Rdd(:,1)' * Rdd(:,3) + 3*Rd(:,1)' * Rddd(:,3) + R(:,1)' * Rdddd(:,3);
  % Wddd = wed_so3(wddd);
  % Rddd(:,1:2) = Rddd * W(:,1:2) + 3*Rdd * Wd(:,1:2) + 3*Rd * Wdd(:,1:2) + R * Wddd(:,1:2);
  
  w_t(:,k)     = w;
  wd_t(:,k)    = wd;
  wdd_t(:,k)   = wdd;
  wddd_t(:,k)  = wddd;
  
  % integrate rotation matrix
  R11 = R(1,1);
  R21 = R(2,1);
  R31 = R(3,1);
  R12 = R(1,2);
  R22 = R(2,2);
  R32 = R(3,2);
  R13 = R(1,3);
  R23 = R(2,3);
  R33 = R(3,3);
  Asq = [0 -R13 R12 R11 / 0.2e1 0 0 R12 / 0.2e1 R13 / 0.2e1 0; 0 -R23 R22 R21 / 0.2e1 0 0 R22 / 0.2e1 R23 / 0.2e1 0; 0 -R33 R32 R31 / 0.2e1 0 0 R32 / 0.2e1 R33 / 0.2e1 0; R13 0 -R11 0 R12 / 0.2e1 0 R11 / 0.2e1 0 R13 / 0.2e1; R23 0 -R21 0 R22 / 0.2e1 0 R21 / 0.2e1 0 R23 / 0.2e1; R33 0 -R31 0 R32 / 0.2e1 0 R31 / 0.2e1 0 R33 / 0.2e1; -R12 R11 0 0 0 R13 / 0.2e1 0 R11 / 0.2e1 R12 / 0.2e1; -R22 R21 0 0 0 R23 / 0.2e1 0 R21 / 0.2e1 R22 / 0.2e1; -R32 R31 0 0 0 R33 / 0.2e1 0 R31 / 0.2e1 R32 / 0.2e1;];
  phi = [R11 ^ 2 + R21 ^ 2 + R31 ^ 2 - 1; R12 ^ 2 + R22 ^ 2 + R32 ^ 2 - 1; R13 ^ 2 + R23 ^ 2 + R33 ^ 2 - 1; R11 * R12 + R21 * R22 + R31 * R32; R11 * R13 + R21 * R23 + R31 * R33; R12 * R13 + R22 * R23 + R32 * R33;];
  
  x_t(:,k+1) = x_t(:,k) + (t(k+1)-t(k)) * Asq * [w; -l*phi ];
  R = [x_t(1:3,k+1) x_t(4:6,k+1) x_t(7:9,k+1)];
  
end

w_t(:,N)     = w;
wd_t(:,N)    = wd;
wdd_t(:,N)   = wdd;
wddd_t(:,N)  = wddd;

end

