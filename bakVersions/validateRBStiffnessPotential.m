clear;

N = 50;

cc = rand(1,N);
hh = rand(3,N);
pp = rand(3,N);

c = sum(cc);
h = (hh*cc')/c;
p = (pp*cc')/c;
P = (cc.*hh)*pp';
Ps = P - c*h*p';
[R0, Ksp] = specialPolarDecomp(Ps');
r0 = p - R0*h;

r = rand(3,1);
R = quat2rotm(unitVec(rand(4,1)));
e = r + R*h-(r0+R0*h);
V01 = .5*c*(e'*e) + trace(Ksp*(eye(3) - R0'*R));

V = 0;
V0 = 0;
dV0 = zeros(6,1);
for k = 1:N
    ck = cc(k);
    hk = hh(:,k);
    pk = pp(:,k);
    ek = r + R*hk - pk;
    V = V + .5*ck*(ek'*ek);
    ek0 = r0 + R0*hk - pk;
    V0 = V0 + .5*ck*(ek0'*ek0);
    dV0 = dV0 + ck*[R0'*ek0; -wed_so3(hk)*R0'*ek0];
end

V0 + V01 - V
dV0