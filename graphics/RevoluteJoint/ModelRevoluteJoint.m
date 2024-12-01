function zd = ModelRevoluteJoint(t, z, param)

% controller
u = CtrlRevoluteJoint(t, z, param);

a = z(1);
w = z(2);

ad = w;
wd = u/param.Jz;

zd = [ad; wd];

end

