function [u, Wc, Pc] = CtrlRevoluteJoint(t, z, param)

% model param
Jz = param.Jz;
% closed loop param
Jcz = param.Jcz;
sigcz = param.sigcz;
kapcz = param.kapcz;

% model state
a = z(1);
w = z(2);
% reference traj
aR  = interp1(param.tR, param.xR', t)';
wR  = interp1(param.tR, param.xiR', t)';
wRd = interp1(param.tR, param.xiRd', t)';

% control law
if (param.CtrlMode == 0)
  u = Jz/Jcz * (Jcz*wRd - sigcz*(w-wR) - kapcz*(a-aR));
  Tc = .5*Jcz*(w-wR)^2;
  Pc = .5*sigcz*(w-wR)^2;
  Vc = .5*kapcz*(a-aR)^2;
elseif (param.CtrlMode == 1)
  u = -Jz / Jcz * (-Jcz * wRd + sigcz * (w - wR) - kapcz * (-cos(aR) * sin(a) + sin(aR) * cos(a)));
  Tc = (w / 0.2e1 - wR / 0.2e1) * Jcz * (w - wR);
  Pc = (w / 0.2e1 - wR / 0.2e1) * sigcz * (w - wR);
  Vc = -kapcz * cos(-aR + a) + kapcz;
elseif (param.CtrlMode == 2)
  u = -Jz / Jcz * (-Jcz * (cos(aR) * sin(a) * wR ^ 2 - sin(aR) * cos(a) * wR ^ 2 + cos(aR) * cos(a) * wRd + sin(aR) * sin(a) * wRd) - sigcz * (0.2e1 * cos(aR) * cos(a) * wR + 0.2e1 * sin(aR) * sin(a) * wR - (2 * w)) / 0.2e1 - kapcz * (-cos(aR) * sin(a) + sin(aR) * cos(a)));
  Tc = (-w * wR * cos(-aR + a) + w ^ 2 / 0.2e1 + wR ^ 2 / 0.2e1) * Jcz;
  Pc = (-w * wR * cos(-aR + a) + w ^ 2 / 0.2e1 + wR ^ 2 / 0.2e1) * sigcz;
  Vc = (-cos(-aR + a) + 0.1e1) * kapcz;
end
Wc = Tc + Vc;


end


