function [dv1, dv2, dt, argpF, aT, eT] = bitangentTransfer(aI, eI, aF, eF, type, argpI, mu)
%BITANGENTTRANSFER Two-impulse bitangent transfer between coaxial ellipses.
%   [dv1, dv2, dt, argpF, aT, eT] = BITANGENTTRANSFER(aI, eI, aF, eF, type,
%   argpI, mu) transfers from the orbit (aI, eI) to the orbit (aF, eF).
%   type = 'pa', 'ap', 'pp' or 'aa' selects the departure and arrival apses
%   (p = periapsis, a = apoapsis). For 'pp' and 'aa' the apse line of the
%   final orbit is rotated by pi, which is returned in argpF.
%   dv1, dv2 [km/s] are signed tangential impulses, dt [s] is the transfer
%   time (half period of the transfer ellipse aT, eT).

switch lower(type)
    case 'pa', r1 = aI*(1 - eI); r2 = aF*(1 + eF); dArgp = 0;
    case 'ap', r1 = aI*(1 + eI); r2 = aF*(1 - eF); dArgp = 0;
    case 'pp', r1 = aI*(1 - eI); r2 = aF*(1 - eF); dArgp = pi;
    case 'aa', r1 = aI*(1 + eI); r2 = aF*(1 + eF); dArgp = pi;
    otherwise, error('bitangentTransfer:badType', 'type must be pa, ap, pp or aa.');
end
aT = (r1 + r2)/2;
eT = abs(r2 - r1)/(r1 + r2);
dv1 = sqrt(mu)*(sqrt(2/r1 - 1/aT) - sqrt(2/r1 - 1/aI));
dv2 = sqrt(mu)*(sqrt(2/r2 - 1/aF) - sqrt(2/r2 - 1/aT));
argpF = mod(argpI + dArgp, 2*pi);
dt = pi*sqrt(aT^3/mu);
end
