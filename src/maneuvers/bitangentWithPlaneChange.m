function m = bitangentWithPlaneChange(oI, oF, type, mu)
%BITANGENTWITHPLANECHANGE Bitangent transfer with the plane change inside the arc.
%   m = BITANGENTWITHPLANECHANGE(oI, oF, type, mu) performs the bitangent
%   transfer of type 'pa', 'ap', 'pp' or 'aa' from orbit oI to the shape of
%   orbit oF (structs with a, e, i, raan, argp) and moves the plane change
%   onto the transfer arc, at the node between the transfer plane and the
%   final plane, where the speed is lower than on the initial orbit.
%
%   Output struct m:
%     dv1, dv2      tangential impulses at the two apses [km/s]
%     dvPlane       plane-change impulse on the transfer arc [km/s]
%     argpBefore    argument of periapsis of the transfer before the plane change
%     argpAfter     ... and after it [rad]
%     argpFinal     argument of periapsis on the final orbit after dv2 [rad]
%     aT, eT        transfer ellipse
%     nuNode        true anomaly of the plane change on the transfer [rad]
%     dtCoast       [coast to the node, coast from the node to arrival] [s]
%     feasible      false if neither node lies on the travelled arc
%
%   Assumption: for pa/pp/aa the departure apsis is lower than the arrival
%   apsis (transfer from periapsis to apoapsis), for ap it is higher.

switch lower(type)
    case 'pa'
        rDep = oI.a*(1 - oI.e); rArr = oF.a*(1 + oF.e);
        nuDep = 0;  nuArr = pi;   argpBefore = oI.argp;              dArgpFinal = 0;
    case 'ap'
        rDep = oI.a*(1 + oI.e); rArr = oF.a*(1 - oF.e);
        nuDep = pi; nuArr = 2*pi; argpBefore = oI.argp;              dArgpFinal = 0;
    case 'pp'
        rDep = oI.a*(1 - oI.e); rArr = oF.a*(1 - oF.e);
        nuDep = 0;  nuArr = pi;   argpBefore = oI.argp;              dArgpFinal = pi;
    case 'aa'
        rDep = oI.a*(1 + oI.e); rArr = oF.a*(1 + oF.e);
        nuDep = 0;  nuArr = pi;   argpBefore = mod(oI.argp + pi, 2*pi); dArgpFinal = 0;
    otherwise
        error('bitangentWithPlaneChange:badType', 'type must be pa, ap, pp or aa.');
end
if xor(strcmpi(type, 'ap'), rDep > rArr)
    error('bitangentWithPlaneChange:apsisOrder', ...
        'Unsupported apsis ordering for type %s (see help).', type);
end

m = struct('dv1', NaN, 'dv2', NaN, 'dvPlane', NaN, 'argpBefore', argpBefore, ...
    'argpAfter', NaN, 'argpFinal', NaN, 'aT', NaN, 'eT', NaN, 'nuNode', NaN, ...
    'dtCoast', [NaN NaN], 'feasible', false);

m.aT = (rDep + rArr)/2;
m.eT = abs(rArr - rDep)/(rArr + rDep);
m.dv1 = sqrt(mu)*(sqrt(2/rDep - 1/m.aT) - sqrt(2/rDep - 1/oI.a));
m.dv2 = sqrt(mu)*(sqrt(2/rArr - 1/oF.a) - sqrt(2/rArr - 1/m.aT));

% Plane change applied to the transfer ellipse
[m.dvPlane, m.argpAfter, m.nuNode] = planeChange(m.aT, m.eT, oI.i, oI.raan, argpBefore, oF.i, oF.raan, mu);

if ~isOnArc(m.nuNode, nuDep, nuArr)
    nuOther = mod(m.nuNode + pi, 2*pi);      % try the other node
    if ~isOnArc(nuOther, nuDep, nuArr)
        return
    end
    m.nuNode = nuOther;
    p = m.aT*(1 - m.eT^2);
    vNode = sqrt(mu/p)*(1 + m.eT*cos(m.nuNode));
    alpha = acos(max(-1, min(1, cos(oI.i)*cos(oF.i) + sin(oI.i)*sin(oF.i)*cos(oF.raan - oI.raan))));
    m.dvPlane = 2*vNode*sin(alpha/2);
end

m.feasible = true;
m.dtCoast = [timeOfFlight(m.aT, m.eT, nuDep, m.nuNode, mu), ...
             timeOfFlight(m.aT, m.eT, m.nuNode, nuArr, mu)];
m.argpFinal = mod(m.argpAfter + dArgpFinal, 2*pi);
end

function tf = isOnArc(nu, nu0, nu1)
% true if nu lies on the arc travelled from nu0 to nu1
nu = mod(nu, 2*pi); nu0 = mod(nu0, 2*pi); nu1 = mod(nu1, 2*pi);
if nu0 <= nu1
    tf = nu >= nu0 && nu <= nu1;
else
    tf = nu >= nu0 || nu <= nu1;
end
end
