function Q = perifocalToInertial(raan, i, argp)
%PERIFOCALTOINERTIAL Rotation matrix from the perifocal to the inertial frame.
%   Q = PERIFOCALTOINERTIAL(raan, i, argp) = R3(raan)*R1(i)*R3(argp) [rad].

cO = cos(raan); sO = sin(raan);
ci = cos(i);    si = sin(i);
cw = cos(argp); sw = sin(argp);
R3raan = [cO -sO 0; sO cO 0; 0 0 1];
R1i    = [1 0 0; 0 ci -si; 0 si ci];
R3argp = [cw -sw 0; sw cw 0; 0 0 1];
Q = R3raan*R1i*R3argp;
end
