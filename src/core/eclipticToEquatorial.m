function vEq = eclipticToEquatorial(vEcl, obliquity)
%ECLIPTICTOEQUATORIAL Rotate a vector from the ecliptic to the equatorial frame.
%   vEq = ECLIPTICTOEQUATORIAL(vEcl, obliquity) rotates about the common
%   x axis (vernal equinox) by the obliquity of the ecliptic [rad].

c = cos(obliquity);
s = sin(obliquity);
vEq = [1 0 0; 0 c -s; 0 s c]*vEcl(:);
end
