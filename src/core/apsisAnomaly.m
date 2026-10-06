function nu = apsisAnomaly(letter)
%APSISANOMALY True anomaly of an apsis: 'p' (periapsis) -> 0, 'a' (apoapsis) -> pi.

switch lower(letter)
    case 'p', nu = 0;
    case 'a', nu = pi;
    otherwise, error('apsisAnomaly:badLetter', 'Use ''p'' or ''a''.');
end
end
