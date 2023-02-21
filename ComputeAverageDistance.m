function [ E ] = ComputeAverageDistance( R, min_dist, alpha )

% Parameters
N = 1e6;
d = linspace(min_dist,2*R,N);

% Variables
p = 4/(pi*R^2) * d .* (acos(d/(2*R)) - d/(2*R) .* sqrt(1 - d.^2/(4*R^2)));
p = p/sum(p);
l = d.^(-alpha);

% Output
E = [sum(p.*d) sum(p.*l)];

end

