function [ I ] = ComputeRealIntegralSG_MISO( radius, aoa, distance, angle, alpha, s )

% Parameters
N       = 1e3;
theta   = linspace(aoa - angle/2,aoa + angle/2,N);
ro      = distance;

% Variables
xmin        = zeros(1,N);
xmax        = zeros(1,N);
fun         = zeros(1,N);
for n   = 1:N
    p       = [1, -2*ro*cos(theta(n) - aoa), ro^2 - radius^2];
    r       = abs(roots(p));
    xmin(n) = min(r);
    xmax(n) = max(r);
    fun(n)  = ComputeIntegralSG(xmin(n),xmax(n),alpha)*s(theta(n));
end

% Output
d       = diff(theta);
if d(1) == 1/N
    I   = abs(mean(fun)); 
else
    I   = d(1)*sum(fun);
end

end