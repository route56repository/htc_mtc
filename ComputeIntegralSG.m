function [ I ] = ComputeIntegralSG( a, b, alpha )

if alpha == 2
    I = log(b/a);
else                           
    I = 1/(2 - alpha) * (b^(2 - alpha) - a^(2 - alpha)); 
end

end