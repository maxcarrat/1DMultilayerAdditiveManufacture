function [ knotVector ] = getOpenKnotVector( layer, p )
%GETOPENKNOTVECTOR generate the knot vector given the layer index and the
%polynomial order


knotVector = zeros( 1, 2 * (p+1) + layer - 1);

for i=1:p
    knotVector(i) = 0.0;
    knotVector(end-(i-1)) = 1.0;
end

knotVector(p+1:end-p) = linspace(0.0, 1.0, layer+1);

end

