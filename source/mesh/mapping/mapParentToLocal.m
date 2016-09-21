function [ X ] = mapParentToLocal( x, X1, X2 )
%MAPPARENTTOLOCAL map quadrature/parent space [-1, 1] onto the 
%BSpline local space [0, 1]
% inputs:
% x = local coordinates
% X1 = left knot location
% X2 = right knot location
% outputs:
% X = parent coordinates

X = 0.5*((X2-X1)*x + (X2+X1));

end

