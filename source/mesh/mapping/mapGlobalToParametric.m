function [ parametricCoordinate ] = mapGlobalToParametric( x, x0, xend )
%MAPGLOBALTOPARAMETRIC map global coordinates onto the parametric space
%[0,1].
%Input:
%x = global coordinate
%x0 = left end of the domain global coordinate
%xend = right end of the domain global coordinate
%Output:
%parametricCoordinate = parametric coordinate

parametricCoordinate = (x - x0) / (xend - x0);

end

