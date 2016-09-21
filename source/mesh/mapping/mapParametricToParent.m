function [ x ] = mapParametricToParent( X, X1, X2 )
%MAPPARAMETRICTOPARENT map parameter space [0, 1] onto the 
%parent space [-1, 1]
% inputs:
% X = parametric coordinates
% X1 = left knot location
% X2 = right knot location
% outputs:
% x = parent coordinates


x = (2.0*X - (X2+X1)) / (X2-X1);

end


