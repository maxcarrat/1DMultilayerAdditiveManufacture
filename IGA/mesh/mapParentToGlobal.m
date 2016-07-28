function [ X ] = mapParentToGlobal( x, Xi1, Xi2, problem, e )
%MAPPARENTTOGLOBAL map quadrature/parent space [-1, 1] onto the 
%Global space 
% inputs:
% x = parent coordinates
% Xi1 = knot coord left
% Xi2 = knot coord right
% problem = IGS physical problem struct
% e = element index
% outputs:
% X = global coordinates

localCoords = mapParentToLocal(x, Xi1, Xi2);

i = findspan(length(problem.coords)-1, problem.p, localCoords, problem.knotVector);
X = 0.0;

[N, ~] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);

for j=1:problem.p+1
    
    X = X + N((e-1)+j) * problem.coords(i - problem.p + j);
    
end

end

