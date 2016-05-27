function [ KE ] = rbNonLinearConductivityMatrix( problem, solution, iMode,...
    enrichedElementCoords, e)
%RBNONLINEARCONDUCTIVITYMATRIX computes the conductivity local to one element
%   problem = definition of the boundary value problem
%   solution = actual temperature value
%   iModes = index of the ith POD mode
%   enrichedElementCoords = coordinates of the refined element
%   e = element index

X1 = problem.coords(e);
X2 = problem.coords(e+1);
ldof = 2;

KE = zeros(ldof, ldof);

for i=1:ldof
    for j=1:ldof
        KE(i,j) = problem.rbB(@(x)problem.rbBasis_fun(x, i, iMode, 1.0, problem, enrichedElementCoords ),...
            @(x)problem.rbBasis_fun(x, j, iMode, 1.0, problem, enrichedElementCoords),...
            @(x, d)evaluateActualTemperature(x, problem, solution, d), X1, X2);
    end
end

end

