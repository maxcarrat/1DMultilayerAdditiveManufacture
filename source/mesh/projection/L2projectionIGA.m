function [ projectedTemperature ] = L2projectionIGA(problem, previousTemperature, controlPoints, previousKnotVector,...
    initialTemperature)
%L2PROJECTIONLINEARDISTRIBUTION project the previous solution onto the updated mesh at the
%new time step
%   previousTemperature = temeprature distribution of the previous mesh
%   controlPoint = control points at new time step
%   previousKnotVector = previous knot vector
%   initialTemperature = initial temperature of the powder

projectedTemperature = zeros(size(controlPoints, 2), 1);

for j=1:size(controlPoints, 2)
    
    x = linspace(0, 1, length(controlPoints));
    globalProjectedValue = globalProjection(x, previousKnotVector, previousTemperature, problem);
    
    if globalProjectedValue(j) ~= 0.0
        projectedTemperature(j) = globalProjectedValue(j);
    else
        projectedTemperature(j) = (previousTemperature(end))/(2^problem.refinementDepth)...
            * (-j + (size(problem.coords,2))) + initialTemperature;
    end
end

end


function [ numericalSolutions ] = globalProjection(x, previousMesh, coefficients, problem)
% numericalSolutions = GLOBALPROJECTION(x, previousMesh, coefficients, problem) evaluates the numerical solution
% x = coordinates to post process
% previousMesh = mesh onto whom I project
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% problem = transient poisson problem struct

previousKnotVector = previousMesh;
numericalSolutions=zeros(size(x));
Xi1 = previousKnotVector(1+problem.p);
Xi2 = previousKnotVector(2+problem.p);

numericalSolutions(x>=Xi1 & x<=Xi2) = localProjection(x(x>=Xi1 & x<=Xi2), 1,...
    coefficients, problem);

for e=2:problem.N
    Xi1 = previousKnotVector(e + problem.p);
    Xi2 = previousKnotVector(e + 1 + problem.p);
    
    numericalSolutions(x>Xi1 & x<=Xi2) = localProjection(x(x>Xi1 & x<=Xi2), e,...
        coefficients, problem);
end

end

function r = localProjection(x, element, coefficients, problem)
% r = LOCALPROJECTION(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   problem = IGA physical problem struct

% Xi1 = problem.knotVector( element + problem.p );
% Xi2 = problem.knotVector( element + problem.p + 1);

r = zeros(size(x));

for k=1:length(x)
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(x(k), problem.p, problem.knotVector);
    
    if isempty(r) == 0
        r(k) = r(k) + N(problem.LM(element,:)) * coefficients(problem.LM(element,:));
    end
    
end

end



