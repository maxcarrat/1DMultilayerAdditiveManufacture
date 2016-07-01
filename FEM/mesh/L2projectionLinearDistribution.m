function [ projectedTemperature ] = L2projectionLinearDistribution(problem, previousTemperature, updatedMesh, previousMesh, initialTemperature, refinementDepth)
%L2PROJECTIONLINEARDISTRIBUTION project the previous solution onto the updated mesh at the
%new time step
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   previousmesh = previous mesh configuration
%   initialTemperature = initial temperature of the powder
%   refinementDepth = refinement dpeth of the last element

projectedTemperature = zeros(size(updatedMesh,2), 1);

for j=1:size(updatedMesh,2)
    x = updatedMesh(j);
    globalProjectedValue = globalProjection(x, previousMesh, previousTemperature, problem);
    if globalProjectedValue ~= 0.0
        projectedTemperature(j) = globalProjectedValue;
    else
        projectedTemperature(j) = (previousTemperature(end))/(2^refinementDepth)...
            * (-j + (size(updatedMesh,2))) + initialTemperature;
    end
end

end


function [ numericalSolutions ] = globalProjection(x, previousMesh, coefficients, problem )
% numericalSolutions = GLOBALPROJECTION(x, previousMesh, coefficients, problem) evaluates the numerical solution
% x = coordinates to post process
% previousMesh = mesh onto whom I project 
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% problem = transient poisson problem struct

    coords = previousMesh;
    numericalSolutions=zeros(size(x));
    X1 = coords(1);
    X2 = coords(2);

    numericalSolutions(x>=X1 & x<=X2) = localProjection(x(x>=X1 & x<=X2), coords, 1, coefficients, problem);
    
    for e=2:size(previousMesh, 2)-1
        X1 = coords(e);
        X2 = coords(e+1);
        
        numericalSolutions(x>X1 & x<=X2) = localProjection(x(x>X1 & x<=X2), coords, e, coefficients, problem);
    end
end

function r = localProjection(x, coords, element, coefficients, problem)
% r = LOCALPROJECTION(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations

    X1 = coords(element);
    X2 = coords(element+1);
    
    [N, ~] = shapeFunctionsAndDerivatives(mapGlobalToLocal(x, X1, X2));

    r = zeros(size(x));
    
    if isempty(r) == 0
        r = r + N * [coefficients(problem.LM(element,1)) coefficients(problem.LM(element,2))]';
    end
end


