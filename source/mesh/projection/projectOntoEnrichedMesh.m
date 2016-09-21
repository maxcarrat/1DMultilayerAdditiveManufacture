function [ projectedTemperature ] = projectOntoEnrichedMesh( problem, temperatureCoefficients,...
    modes, refinedMesh,  previousMesh, PODRefinementDepth, initialTemperature )
%PROJECTONTOENRICHEDMESH 

projectedTemperature = zeros(numel(refinedMesh) + (2.^PODRefinementDepth)*modes, 1);

for j=1:numel(refinedMesh)
    x = refinedMesh(j);
    globalProjectedValue = globalProjection(x, previousMesh, temperatureCoefficients, problem);
    if globalProjectedValue ~= 0.0
        projectedTemperature(j) = globalProjectedValue;
    else
        projectedTemperature(j) = initialTemperature;
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
    numericalSolutions = zeros(size(x));
    X1 = coords(1);
    X2 = coords(2);

    numericalSolutions(x>=X1 & x<=X2) = localProjection(x(x>=X1 & x<=X2), coords, 1, coefficients, problem);
    
    for e=2:numel(previousMesh)-1
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