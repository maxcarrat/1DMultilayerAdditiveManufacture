function [ projectedTemperature ] = L2projection(problem, previousTemperature, updatedMesh, previousMesh)
%L2PROJECTION project the previous solution onto the updated mesh at the
%new time step
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   previousmesh = previous mesh configuration

projectedTemperature = zeros(size(updatedMesh,2), 1);

for j=1:size(updatedMesh,2)
    x = updatedMesh(j);
    projectedTemperature(j) = globalProjection(x, previousMesh, previousTemperature, problem);
end

end


function [ numericalSolutions ] = globalProjection(x, previousMesh, coefficients, problem )
% numericalSolutions = GLOBALPROJECTION(postProcessingCoords,, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations

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
    ldof = 2;

    r = zeros(size(x));
    for i=1:ldof
        r=r+coefficients(problem.LM(element,i)).* problem.basis_fun(mapGlobalToLocal(x, X1, X2), i, 0);
    end

    
end


