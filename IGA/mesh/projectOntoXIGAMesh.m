function [ projectedTemperature ] = projectOntoXIGAMesh( problem, temperatureCoefficients,...
    modes, knotVector,  previousKnotVector, PODRefinementDepth, initialTemperature )
%PROJECTONTOXIGAMESH

projectedTemperature = zeros(numel(knotVector) + (2.^PODRefinementDepth)*modes, 1);

for j=1:numel(problem.coords)
    x = problem.coords(j);
    globalProjectedValue = globalProjection(x, previousKnotVector, temperatureCoefficients, problem);
    if globalProjectedValue ~= 0.0
        projectedTemperature(j) = globalProjectedValue;
    else
        projectedTemperature(j) = initialTemperature;
    end
end

end

function [ numericalSolutions ] = globalProjection(x, previousKnotVector, coefficients, problem )
% numericalSolutions = GLOBALPROJECTION(x, previousMesh, coefficients, problem) evaluates the numerical solution
% x = coordinates to post process
% previousKnotVector = mesh onto whom I project
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% problem = transient poisson problem struct

numericalSolutions = zeros(size(x));
Xi1 = previousKnotVector( 1 + problem.p );
Xi2 = previousKnotVector( 2 + problem.p );
numericalSolutions(x>=Xi1 & x<=Xi2) = localProjection(x(x>=Xi1 & x<=Xi2), 1, coefficients, problem);

for e=2:numel(previousKnotVector)-problem.p-1
    Xi1 = previousKnotVector( e + problem.p );
    Xi2 = previousKnotVector( e + 1 + problem.p );
    
    numericalSolutions(x>Xi1 & x<=Xi2) = localProjection(x(x>Xi1 & x<=Xi2), e, coefficients, problem);
end
end

function r = localProjection(x, element, coefficients, problem)
% r = LOCALPROJECTION(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
    
    r = zeros(size(x));
    
    if isempty(r) == 0
        [N, ~] = BsplinesShapeFunctionsAndDerivatives(x, problem.p, problem.knotVector);
        r = r + N(problem.LM(element,:)) * coefficients(problem.LM(element,:));
    end

    
end