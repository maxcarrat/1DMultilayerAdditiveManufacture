function [ projectedTemperature ] = projectOntoXIGAMesh( problem, previousProblem, temperatureCoefficients,...
    modes, knotVector,  previousKnotVector, numberOfEnrichedControlPoints, initialTemperature )
%PROJECTONTOXIGAMESH

projectedTemperature = zeros( size(knotVector, 2) - problem.p - 1 + numberOfEnrichedControlPoints*modes, 1);

for j=1:numel(problem.coords)
    x = problem.coords(j);
    globalProjectedValue = globalProjection(x, previousKnotVector, temperatureCoefficients, previousProblem);
    
    if globalProjectedValue ~= 0.0
        projectedTemperature(j) = globalProjectedValue;
    else
        projectedTemperature(j) = 0.0; %initialTemperature;
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

X1 = mapParametricToGlobal(Xi1, problem);
X2 = mapParametricToGlobal(Xi2, problem);

parentCoord = mapGlobalToLocal(x, X1, X2);
parametricCoord = mapParentToLocal(parentCoord, Xi1, Xi2);

numericalSolutions(x>=X1 & x<=X2) = localProjection(parametricCoord(x>=X1 & x<=X2), 1, coefficients, problem);

for e=2:numel(previousKnotVector)-problem.p-1
    Xi1 = previousKnotVector( e + problem.p );
    Xi2 = previousKnotVector( e + 1 + problem.p );
    
    X1 = mapParametricToGlobal(Xi1, problem);
    X2 = mapParametricToGlobal(Xi2, problem);
    
    parentCoord = mapGlobalToLocal(x, X1, X2);
    parametricCoord = mapParentToLocal(parentCoord, Xi1, Xi2);

    numericalSolutions(x>X1 & x<=X2) = localProjection(parametricCoord(x>X1 & x<=X2), e, coefficients, problem);
end

end

function r = localProjection(parametricCoord, element, coefficients, problem)
% r = LOCALPROJECTION(parametricCoord, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   parametricCoord = points where the element numerical solution has to be evaluated
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations

r = zeros(size(parametricCoord));

if isempty(r) == 0
    [N, ~] = BsplinesShapeFunctionsAndDerivatives(parametricCoord, problem.p, problem.knotVector);
    r = r + N(problem.LM(element,:)) * coefficients(problem.LM(element,:));
end

end