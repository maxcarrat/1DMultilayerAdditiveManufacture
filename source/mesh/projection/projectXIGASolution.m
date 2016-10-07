function [ projectedTemperature ] = projectXIGASolution( problem, temperatureCoefficients, knotVector,...
    previousKnotVector, numberOfEnrichedControlPoints, initialTemperature )
%PROJECTXIGASOLUTION This function 1st project the coefficients of CPs  of
%the previous mesh onto a set x, of equally distributed evaluation points,
%in the new mesh. Secondly, the evaluation points values are used to
%calculate the coefficients of the CPs of the new mesh.

% initialize vector of the coefficients at the eveluated points
evaluationPointsGlobal = linspace(problem.coords(1), problem.coords(end), length(problem.coords)*10);

evaluatedProjectedTemperature = zeros( length(evaluationPointsGlobal), 1);

for j=1:numel(evaluationPointsGlobal)
%     x = problem.coords(j);

    % evaluate the value at x, where x is a set of points along the new
    % domain mesh. 
    globalProjectedValue = globalProjection(evaluationPointsGlobal(j), previousKnotVector, temperatureCoefficients, problem);
    
    if globalProjectedValue ~= 0.0
        evaluatedProjectedTemperature(j) = globalProjectedValue;
    else
        evaluatedProjectedTemperature(j) = initialTemperature;
    end
end

% evaluate teh coefficients of the control points on the new mesh
N = assembleBSplineMatrix(problem, evaluationPointsGlobal, knotVector);
projectedTemperature = N \ evaluatedProjectedTemperature;

% set eXtended value to zero
projectedTeperatureExtended = zeros(problem.modes*numberOfEnrichedControlPoints, 1);

projectedTemperature = [projectedTemperature; projectedTeperatureExtended];

end

function [ N ] = assembleBSplineMatrix(problem, x, knotVector)
%N = ASSEMBLEBSPLINEMATRIX(problem, x)
%input:
% x = vector of evaluated global points coordinates
% problem = eXtended Poisson problem struct
%output:
% N = matrix of BSpline at x position

N = [];
evaluationPointsParametric = linspace(0.0, 1.0, length(x));

for j=1:length(x)
    [N_j, ~] = BsplinesShapeFunctionsAndDerivatives(evaluationPointsParametric(j), problem.p, knotVector);
    N = [N; N_j];
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

for e=2:numel(previousKnotVector)- 2 * problem.p - 1
    
    if e > problem.N - problem.XN % element of the previous mesh is enriched
        
        elementEnrichedIndex = e - (problem.N - problem.XN);
        
        if elementEnrichedIndex == 1
            indexLocalEnrichedNodes = problem.IGAdof;
        else
            indexLocalEnrichedNodes = [1, 2];
        end
        
        % On active elements use the refined domain as integration domain
        refinedNodes = 2^problem.refinementDepth+problem.p;
        
        Xi1 = previousKnotVector( e + problem.p-1 );
        Xi2 = previousKnotVector( e + 1 + problem.p-1 );
        
        X1 = mapParametricToGlobal(Xi1, problem);
        X2 = mapParametricToGlobal(Xi2, problem);
        
        parentCoord = mapGlobalToLocal(x, X1, X2);
        parametricCoord = mapParentToLocal(parentCoord, Xi1, Xi2);
        
        integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
        Xiparent1 = integrationDomain(1);
        Xiparent2 = integrationDomain(2);
        
        numericalSolutions(parentCoord>=Xiparent1 & parentCoord<=Xiparent2) = localSubElementsProjection(1, elementEnrichedIndex, ...
            parametricCoord(parentCoord>=Xiparent1 & parentCoord<=Xiparent2), problem, e, coefficients, integrationDomain,...
            indexLocalEnrichedNodes, previousKnotVector);
        
        for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
            
            Xiparent1 = integrationDomain(integrationSubDomainIndex);
            Xiparent2 = integrationDomain(integrationSubDomainIndex + 1);
        
            numericalSolutions(parentCoord>Xiparent1 & parentCoord<=Xiparent2) = localSubElementsProjection(integrationSubDomainIndex, elementEnrichedIndex,...
                parametricCoord(parentCoord>Xiparent1 & parentCoord<=Xiparent2), problem, e, coefficients,...
                integrationDomain, indexLocalEnrichedNodes, previousKnotVector);
        end

    else % element is not enriched
        
        Xi1 = previousKnotVector( e + problem.p );
        Xi2 = previousKnotVector( e + 1 + problem.p );
        
        X1 = mapParametricToGlobal(Xi1, problem);
        X2 = mapParametricToGlobal(Xi2, problem);
        
        parentCoord = mapGlobalToLocal(x, X1, X2);
        parametricCoord = mapParentToLocal(parentCoord, Xi1, Xi2);
        
        numericalSolutions(x>X1 & x<=X2) = localProjection(parametricCoord(x>X1 & x<=X2), e, coefficients, problem);
        
    end
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

function [ projectedCoefficients ] = localSubElementsProjection(subDomainIndex,...
    elementEnrichedIndex, x, problem, element, solutionCoefficients,...
    integrationDomain, indexLocalEnrichedNodes, previousKnotVector)
% LOCALSUBELEMENTSPROJECTION project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   problem.modes = number of enrichment problem.modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, length(problem.knotVector)-problem.p-2 + problem.modes * length(indexLocalEnrichedNodes) );

localCoords = x;
refinedNodes = 2^problem.refinementDepth + problem.p;

N = zeros(length(x), length(problem.knotVector)-problem.p-1);

F = zeros(length(x), problem.modes * length(indexLocalEnrichedNodes));

if isempty(localCoords) == 0
    
    PODCoefficients = problem.reductionOperator(...
        (elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN))+1:(elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN)) +...
        ceil(refinedNodes/problem.XN),:);
    
    for k=1:length(x)
        
        [N(k,:), ~] = BsplinesShapeFunctionsAndDerivatives(localCoords(k), problem.p, previousKnotVector);
        [F(k,:), ~] = PODModesAndDerivativesIGA( problem, localCoords(k), problem.modes, PODCoefficients,...
            integrationDomain, subDomainIndex, indexLocalEnrichedNodes, element, previousKnotVector );
        
    end
    
    if length(indexLocalEnrichedNodes) == 1 %lhs
        coefficients = [solutionCoefficients((problem.N - problem.XN) : (problem.N - problem.XN) + elementEnrichedIndex + 1); ...
            solutionCoefficients(problem.N + problem.p + 1 + (elementEnrichedIndex - 1) * length(indexLocalEnrichedNodes) * problem.modes : end )];
    else
        coefficients = [solutionCoefficients((problem.N - problem.XN) + elementEnrichedIndex : (problem.N - problem.XN) + elementEnrichedIndex + 1); ...
            solutionCoefficients(problem.N  + (elementEnrichedIndex - 2) * problem.modes + 2 : problem.N...
            + (elementEnrichedIndex - 2) * problem.modes + length(indexLocalEnrichedNodes)  * problem.modes + problem.p )];
    end
    
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)-1) = [N(i,1:end-1), F(i,:)];
    end
    
    projectedCoefficients = projectionOperator(:, end-problem.p-problem.modes:end) * coefficients ;
    
else
    projectedCoefficients = [];
end

end
