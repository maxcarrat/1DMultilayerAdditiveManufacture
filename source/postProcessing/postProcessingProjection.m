function [ projectedCoefficients ] = postProcessingProjection(x, t, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTION project the previous solution onto the post-processing mesh.
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives


projectedCoefficients=zeros(size(x));
X1 = problem.coords(1);
X2 = problem.coords(2);

projectedCoefficients(x>=X1 & x<=X2) = postProcessingProjectionElement(1, x(x>=X1 & x<=X2), t, problem,... % first element
    solutionCoefficients, modes, derivative);

% loop over elements
for e=2:problem.N
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    projectedCoefficients(x>X1 & x<=X2) = postProcessingProjectionElement(e, x(x>X1 & x<=X2), t, problem,...
        solutionCoefficients, modes, derivative);
    
end

end

function [ projectedCoefficients ] = postProcessingProjectionElement(e, x, t, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTIONELEMENT project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);

X1 = problem.coords(e);
X2 = problem.coords(e+1);

if e > problem.N - problem.XN  % element is active
    
    elementEnrichedIndex = e - (problem.N - problem.XN);

    if elementEnrichedIndex == 1
        indexLocalEnrichedNodes = 2; %rhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    refinedNodes = 2^problem.refinementDepth+1;
    integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    subDomainShapeFunctionCoefficients = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    localCoordinates = mapGlobalToLocal( x, X1, X2);
    projectedCoefficients=zeros(size(localCoordinates));
    
    projectedCoefficients(localCoordinates>=Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, 1, elementEnrichedIndex,...
        localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
        solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
    
    for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
        
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
        
        projectedCoefficients(localCoordinates>Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, integrationSubDomainIndex, elementEnrichedIndex,...
            localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
            solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
          
  
    end
    
    projectedCoefficients = projectedCoefficients .* (2/(X2-X1)) ^ derivative;

    
    
else     % Element not enriched
    
    projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
    
    N = zeros(length(x), 2);
    B = zeros(length(x), 2);
    
    localCoords = mapGlobalToLocal( x, X1, X2);
    
    for k=1:length(x)
        [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    end
    
    if derivative == 0
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
    else
        for i=1:length(x)
            projectionOperator(i,1:size(B,2)) = B(i,:);
        end
    end
    
    projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:)) .* (2/(X2-X1)) ^ derivative;
    
    if derivative == 1
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
        temperature = projectionOperator * solutionCoefficients(problem.LM(e,:));
        for i=1:numel(projectedCoefficients)
            projectedCoefficients(i) = projectedCoefficients(i) * ...
            problem.k(x, t, temperature(i));
        end

    end
end


end


function [ projectedCoefficients ] = postProcessingProjectionSubElements(globalCoords, subDomainIndex, e, x, t, problem, integrationDomain, solutionCoefficients,...
    modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes)
% POSTPROCESSINGPROJECTIONSUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, 2 + modes * length(indexLocalEnrichedNodes) );

localCoords = x;
refinedNodes = 2^problem.refinementDepth + 1;
N = zeros(length(x), 2);
B = zeros(length(x), 2);

F = zeros(length(x), modes *length(indexLocalEnrichedNodes));
G = zeros(length(x), modes *length(indexLocalEnrichedNodes));

PODCoefficients = problem.reductionOperator(...
    (e-1)*(floor(refinedNodes/problem.XN))+1:(e-1)*(floor(refinedNodes/problem.XN)) +...
     ceil(refinedNodes/problem.XN),:);

for k=1:length(x)
    
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    
    [F(k,:), G(k,:)] = PODModesAndDerivativesGaussIntegration(localCoords(k), modes, PODCoefficients,...
        subDomainShapeFunctionCoefficients, subDomainIndex, indexLocalEnrichedNodes);
end

if derivative == 0
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
else
    for i=1:length(localCoords)
        projectionOperator(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
    end
end

if length(indexLocalEnrichedNodes) == 1 %lhs
    coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
        solutionCoefficients(problem.N  + (e - 1) * length(indexLocalEnrichedNodes) * modes + 2 : problem.N...
        + (e) * length(indexLocalEnrichedNodes) * modes + 1 )];
else
    coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
        solutionCoefficients(problem.N  + (e - 2) * modes + 2 : problem.N...
        + (e - 2) * modes + 2 + length(indexLocalEnrichedNodes)  * modes - 1 )];
end

projectedCoefficients = projectionOperator * coefficients ;

if derivative == 1
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
    
    temperature = projectionOperator * coefficients ;
    for i=1:numel(x)
        projectedCoefficients(i) = projectedCoefficients(i) .* ...
            problem.k(globalCoords(i), t, temperature(i));
    end
end

end