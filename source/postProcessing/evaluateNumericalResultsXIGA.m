function [ numericalSolutions ] = evaluateNumericalResultsXIGA( x, t, modes, problem,...
    coefficients, layer, numberOfLayers, derivative )
%EVALUATEUMERICALRESULTSXIGA evaluates the numerical solution of X-PODIGA
%Input:
%x = coordinates to post process
%t = actual time
%modes = number of POD modes
%problem = struct that defines the boundary value problem
%coefficients = coefficients of the basis function obtained by solving the 
%mass matrix-load vector system of equations
%layer = layer index
%numberOfLayers = total number of layer in the process
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

numericalSolutions=zeros(size(x));

%knot span lft and right end
Xi1 = problem.knotVector( 1 + problem.p );
Xi2 = problem.knotVector( 2 + problem.p );

%local coordinates in post-processing
localCoordinates = linspace(0, 1, length(x) / numberOfLayers * layer );

%evaluate numerical solution at the points in the 1st element
numericalSolutions(localCoordinates>=Xi1 & localCoordinates<=Xi2) = element_num_sol...
    ( localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), t, modes, problem, 1, coefficients, derivative);

%loop over remaining elements
for e=2:size(problem.LM, 1)
    Xi1 = problem.knotVector( e + problem.p );
    Xi2 = problem.knotVector( e + 1 + problem.p );
    
    numericalSolutions(localCoordinates>Xi1 & localCoordinates<=Xi2) = element_num_sol...
        ( localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), t, modes, problem, e, coefficients, derivative);
    
end

end

function r = element_num_sol(x, t, modes, problem, e, coefficients, derivative)
%ELEMENT_NUM_SOL evaluates the numerical solution associated with a single specific element
%Input:
%x = points where the element numerical solution has to be evaluated
%t = actual time
%modes = number of POD modes
%problem = struct that defines the boundary value problem
%e = index of the element where to evaluate the element numerical solution
%coefficients = coefficients of the basis function obtained by solving the 
%mass matrix-load vector system of equations
%derivative = index of the deriative of the element numerical derivative to
%be evaluated

numberOfProjectionPoints = length(x);
m = length(problem.knotVector);

%% Initialize variables

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(x), m - 1 - problem.p);
B = zeros(length(x), m - 1 - problem.p);
JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

%check if element is enriched
if e > problem.N - problem.XN  % element is enriched
    
    elementEnrichedIndex = e - (problem.N - problem.XN);
    
    if elementEnrichedIndex == 1 %if last element
        indexLocalEnrichedNodes = problem.IGAdof;
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    refinedNodes = length(problem.reductionOperator)-1;    
    integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    
    %left and right end of the integration sub-domain
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    %left and right end of the knot vector
    Xp1 = problem.knotVector(e + problem.p);
    Xp2 = problem.knotVector(e + problem.p + 1);
    
    %map integration domain ends from parent to sub-domain space
    xLocal1 = mapParentToLocal(Xi1, Xp1, Xp2);
    xLocal2 = mapParentToLocal(Xi2, Xp1, Xp2);
    
    localCoordinates = x;
    projectedCoefficients=zeros(size(localCoordinates));
    
    %post-process on 1st sub-domain
    projectedCoefficients(localCoordinates>=xLocal1 & localCoordinates<=xLocal2) =...
        postProcessingProjectionSubElements(1, elementEnrichedIndex,...
        localCoordinates(localCoordinates>=xLocal1 & localCoordinates<=xLocal2), t, problem, e,...
        coefficients, modes, derivative, integrationDomain, indexLocalEnrichedNodes);
 
    %loop over remaining sub-domains
    for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
        
        %left and right end of the integration sub-domain
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
 
        %map integration domain ends from parent to sub-domain space
        xLocal1 = mapParentToLocal(Xi1, Xp1, Xp2);
        xLocal2 = mapParentToLocal(Xi2, Xp1, Xp2);
        
        projectedCoefficients(localCoordinates>xLocal1 & localCoordinates<=xLocal2) = postProcessingProjectionSubElements(integrationSubDomainIndex, elementEnrichedIndex,...
            localCoordinates(localCoordinates>xLocal1 & localCoordinates<=xLocal2), t, problem, e,...
            coefficients, modes, derivative, integrationDomain, indexLocalEnrichedNodes);
    end
    
    r = projectedCoefficients;
    
else % element is not enriched
    
    %evaluate temperature
    if derivative == 0
        %loop over post-processing points in the eth element
        for i=1:length(x)
            [N(i,:), ~] = BsplinesShapeFunctionsAndDerivatives(x(i), problem.p, problem.knotVector);
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
        r = projectionOperator(:,problem.LM(e,:)) * coefficients(problem.LM(e,:)) ;
        
    %evaluate heat fluxes
    else
        
        %loop over post-processing points in the eth element
        for i=1:length(x)
            [N(i,:), B(i,:)] = BsplinesShapeFunctionsAndDerivatives(x(i), problem.p, problem.knotVector);
            projectionOperator(i,1:size(N,2)) = N(i,:);
            projectionOperator_der(i,1:size(B,2)) = B(i,:);
            
            JacobianX_Xi(i) = B(i,problem.LM(e, :)) *...
                problem.coords(problem.LM(e, :))';
            inverseJacobianX_Xi(i) = 1 / JacobianX_Xi(i);
        end
        
        r = (projectionOperator_der(:,problem.LM(e,:)) * coefficients(problem.LM(e,:))) .* inverseJacobianX_Xi .*  ...
            problem.k(mapParametricToGlobal(x, problem), t, projectionOperator(:,problem.LM(e,:)) * coefficients(problem.LM(e,:))) ;
    end
end


end



function [ projectedCoefficients ] = postProcessingProjectionSubElements(subDomainIndex,...
    elementEnrichedIndex, x, t, problem, element, solutionCoefficients,...
    modes, derivative, integrationDomain, indexLocalEnrichedNodes)
%POSTPROCESSINGPROJECTIONSUBELEMENTS project the previous solution onto the element.
%Input:
%subDomainIndex = index of the actual sub-domain
%elementEnrichedIndex = index of the element enriched
%x = post-processing mesh in local coordinates of the integration domain
%t = actual time
%problem = struct of teh poisson initial boundary value problem
%element = element index
%solutionCoefficients = coefficients of the solutions
%modes = number of enrichment modes
%derivative = order of derivatives
%integrationDomain = integration domain of the enriched element in local
%indexOfLocalEnrichedNodes = index of the enriched degrees of freedom
%Output:
%projectedCoefficients = value of the solution at the post-processing
%points

%% Initialize variables
numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, problem.p+1 + modes * length(indexLocalEnrichedNodes) );
projectionOperator_der = zeros(numberOfProjectionPoints, problem.p+1 + modes * length(indexLocalEnrichedNodes) );

localCoords = x;
refinedNodes = length(problem.reductionOperator)-1;

N = zeros(length(x), length(problem.knotVector)-problem.p-1);
B = zeros(length(x), length(problem.knotVector)-problem.p-1);

JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

F = zeros(length(x), modes * length(indexLocalEnrichedNodes));
G = zeros(length(x), modes * length(indexLocalEnrichedNodes));

PODCoefficients = problem.reductionOperator(...
    (elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN))+1:(elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN)) +...
    ceil(refinedNodes/problem.XN),:);

%loop over post-processing points in the sub-domain
for k=1:length(x)
    [N(k,:), B(k,:)] = BsplinesShapeFunctionsAndDerivatives(localCoords(k), problem.p, problem.knotVector);
    [F(k,:), G(k,:)] = PODModesAndDerivativesIGA( problem, localCoords(k), modes, PODCoefficients,...
        integrationDomain, subDomainIndex, indexLocalEnrichedNodes, element, problem.knotVector );
end

if length(indexLocalEnrichedNodes) == 1 %lhs
    coefficients = [solutionCoefficients((problem.N - problem.XN) + elementEnrichedIndex : (problem.N - problem.XN) + elementEnrichedIndex + 1); ...
        solutionCoefficients(problem.N  + (elementEnrichedIndex - 1) * length(indexLocalEnrichedNodes) * modes + 2 : problem.N...
        + (elementEnrichedIndex) * length(indexLocalEnrichedNodes) * modes + problem.p )];
else
    coefficients = [solutionCoefficients((problem.N - problem.XN) + elementEnrichedIndex : (problem.N - problem.XN) + elementEnrichedIndex + 1); ...
        solutionCoefficients(problem.N  + (elementEnrichedIndex - 2) * modes + 2 : problem.N...
        + (elementEnrichedIndex - 2) * modes + length(indexLocalEnrichedNodes)  * modes + problem.p )];
end

%evaluate temperature values in te sub-domain
if derivative == 0   
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
    
    projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * coefficients ;

%evaluate heat fluxes in the sub-domain
else    
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
        projectionOperator_der(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
    end
    
    for i=1:length(localCoords)
        JacobianX_Xi(i) = B(i,problem.LM(element, :)) *...
            problem.coords(problem.LM(element, :))';
        inverseJacobianX_Xi(i) = 1 / norm(JacobianX_Xi(i));
    end
    
    projectedCoefficients_der = projectionOperator_der(:, end-problem.p-modes:end) * coefficients ;
    projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * coefficients ;
    
    globalCoords = mapParametricToGlobal(localCoords, problem);
    
    for i=1:numel(projectedCoefficients_der)
        projectedCoefficients_der(i) = projectedCoefficients_der(i) * inverseJacobianX_Xi(i) * ...
            problem.k(globalCoords(i), t, projectedCoefficients(i));
    end
    
    projectedCoefficients = projectedCoefficients_der;
    
end

end