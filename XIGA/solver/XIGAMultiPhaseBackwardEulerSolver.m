function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]=...
    XIGAMultiPhaseBackwardEulerSolver( p, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance,...
    maxIterations, numberOfRefinedElementsToBeKept,...
    TotalNumberOfControlPoints, refinementDepth, numberOfEnrichedControlPoints,...
    numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes,...
    integrationOrder, integrationModesOrder, layerThickness)
% IGAMultiPhaseBackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   CPs = control points
%   rhs = external heat source
%   leftDirichletBoundaryConditionValue
%   rightDirichletBoundaryConditionValue
%   k = conductivity of the material
%   heatCapacity = heat capacity of the material
%   timeVector = vector of time steps for Backward Euler implicit scheme
%   refinementDepth = depth of the refined mesh
%   numberOfTrainingTimeSteps = number of time steps in the training phase

timeSteps=size(timeVector,2);

timeStepSize=max(timeVector)/( timeSteps );
temperaturePostProcessing = zeros(length(postProcessingCoords), timeSteps);

convergence = 0.0;
computationalTime = [];

heatFluxes = zeros(length(postProcessingCoords), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for layer = 1:numberOfTrainingLayers
    
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %% Refinement
    % In order to get POD modes corresponding to the last layer the
    % interface between the last ( refined ) layer and the other has to be
    % C0 continuous, i.e. we have to insert an additional knot at the layer
    % interface.

    newKnots = linspace( knotVector(end-(p+1)),  knotVector(end-p), 2^refinementDepth + 1);
    
    if layer == 1 || p == 1
        [ CPs, refinedKnotVector] = refineKnotVector( p, CPs, knotVector, newKnots(2:end-1));
    else
        [ CPs, refinedKnotVector] = refineKnotVector( p, CPs, knotVector, newKnots(1:end-1));
    end
    
    %Project old solution onto the new mesh
    if layer > 1
        refinedTemperatureSolutions = L2projectionIGA(poissonProblem,...
            refinedTemperatureSolutions, CPs, previousKnotVector, initialTemperature);
    else
        refinedTemperatureSolutions = zeros(length(CPs), 1);
        localRefinedTemperatureSolutions = zeros(length(CPs), timeSteps);
    end
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    % foor loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonProblem = poissonProblemTransientIGA(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            refinedKnotVector, p, refinementDepth);
        
        %Update and merge temperature into global domain
        [refinedTemperatureSolutions, convergenceFlag] = solveMultiPhaseProblemIGA( poissonProblem,...
            currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, TotalNumberOfControlPoints);
        
        %Cache the local solution for POD
        localRefinedTemperatureSolutions(:, t+1) = getIGASolution(refinedTemperatureSolutions,...
            layer, poissonProblem, numberOfRefinedElementsToBeKept);
        
        %Cache knot vector
        previousKnotVector = poissonProblem.knotVector;
        
        %Update convergence flag and register time to solve the time step
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 1);
        
        %Register number of Dofs
        DOFs = numel(refinedTemperatureSolutions);
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

%% Generate the reduced basis
% POD on the local/layer solution snapshots, omitt the first zero-solution
% vector
[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,5*numberOfLayersTimeSteps+2:numberOfTrainingLayers*numberOfLayersTimeSteps+1), numberOfPODModes);

%% Enriched mesh using RB
% for loop ROM phase layers
for layer = (numberOfTrainingLayers+1):numberOfLayers
    
%     p = 2;
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %Enrich only the bar tip CP
    activeNumberOfCPs = numberOfEnrichedControlPoints;
    
    %Generate the XIGA Poisson problem at layer for projection
    poissonXProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
        leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
        neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
        knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
    
    %Project global solution onto the enriched modal space
    temperatureSolutions = projectOntoXIGAMesh(poissonXProblem, poissonProblem, refinedTemperatureSolutions,...
        modes, knotVector, previousKnotVector, numberOfEnrichedControlPoints, initialTemperature);
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    % for loop time integration @layer    
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate Local problem
        poissonXProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        [temperatureSolutions, convergenceFlag] = solveXIGAMultiPhase( poissonXProblem, currentTime,...
            timeStepSize,integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );

        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
        poissonXProblem, temperatureSolutions,...
            layer, numberOfLayers, 0);
        heatFluxes(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
            poissonXProblem,temperatureSolutions,...
            layer, numberOfLayers, 1);
        
%         temperaturePostProcessing(:, t+1) = postProcessingProjection(postProcessingCoords,...
%             currentTime, poissonXProblem, temperatureSolutions, modes, 0.0);
%         heatFluxes(:, t+1) = postProcessingProjection(postProcessingCoords, currentTime,...
%             poissonXProblem,temperatureSolutions, modes, 1);
        
        DOFs = numel(temperatureSolutions);
        
    end
    
    poissonProblem = poissonXProblem;
    previousKnotVector = knotVector;
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IGA POST-PROCESSING

function [ numericalSolutions ] = evaluateNumericalResultsXIGA( x, t, modes, problem,...
    coefficients, layer, numberOfLayers, derivative )
% numericalSolutions = evaluateNumerialResultsIGA(postProcessingCoords, problem, coefficients, derivative) evaluates the numerical solution
% x = coordinates to post process
% problem = struct that defines the boundary value problem
% coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
% derivative = index of the deriative of the element numerical derivative to be evaluated

numericalSolutions=zeros(size(x));
Xi1 = problem.knotVector( 1 + problem.p );
Xi2 = problem.knotVector( 2 + problem.p );

x = linspace(0, 1, length(x) / numberOfLayers * layer );
numericalSolutions(x>=Xi1 & x<=Xi2) = element_num_sol( x(x>=Xi1 & x<=Xi2), t, modes, problem, 1, coefficients, derivative);

for e=2:size(problem.LM, 1)
    Xi1 = problem.knotVector( e + problem.p );
    Xi2 = problem.knotVector( e + 1 + problem.p );
      
numericalSolutions(x>Xi1 & x<=Xi2) = element_num_sol( x(x>Xi1 & x<=Xi2), t, modes, problem, e, coefficients, derivative);

end

end

function r = element_num_sol(x, t, modes, problem, e, coefficients, derivative)
% r = ELEMENT_NUM_SOL(x, coords, p, problem, element, coefficients, derivative) evaluates the numerical solution associated with a single specific element
%   x = points where the element numerical solution has to be evaluated
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   element = index of the element where to evaluate the element numerical solution
%   coefficients = coefficients of the basis function obtained by solving the mass matrix-load vector system of equations
%   derivative = index of the deriative of the element numerical derivative to be evaluated

numberOfProjectionPoints = length(x);
m = length(problem.knotVector);

projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
projectionOperator_der = zeros(numberOfProjectionPoints, size(problem.LM, 2));

N = zeros(length(x), m - 1 - problem.p);
B = zeros(length(x), m - 1 - problem.p);
JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

for k=1:length(x)
    [N(k,:), B(k,:)] = BsplinesShapeFunctionsAndDerivatives(x(k), problem.p, problem.knotVector);
end

if e > problem.N - problem.XN  % element is enriched
 
    elementEnrichedIndex = e - (problem.N - problem.XN);
    
    if elementEnrichedIndex == 1
        indexLocalEnrichedNodes = problem.IGAdof;
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    refinedNodes = 2^problem.refinementDepth+1;
    integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    subDomainShapeFunctionCoefficients = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    Xp1 = problem.knotVector(e + problem.p);
    Xp2 = problem.knotVector(e + problem.p + 1);
    
    xLocal1 = mapParentToLocal(Xi1, Xp1, Xp2);
    xLocal2 = mapParentToLocal(Xi2, Xp1, Xp2);
    
    localCoordinates = x;
    projectedCoefficients=zeros(size(localCoordinates));
    
    projectedCoefficients(localCoordinates>=xLocal1 & localCoordinates<=xLocal2) = postProcessingProjectionSubElements(1, elementEnrichedIndex,...
        localCoordinates(localCoordinates>=xLocal1 & localCoordinates<=xLocal2), t, problem, e,...
        coefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);    
    
    for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
        
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
        
        xLocal1 = mapParentToLocal(Xi1, Xp1, Xp2);
        xLocal2 = mapParentToLocal(Xi2, Xp1, Xp2);
        
        projectedCoefficients(localCoordinates>xLocal1 & localCoordinates<=xLocal2) = postProcessingProjectionSubElements(integrationSubDomainIndex, elementEnrichedIndex,...
            localCoordinates(localCoordinates>xLocal1 & localCoordinates<=xLocal2), t, problem, e,...
            coefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
    end    
    
    r = projectedCoefficients;% .* (2/(Xp2-Xp1)) ^ derivative;


else % element is not enriched
    
    if derivative == 0
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
        r = projectionOperator(:,problem.LM(e,:)) * coefficients(problem.LM(e,:)) ;
    else
        for i=1:length(x)
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
projectionOperator_der = zeros(numberOfProjectionPoints, 2 + modes * length(indexLocalEnrichedNodes) );
    
localCoords = x;
refinedNodes = 2^problem.refinementDepth + problem.p;
% refinedNodes = 2^problem.refinementDepth + 1;

N = zeros(length(x), length(problem.knotVector)-problem.p-1);
B = zeros(length(x), length(problem.knotVector)-problem.p-1);

Xp1 = problem.knotVector( element + problem.p);
Xp2 = problem.knotVector( element + problem.p + 1);

JacobianX_Xi = zeros(length(x), 1);
inverseJacobianX_Xi = zeros(length(x), 1);

F = zeros(length(x), modes * length(indexLocalEnrichedNodes));
G = zeros(length(x), modes * length(indexLocalEnrichedNodes));

PODCoefficients = problem.reductionOperator(...
    (elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN))+1:(elementEnrichedIndex-1)*(floor(refinedNodes/problem.XN)) +...
    ceil(refinedNodes/problem.XN),:);

for k=1:length(x)
    
    [N(k,:), B(k,:)] = BsplinesShapeFunctionsAndDerivatives(localCoords(k), problem.p, problem.knotVector);

    [F(k,:), G(k,:)] = PODModesAndDerivativesIGA( problem, localCoords(k), modes, PODCoefficients,...
        subDomainShapeFunctionCoefficients, subDomainIndex, indexLocalEnrichedNodes, element );
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

if derivative == 0
    
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end

    projectedCoefficients = projectionOperator(:, end-problem.p-modes:end) * coefficients ;
    
else
    
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
        projectionOperator_der(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
    end
    
    for i=1:length(localCoords)
        JacobianX_Xi(i) = B(i,problem.LM(element, :)) *...
            problem.coords(problem.LM(element, :))';
        inverseJacobianX_Xi(i) = 1 / JacobianX_Xi(i);
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% LINEAR POST-PROCESSING
% function [ projectedCoefficients ] = postProcessingProjection(x, t, problem, solutionCoefficients,...
%     modes, derivative)
% % POSTPROCESSINGPROJECTION project the previous solution onto the post-processing mesh.
% %   x = post-processing mesh
% %   problem
% %   solutionCoefficients = temeprature distribution of the previous mesh
% %   modes = number of enrichment modes
% %   derivative = order of derivatives
% 
% 
% projectedCoefficients=zeros(size(x));
% X1 = mapParametricToGlobal(problem.knotVector(problem.p+1), problem);
% X2 = mapParametricToGlobal(problem.knotVector(problem.p+2), problem);
% 
% projectedCoefficients(x>=X1 & x<=X2) = postProcessingProjectionElement(1, x(x>=X1 & x<=X2), t, problem,... % first element
%     solutionCoefficients, modes, derivative);
% 
% % loop over elements
% for e=2:problem.N
%     
%     X1 = mapParametricToGlobal(problem.knotVector(problem.p+e), problem);
%     X2 = mapParametricToGlobal(problem.knotVector(problem.p+e+1), problem);
%     
%     projectedCoefficients(x>X1 & x<=X2) = postProcessingProjectionElement(e, x(x>X1 & x<=X2), t, problem,...
%         solutionCoefficients, modes, derivative);
%     
% end
% 
% end
% 
% function [ projectedCoefficients ] = postProcessingProjectionElement(e, x, t, problem, solutionCoefficients,...
%     modes, derivative)
% % POSTPROCESSINGPROJECTIONELEMENT project the previous solution onto the element.
% %   e = element index
% %   x = post-processing mesh
% %   problem
% %   solutionCoefficients = temeprature distribution of the previous mesh
% %   modes = number of enrichment modes
% %   derivative = order of derivatives
% 
% numberOfProjectionPoints = length(x);
% 
% X1 = mapParametricToGlobal(problem.knotVector(problem.p+e), problem);
% X2 = mapParametricToGlobal(problem.knotVector(problem.p+e+1), problem);
% 
% if e > problem.N - problem.XN  % element is active
%     
%     elementEnrichedIndex = e - (problem.N - problem.XN);
% 
%     if elementEnrichedIndex == 1
%         indexLocalEnrichedNodes = 2; %rhs node
%     else
%         indexLocalEnrichedNodes = [1, 2];
%     end
%     
%     % On active elements use the refined domain as integration domain
%     refinedNodes = 2^problem.refinementDepth+1;
%     integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
%     subDomainShapeFunctionCoefficients = linspace(-1, +1, ceil(refinedNodes/problem.XN));
%     
%     Xi1 = integrationDomain(1);
%     Xi2 = integrationDomain(2);
%     
%     localCoordinates = mapGlobalToLocal( x, X1, X2);
%     projectedCoefficients=zeros(size(localCoordinates));
%     
%     projectedCoefficients(localCoordinates>=Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, 1, elementEnrichedIndex,...
%         localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
%         solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
%     
%     for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
%         
%         Xi1 = integrationDomain(integrationSubDomainIndex);
%         Xi2 = integrationDomain(integrationSubDomainIndex+1);
%         
%         projectedCoefficients(localCoordinates>Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, integrationSubDomainIndex, elementEnrichedIndex,...
%             localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
%             solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
% 
%     end
%     
%     projectedCoefficients = projectedCoefficients .* (2/(X2-X1)) ^ derivative;
%    
% else     % Element not enriched
%     
%     projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
%     
%     N = zeros(length(x), 2);
%     B = zeros(length(x), 2);
%     
%     localCoords = mapGlobalToLocal( x, X1, X2);
%     
%     for k=1:length(x)
%         [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
%     end
%     
%     if derivative == 0
%         for i=1:length(x)
%             projectionOperator(i,1:size(N,2)) = N(i,:);
%         end
%     else
%         for i=1:length(x)
%             projectionOperator(i,1:size(B,2)) = B(i,:);
%         end
%     end
%     
%     projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:)) .* (2/(X2-X1)) ^ derivative;
%     
%     if derivative == 1
%         for i=1:length(x)
%             projectionOperator(i,1:size(N,2)) = N(i,:);
%         end
%         temperature = projectionOperator * solutionCoefficients(problem.LM(e,:));
%         for i=1:numel(projectedCoefficients)
%             projectedCoefficients(i) = projectedCoefficients(i) * ...
%             problem.k(x, t, temperature(i));
%         end
% 
%     end
% end
% 
% 
% end
% 
% 
% function [ projectedCoefficients ] = postProcessingProjectionSubElements(globalCoords, subDomainIndex, e, x, t, problem, integrationDomain, solutionCoefficients,...
%     modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes)
% % POSTPROCESSINGPROJECTIONSUBELEMENTS project the previous solution onto the element.
% %   e = element index
% %   x = post-processing mesh in local coordinates of the integration domain
% %   problem
% %   integrationDomain = integration domain of the enriched element in local
% %   coords
% %   solutionCoefficients = temeprature distribution of the previous mesh
% %   modes = number of enrichment modes
% %   derivative = order of derivatives
% 
% numberOfProjectionPoints = length(x);
% projectionOperator = zeros(numberOfProjectionPoints, 2 + modes * length(indexLocalEnrichedNodes) );
% 
% localCoords = x;
% refinedNodes = 2^problem.refinementDepth + 1;
% N = zeros(length(x), 2);
% B = zeros(length(x), 2);
% 
% F = zeros(length(x), modes * length(indexLocalEnrichedNodes));
% G = zeros(length(x), modes * length(indexLocalEnrichedNodes));
% 
% PODCoefficients = problem.reductionOperator(...
%     (e-1)*(floor(refinedNodes/problem.XN))+1:(e-1)*(floor(refinedNodes/problem.XN)) +...
%      ceil(refinedNodes/problem.XN),:);
% 
% for k=1:length(x)
%     
%     [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
%     
%     [F(k,:), G(k,:)] = PODModesAndDerivativesGaussIntegration(localCoords(k), modes, PODCoefficients,...
%         subDomainShapeFunctionCoefficients, subDomainIndex, indexLocalEnrichedNodes);
% end
% 
% if derivative == 0
%     for i=1:length(localCoords)
%         projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
%     end
% else
%     for i=1:length(localCoords)
%         projectionOperator(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
%     end
% end
% 
% % if length(indexLocalEnrichedNodes) == 1 %lhs
% %     coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
% %         solutionCoefficients(problem.N  + (e - 1) * length(indexLocalEnrichedNodes) * modes + 2 : problem.N...
% %         + (e) * length(indexLocalEnrichedNodes) * modes + 1 )];
% % else
% %     coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
% %         solutionCoefficients(problem.N  + (e - 2) * modes + 2 : problem.N...
% %         + (e - 2) * modes + 2 + length(indexLocalEnrichedNodes)  * modes - 1 )];
% % end
% 
% if length(indexLocalEnrichedNodes) == 1 %lhs
%     coefficients = [solutionCoefficients((problem.N - problem.XN) + e + problem.p-1 : (problem.N - problem.XN) + e + 1); ...
%         solutionCoefficients(problem.N  + (e - 1) * length(indexLocalEnrichedNodes) * modes + 2 : problem.N...
%         + (e) * length(indexLocalEnrichedNodes) * modes + problem.p )];
% else
%     coefficients = [solutionCoefficients((problem.N - problem.XN) + e + problem.p-1: (problem.N - problem.XN) + e + 1); ...
%         solutionCoefficients(problem.N  + (e - 2) * modes + 2 : problem.N...
%         + (e - 2) * modes + length(indexLocalEnrichedNodes)  * modes + problem.p )];
% end
% 
% projectedCoefficients = projectionOperator * coefficients ;
% 
% if derivative == 1
%     for i=1:length(localCoords)
%         projectionOperator_Temp(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
%     end
%     
%     temperature = projectionOperator_Temp * coefficients ;
%     for i=1:numel(x)
%         projectedCoefficients(i) = projectedCoefficients(i) .* ...
%             problem.k(globalCoords(i), t, temperature(i));
%     end
% end
% 
% end