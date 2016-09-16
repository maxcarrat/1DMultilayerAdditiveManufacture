function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime]= nonLinearBackwardEulerGaussIntegrationSolver(coords, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance, maxIterations, numberOfRefinedElementsToBeKept,...
    refinementDepth, PODRefinementDepth, numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes, integrationOrder, integrationModesOrder)
% BackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
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
temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

refinedTemperatureSolutions = linspace(0.0,...
    0.0, 2^refinementDepth);

convergence = 0.0;

computationalTime = [];
localRefinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

previousMesh = coords;

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for layer = 1:numberOfTrainingLayers
    
    %Refinement
    %     refinedMesh = refineLayer(coords, refinementDepth, layer, numberOfLayers);
    refinedMesh = refineLayerNoCoarse(coords, refinementDepth, layer, numberOfLayers, numberOfRefinedElementsToBeKept);
    
    %Project old solution onto the new mesh
    if layer > 1
        refinedTemperatureSolutions = L2projectionLinearDistribution(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature, refinementDepth);
    else
        refinedTemperatureSolutions = zeros(size(refinedMesh, 2), 1);
        localRefinedTemperatureSolutions = zeros(size(refinedMesh, 2), timeSteps);
    end
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime);
        
        %Update and merge temperature into global domain
        [refinedTemperatureSolutions, convergenceFlag] = solveNonLinearProblemGaussIntegration( poissonTransientProblem,...
            currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
        
%         localRefinedTemperatureSolutions(:, t+1) = getLayerSolution(refinedTemperatureSolutions, layer, numberOfLayers, coords);
        localRefinedTemperatureSolutions(:, t+1) = getSolution(refinedTemperatureSolutions,...
            layer, numberOfLayers, coords, refinementDepth, numberOfRefinedElementsToBeKept);

        previousMesh = refinedMesh;
        
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 1);
        %         internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,1:numberOfTrainingLayers*numberOfLayersTimeSteps), numberOfPODModes);

%% Enriched mesh using RB
if layer ~= numberOfLayers
%Project old solution onto the new mesh
refinedDOFs = 2^refinementDepth;
% temperatureSolutions = [refinedTemperatureSolutions(1:end-refinedDOFs);...
%     refinedTemperatureSolutions(end); initialTemperature];

% temperatureSolutions = [refinedTemperatureSolutions; initialTemperature];

% project mesh onto coarse mesh
switch numberOfRefinedElementsToBeKept
    case(1)
        refinedMeshEnriched = [refinedMesh(1:end-refinedDOFs)';...
            refinedMesh(end)'; coords(layer+2)];
    case(5)
        refinedMeshEnriched = [refinedMesh(1);...
            refinedMesh(end-4*refinedDOFs:end)'; coords(layer+2)];
    otherwise
        disp('Case not implemented yet!');
end
% [previousMesh, ~] = getLayerActiveCoords(coords, layer, numberOfLayers);
previousMesh = refinedMeshEnriched;

temperatureSolutions = projectOntoEnrichedMesh(poissonTransientProblem, refinedTemperatureSolutions,...
    modes, refinedMeshEnriched, refinedMesh, PODRefinementDepth, initialTemperature);
end
for layer = (numberOfTrainingLayers+1):numberOfLayers
    
    [activeMesh, numberOfActiveElementsLayer] = getLayerActiveCoords(coords, layer, numberOfLayers);
    activeMesh = refineLayerEnriched(activeMesh, numberOfActiveElementsLayer, refinementDepth, layer,...
        numberOfTrainingLayers, PODRefinementDepth,  numberOfRefinedElementsToBeKept-1);

    poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth, PODRefinementDepth, solutionReductionOperator);
    
    %Project global solution onto the enriched modal space
%     temperatureSolutions = eXtendedProjection(poissonTransientProblemEnriched,temperatureSolutions,...
%         modes, initialTemperature);
    temperatureSolutions = eXtendedProjectionNoCoarse(poissonTransientProblemEnriched,temperatureSolutions,...
        modes, activeMesh, previousMesh, PODRefinementDepth, initialTemperature);
    
    timeToGenerateAndSolveTheSystem = 0.0;

    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate Local problem
        poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
            rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth, PODRefinementDepth, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        [temperatureSolutions, convergenceFlag] = solveNonLinearEnrichedProblem( poissonTransientProblemEnriched, currentTime,...
            timeStepSize,integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
        
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;

        %Post-Processing
        temperaturePostProcessing(:, t+1) = postProcessingProjection(postProcessingCoords, poissonTransientProblemEnriched, temperatureSolutions,...
            modes, 0.0);
        heatFluxes(:, t+1) = postProcessingProjection(postProcessingCoords, poissonTransientProblemEnriched,temperatureSolutions,...
            modes, 1);
        
        %     [~, K, ~] = assembly(poissonTransientProblem);
        %     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
        
    end
    
    previousMesh = activeMesh;
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

    if convergence <= 1
        disp('The analysis always converged')
    else
        disp('The analysis did not always converged !!!')
    end
end

function [ projectedCoefficients ] = eXtendedProjection(problem, previousTemperature,...
    modes, initialTemperature)
% EXTENDEDPROJECTION project the previous solution onto the updated mesh at the
% new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   coarseMesh = initial/coarse mesh configuration
%   modes = number of enrichment modes
%   initialTemperature = initial temperature of the powder

projectedCoefficients = zeros(problem.N + 1 +...
    modes * (2*problem.XN - 1), 1);

for e=1:problem.N

    if e == problem.N % last element
        % lhs node
        projectedCoefficients(e) = previousTemperature(e);
        % rhs node
        projectedCoefficients(e+1) = initialTemperature;
        
    else     
        projectedCoefficients(e) = previousTemperature(e);
    end
end


end

function [ projectedCoefficients ] = eXtendedProjectionNoCoarse(problem, previousTemperature,...
    modes, updatedMesh, previousMesh, PODRefinementDepth, initialTemperature)
% EXTENDEDPROJECTIONNOCOARSE project the previous solution onto the updated mesh at the
% new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   coarseMesh = initial/coarse mesh configuration
%   modes = number of enrichment modes
%   initialTemperature = initial temperature of the powder

projectedCoefficients = zeros(problem.N + 1 +...
    modes * (2*problem.XN - 1), 1);

% projectedCoefficients(1:problem.N+1) = L2projectionLinearDistribution(problem, previousTemperature, updatedMesh,...
%     previousMesh, initialTemperature, 0);

projectedCoefficients = projectOntoEnrichedMesh(problem, previousTemperature, modes, updatedMesh,...
    previousMesh, PODRefinementDepth, initialTemperature);

% for e=1:problem.N
% 
%     if e == problem.N % last element
%         % lhs node
%         projectedCoefficients(e) = previousTemperature(e);
%         % rhs node
%         projectedCoefficients(e+1) = initialTemperature;        
%     else     
%         projectedCoefficients(e) = previousTemperature(e);
%     end
% end

end

function [ projectedCoefficients ] = postProcessingProjection(x, problem, solutionCoefficients,...
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

projectedCoefficients(x>=X1 & x<=X2) = postProcessingProjectionElement(1, x(x>=X1 & x<=X2), problem,... % first element
    solutionCoefficients, modes, derivative);

% loop over elements
for e=2:problem.N
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    projectedCoefficients(x>X1 & x<=X2) = postProcessingProjectionElement(e, x(x>X1 & x<=X2), problem,...
        solutionCoefficients, modes, derivative);
    
end

end

function [ projectedCoefficients ] = postProcessingProjectionElement(e, x, problem, solutionCoefficients,...
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
    
    if e == problem.N - problem.XN +1
        indexLocalEnrichedNodes = 2; %rhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    integrationDomain = linspace(-1, +1, 2^problem.refinementDepth + 1);
    subDomainShapeFunctionCoefficients = linspace(-1, 1, 2^problem.refinementDepth + 1);
    
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    localCoordinates = mapGlobalToLocal( x, X1, X2);
    projectedCoefficients=zeros(size(localCoordinates));
    
    projectedCoefficients(localCoordinates>=Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(1,...
        localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), problem, integrationDomain,...
        solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
    
    for integrationSubDomainIndex=2:length(integrationDomain)-1
        
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
        
        projectedCoefficients(localCoordinates>Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(integrationSubDomainIndex,...
            localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), problem, integrationDomain,...
            solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
          
  
    end
    
%     projectedCoefficients = projectedCoefficients .* (2/(Xi2-Xi1)) ^ derivative;
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
    
    projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:));
    
    projectedCoefficients = projectedCoefficients .* (2/(X2-X1)) ^ derivative;
    
end


end


function [ projectedCoefficients ] = postProcessingProjectionSubElements(e, x, problem, integrationDomain, solutionCoefficients,...
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
projectionOperator = zeros(numberOfProjectionPoints, length(solutionCoefficients(problem.N:end)) );

% X1 = integrationDomain(e);
% X2 = integrationDomain(e+1);
% localCoords =  mapGlobalToLocal( x, X1, X2 );
localCoords = x;

N = zeros(length(x), 2);
B = zeros(length(x), 2);

F = zeros(length(x), modes);
G = zeros(length(x), modes);

for k=1:length(x)
    
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    
    [F(k,:), G(k,:)] = PODModesAndDerivativesGaussIntegration(localCoords(k), modes, problem.reductionOperator,...
        subDomainShapeFunctionCoefficients, e, indexLocalEnrichedNodes);
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

projectedCoefficients = projectionOperator * solutionCoefficients(problem.N:end);
 
end