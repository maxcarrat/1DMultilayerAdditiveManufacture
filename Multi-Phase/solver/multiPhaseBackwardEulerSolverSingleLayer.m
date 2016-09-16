function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]...
    = multiPhaseBackwardEulerSolverSingleLayer(coords, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance, maxIterations, numberOfRefinedElementsToBeKept,...
    refinementDepth, PODRefinementDepth, numberOfTrainingTimeStep, numberOfPODModes,...
    integrationOrder, integrationModesOrder)
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

convergence = 0.0;

computationalTime = [];
localHeatFluxes = zeros(2^refinementDepth+1, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
refinedMesh = linspace(coords(end-1), coords(end), 2^(refinementDepth)+ 1);
refinedTemperatureSolutions = zeros(size(refinedMesh, 2), 1);
localRefinedTemperatureSolutions = zeros(size(refinedMesh, 2), 1);

timeToGenerateAndSolveTheSystem = 0.0;

for iTime = 1:numberOfTrainingTimeStep
    tic
    t = iTime;
    
    formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
    fprintf(formatSpec, t)
    
    currentTime = timeStepSize * t;
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
        leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
        neumannBoundaryconditionValue, k, heatCapacity, currentTime);
    
    %Update and merge temperature into global domain
    [refinedTemperatureSolutions, convergenceFlag] = solveMultiPhaseProblem( poissonTransientProblem,...
        currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
    mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
    
    localRefinedTemperatureSolutions(:, t+1) = getSolution(refinedTemperatureSolutions,...
        1, 1, coords, refinementDepth, numberOfRefinedElementsToBeKept);
    
    localHeatFluxes(:, t+1) = getSolution(evaluateNumericalResults(refinedMesh, currentTime,...
        poissonTransientProblem, mergedTemperature, 1),...
        1, 1, coords, refinementDepth, numberOfRefinedElementsToBeKept);
    
       
    convergence = convergence + convergenceFlag;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    
    %Post-Processing
    temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 1);
    %         internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
    
    DOFs = numel(refinedTemperatureSolutions);
    
end
computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,1:numberOfTrainingTimeStep), numberOfPODModes);

%% Enriched mesh using RB
%Project old solution onto the new mesh
refinedDOFs = 2^refinementDepth + 1;

%Project mesh onto coarse mesh
switch numberOfRefinedElementsToBeKept
    case(1)
        refinedMeshEnriched = refinedMesh;
    case(5)
        refinedMeshEnriched = [refinedMesh(1);...
            refinedMesh(end-4*refinedDOFs:end)'; 1];
    otherwise
        disp('Case not implemented yet!');
end

temperatureSolutions = projectOntoEnrichedMesh(poissonTransientProblem, refinedTemperatureSolutions,...
    modes, refinedMeshEnriched, coords, PODRefinementDepth, initialTemperature);

activeMesh = refinedMeshEnriched;

timeToGenerateAndSolveTheSystem = 0.0;

for iTime = (iTime+1):timeSteps
    tic
    t = iTime;
    
    formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
    fprintf(formatSpec, t)
    
    currentTime = timeStepSize * t;
    
    %Generate Local problem
    activeMeshSize = 2.^PODRefinementDepth;
    poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh, activeMeshSize, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth,...
        PODRefinementDepth, solutionReductionOperator);
    
    disp(' Solve Local Enriched Problem ');
    
    %Solve Local problem enriched
    [temperatureSolutions, convergenceFlag] = solveMultiPhaseXFEMProblem( poissonTransientProblemEnriched, currentTime,...
        timeStepSize, integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
    
    convergence = convergence + convergenceFlag;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    
    %Post-Processing
    temperaturePostProcessing(:, t+1) = postProcessingProjection(postProcessingCoords, currentTime, poissonTransientProblemEnriched, temperatureSolutions,...
        modes, 0.0);
    heatFluxes(:, t+1) = postProcessingProjection(postProcessingCoords, currentTime, poissonTransientProblemEnriched,temperatureSolutions,...
        modes, 1);
    
    %     [~, K, ~] = assembly(poissonTransientProblem);
    %     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
    
    DOFs = numel(temperatureSolutions);
    
end

previousMesh = activeMesh;

computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];


if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
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

projectedCoefficients = projectOntoEnrichedMesh(problem, previousTemperature, modes, updatedMesh,...
    previousMesh, PODRefinementDepth, initialTemperature);

end

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
    
    [F(k,:), G(k,:)] = PODModesAndDerivativesGaussIntegration(problem, localCoords(k), modes, PODCoefficients,...
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