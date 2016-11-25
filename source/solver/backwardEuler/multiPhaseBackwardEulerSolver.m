function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs, ...
    localRefinedTemperatureSolutions, solutionReductionOperator]= multiPhaseBackwardEulerSolver(coords, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance, maxIterations, numberOfRefinedElementsToBeKept,...
    refinementDepth, PODRefinementDepth, numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes,...
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

refinedTemperatureSolutions = linspace(0.0,...
    0.0, 2^refinementDepth);

convergence = 0.0;

computationalTime = [];
localRefinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);
localHeatFluxes = zeros(2^refinementDepth+1, timeSteps);

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
        [refinedTemperatureSolutions, convergenceFlag] = solveMultiPhaseProblem( poissonTransientProblem,...
            currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
        
%         localRefinedTemperatureSolutions(:, t+1) = getLayerSolution(refinedTemperatureSolutions, layer, numberOfLayers, coords);
        localRefinedTemperatureSolutions(:, t+1) = getSolution(refinedTemperatureSolutions,...
            layer, numberOfLayers, coords, refinementDepth, numberOfRefinedElementsToBeKept);
        localHeatFluxes(:, t+1) = getSolution(evaluateNumericalResults(refinedMesh, currentTime,...
            poissonTransientProblem, mergedTemperature, 1),...
            layer, numberOfLayers, coords, refinementDepth, numberOfRefinedElementsToBeKept);


        previousMesh = refinedMesh;
        
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, mergedTemperature, 1);
        %         internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
        
        DOFs = numel(refinedTemperatureSolutions);
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,3*numberOfLayersTimeSteps+2:numberOfTrainingLayers*numberOfLayersTimeSteps+1), numberOfPODModes);
% [fluxesReductionOperator, modes] = properOrthogonalDecomposition(localHeatFluxes(:,1:numberOfTrainingLayers*numberOfLayersTimeSteps), numberOfPODModes);

%% Enriched mesh using RB
if layer ~= numberOfLayers
%Project old solution onto the new mesh
refinedDOFs = 2^refinementDepth + 1;
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
    [activeMesh, numberOfActiveElementsLayer] = refineLayerEnriched(activeMesh, numberOfActiveElementsLayer,...
        refinementDepth, layer,...
        numberOfTrainingLayers, PODRefinementDepth, numberOfRefinedElementsToBeKept);

    poissonTransientProblemEnriched = poissonProblemXFEM(activeMesh,numberOfActiveElementsLayer, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth,...
        PODRefinementDepth, solutionReductionOperator);
    
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
            rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue, k, heatCapacity, currentTime, refinementDepth,...
            PODRefinementDepth, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        [temperatureSolutions, convergenceFlag] = solveMultiPhaseXFEMProblem( poissonTransientProblemEnriched, currentTime,...
            timeStepSize,integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
        
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
end

    if convergence < 1
        disp('The analysis always converged')
    else
        disp('The analysis did not always converged !!!')
        fprintf(num2str(convergence));
    end
end