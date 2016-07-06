function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime] = nonLinearBackwardEulerNoCoarseningGaussIntegrationSolver(coords, postProcessingCoords, rhs, initialTemperature, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, neumannBoundaryConditionValue, k, heatCapacity, timeVector, tolerance, maxIterations,...
    refinementDepth, numberOfLayersTimeSteps, numberOfLayers, integrationOrder)
% NonLinearBackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   rhs = external heat source
%   leftDirichletBoundaryConditionValue
%   rightDirichletBoundaryConditionValue
%   k = conductivity of the material
%   heatCapacity = heat capacity of the material
%   timeVector = vector of time steps for Backward Euler implicit scheme
%   refinementDepth = depth of the refined mesh

timeSteps=size(timeVector,2);

timeStepSize=max(timeVector)/( timeSteps );
temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

computationalTime = [];

refinedTemperatureSolutions = linspace(0.0,...
    0.0,2^refinementDepth);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

previousMesh = coords;
convergence = 0.0;

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Backward Euler
for layer = 1:numberOfLayers
    
    %Refinement
    %generate refined mesh for the new active configuration
    [activeMesh, numberOfElementsPerLayer] = getLayerActiveCoords(coords, layer, numberOfLayers);
    refinedActiveElement = linspace(activeMesh(end-numberOfElementsPerLayer), activeMesh(end), numberOfElementsPerLayer*2^(refinementDepth) + 1);
    
    %Project old solution onto the new mesh
    if layer > 1
        refinedMesh = [refinedMesh(1:end-1) refinedActiveElement];
        refinedTemperatureSolutions = L2projectionLinearDistribution(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature, refinementDepth);
    else
        refinedMesh = refinedActiveElement;
        refinedTemperatureSolutions = linspace(0.0,...
    0.0,size(refinedMesh, 2))';
    end
    
    timeToGenerateAndSolveTheSystem = 0.0;

    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryConditionValue, k, heatCapacity, currentTime);
        
        %Update and merge temperature into global domain
        [refinedTemperatureSolutions, convergenceFlag] = solveNonLinearProblemGaussIntegration( poissonTransientProblem,...
            currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
        previousMesh = refinedMesh;
        
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
        %         internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
    
    if convergence <= 1
        disp('The analysis always converged')
    else
        disp('The analysis did not always converged !!!')
    end
end

