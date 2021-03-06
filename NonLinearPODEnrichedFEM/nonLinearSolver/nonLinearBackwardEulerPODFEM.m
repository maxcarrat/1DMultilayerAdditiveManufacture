function [temperaturePostProcessing, heatFluxes, internalEnergy, modes]= nonLinearBackwardEulerPODFEM(coords, postProcessingCoords, rhs, initialTemperature, leftDirichletBoundaryConditionValue,...
    rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refinementDepth, numberOfTrainingTimeSteps, maxIterations, tolerance, numberOfModes)
% nonLinearBackwardEulerPODFEM computes the 1D h-FEM numerical solution of a non-linear boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   rhs = external heat source
%   initialTemperature = initial temeprature of each new layer
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

refinedTemperatureSolutions = zeros(2^refinementDepth+1, 1);
refinedTemperatureSolutions(:) = initialTemperature;

localRefinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

previousMesh = coords;

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for t = 2:numberOfTrainingTimeSteps
    formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Refinement
    refinedMesh = refineMesh(coords, refinementDepth, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonNonLinearProblemTransient(refinedMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    
    %Project old solution onto the new mesh
    if (norm(refinedTemperatureSolutions - initialTemperature))~=0.0
        refinedTemperatureSolutions = L2projection(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature);
    end
    
    %Solve the problem using Newton-Raphson scheme
    refinedTemperatureSolutions = newtonRaphsonGlobalProblem(refinedTemperatureSolutions, poissonTransientProblem, tolerance, maxIterations, timeStepSize);
    
    %Update and merge temperature into global domain
    mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
    localRefinedTemperatureSolutions(:, t) = refinedTemperatureSolutions(t-1:end);
    previousMesh = refinedMesh;
    
    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
    %     internalEnergy(t) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
end

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,5:numberOfTrainingTimeSteps), numberOfModes);

%% Enriched mesh using RB

for t = (numberOfTrainingTimeSteps+1):timeSteps
    
    formatSpec = 'Backward Euler Time Step (Reduced Basis): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Get active mesh
    activeMesh = getActiveCoordinates(coords, t, timeSteps);
    activeMeshSize = numel(activeMesh);
    
    %Generate the Poisson problem at timeStep t on the coarse active mesh
    poissonTransientProblem = poissonNonLinearProblemTransient(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    
    %Project old solution onto the new mesh
    previousSolution = [refinedTemperatureSolutions(1:activeMeshSize-2);refinedTemperatureSolutions(end); 0.0];
    
    disp(' Solve Global Problem');
    
    %Solve the global problem using Newton-Raphson scheme on the coarse
    %mesh
    temperatureSolutionsGlobal = newtonRaphsonGlobalProblem(previousSolution, poissonTransientProblem, tolerance, maxIterations, timeStepSize);
    
    %Generate the Poisson problem at timeStep t on the coarse active mesh
    poissonTransientProblemEnriched = poissonNonLinearProblemTransientPODEnrichment(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
    
    %Project old solution onto the new mesh
    temperatureSolutionsProjected = L2projectionEnriched(poissonTransientProblemEnriched, previousSolution,...
        activeMesh, coords, modes, initialTemperature);
    temperatureSolutionsProjected(1) = temperatureSolutionsGlobal(end);
    
    disp(' Solve Local Enriched Problem');
    
    %Solve the local problem using Newton-Raphson scheme using POD reduced
    %basis
    temperatureSolutionsEnriched = newtonRaphsonEnrichedProblem(temperatureSolutionsProjected,...
        poissonTransientProblemEnriched, modes, tolerance, maxIterations, timeStepSize);
    
    %Update solutions
    %Update solutions
    modesOffset = modes*2+2;
    temperatureSolutions = zeros(size(temperatureSolutionsGlobal,1)-2+modesOffset,1);
    temperatureSolutions(1:end-modesOffset+2) = temperatureSolutionsGlobal(1:end);
    temperatureSolutions(end-modesOffset+2:end) = temperatureSolutionsEnriched(2:end);
        
    refinedTemperatureSolutions = temperatureSolutionsGlobal;
    
    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh,...
        poissonTransientProblemEnriched, temperatureSolutions, temperatureSolutionsGlobal, 0) ;
    heatFluxes(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh, poissonTransientProblemEnriched,...
        temperatureSolutions, temperatureSolutionsGlobal, 1);
    
    %     [~, K, ~] = assembly(poissonTransientProblem);
    %     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
    
end

end






