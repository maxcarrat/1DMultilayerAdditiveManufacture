function [temperaturePostProcessing, heatFluxes, internalEnergy, modes]= backwardEulerXFEM(coords, postProcessingCoords, rhs, initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector,...
    refinementDepth, numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes)
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

refinedTemperatureSolutions = zeros(2^refinementDepth+1, 1);
% refinedTemperatureSolutions(:) = initialTemperature;

localRefinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

previousMesh = coords;

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for layer = 1:numberOfTrainingLayers
 
    %Refinement
    refinedMesh = refineLayer(coords, refinementDepth, layer, numberOfLayers);
%     refinedMesh = refineMesh(coords, refinementDepth, layer, numberOfTrainingLayers);
    
    %Project old solution onto the new mesh
    if layer > 1
        refinedTemperatureSolutions = L2projection(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature);
    else
       refinedTemperatureSolutions = zeros(size(refinedMesh, 2), 1);
       localRefinedTemperatureSolutions = zeros(size(refinedMesh, 2), timeSteps);
    end
    
    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            k, heatCapacity, currentTime);
        
        %Backward Euler Scheme
        [M, K, f] = assembly(poissonTransientProblem);
        [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, refinedTemperatureSolutions, timeStepSize);
        temperatureIncrement = LHS\RHS;
        
        %Update and merge temperature into global domain
        refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
        
        localRefinedTemperatureSolutions(:, t+1) = getLayerSolution(refinedTemperatureSolutions, layer, numberOfLayers, coords);
        previousMesh = refinedMesh;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
        internalEnergy(t+1) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
        
    end
end

%% Generate the reduced basis

[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,5:numberOfTrainingLayers), numberOfPODModes);

%% Enriched mesh using RB


for layer = (numberOfTrainingLayers+1):numberOfLayers
    
    [activeMesh, ~] = getLayerActiveCoords(coords, layer, numberOfLayers);
    activeMeshSize = numel(activeMesh);
    
    %Project old solution onto the new mesh
    previousSolution = [refinedTemperatureSolutions(1:activeMeshSize-2);refinedTemperatureSolutions(end); 0.0];
        

    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * (t-1);
        
        %Generate the Poisson problem at timeStep t on the coarse active mesh
        poissonTransientProblem = poissonProblemTransient(activeMesh, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            k, heatCapacity, currentTime);
        
        disp(' Solve Global Problem ');
        
        %Solve Global/Coarse problem
        temperatureSolutionsGlobal = solveGlobalProblem(previousSolution, poissonTransientProblem, timeStepSize);
        
        %Generate Local problem
        poissonTransientProblemEnriched = poissonProblemTransientEnriched(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
            rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
        
        %Project global solution onto the enriched modal space
        temperatureSolutionsProjected = L2projectionEnriched(poissonTransientProblemEnriched,previousSolution,...
            activeMesh, coords, modes, initialTemperature);
        temperatureSolutionsProjected(1) = temperatureSolutionsGlobal(end-1);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        temperatureSolutionsEnriched = solveLocalProblem(temperatureSolutionsProjected, poissonTransientProblemEnriched, timeStepSize, modes);
        
        %Update solutions
        modesOffset = modes*2+2;
        temperatureSolutions = zeros(size(temperatureSolutionsGlobal,1)-2+modesOffset,1);
        temperatureSolutions(1:end-modesOffset+2) = temperatureSolutionsGlobal(1:end);
        temperatureSolutions(end-modesOffset+2:end) = temperatureSolutionsEnriched(2:end);
        
        refinedTemperatureSolutions = temperatureSolutionsGlobal;
        
        temperaturePostProcessing(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh,...
            poissonTransientProblemEnriched, temperatureSolutions, temperatureSolutionsGlobal, 0) ;
        heatFluxes(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh, poissonTransientProblemEnriched,...
            temperatureSolutions, temperatureSolutionsGlobal, 1);
        
        %     [~, K, ~] = assembly(poissonTransientProblem);
        %     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
        
    end

end






