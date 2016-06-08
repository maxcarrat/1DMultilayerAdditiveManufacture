function [temperaturePostProcessing, heatFluxes, internalEnergy, modes]= backwardEulerXFEM(coords, postProcessingCoords, rhs, initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector,...
    refinementDepth, numberOfTrainingTimeSteps, numberOfPODModes)
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
  for t = 2:numberOfTrainingTimeSteps    
      
    formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Refinement
    refinedMesh = refineMesh(coords, refinementDepth, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
        leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
        k, heatCapacity, currentTime);
    
    %Project old solution onto the new mesh
%     previousSolution = refinedTemperatureSolutions;
    if (norm(refinedTemperatureSolutions))~=0.0
        refinedTemperatureSolutions = L2projection(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh, initialTemperature);
    end
    
    %Backward Euler Scheme
    [M, K, f] = assembly(poissonTransientProblem);
    [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, refinedTemperatureSolutions, timeStepSize);
    temperatureIncrement = LHS\RHS;
    
    %Update and merge temperature into global domain
    refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
    mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));

    localRefinedTemperatureSolutions(:, t) = refinedTemperatureSolutions(t-1:end);    
    previousMesh = refinedMesh;

    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
    internalEnergy(t) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
  end
  
 %% Generate the reduced basis
 
  [solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,4:numberOfTrainingTimeSteps), numberOfPODModes);
  
 %% Enriched mesh using RB


  for t = (numberOfTrainingTimeSteps+1):timeSteps
    
    formatSpec = 'Backward Euler Time Step (Reduced Basis): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Get active mesh
    activeMesh = getActiveCoordinates(coords, t, timeSteps);
    
    %Project old solution onto the new mesh
    previousSolution = refinedTemperatureSolutions;
    refinedTemperatureSolutions = L2projection(poissonTransientProblem, previousSolution, activeMesh, refinedMesh, initialTemperature);

    %Solve Global/Coarse problem
    poissonTransientProblem = poissonProblemTransient(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    temperatureSolutionsGlobal = solveGlobalProblem(refinedTemperatureSolutions, poissonTransientProblem, timeStepSize);
    
    %Generate Local problem
    poissonTransientProblemEnriched = poissonProblemTransientEnriched(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
   
    %Project global solution onto the enriched modal space
    temperatureSolutionsProjected = L2projectionEnriched(poissonTransientProblemEnriched,refinedTemperatureSolutions,...
        activeMesh, coords, modes, initialTemperature);
    
    %Solve Local problem enriched
    temperatureSolutionsEnriched = solveLocalProblem(temperatureSolutionsProjected, poissonTransientProblemEnriched, timeStepSize, modes);
    
    %Update solutions
    modesOffset = modes*2+2;
    temperatureSolutions = zeros(size(temperatureSolutionsGlobal,1)-2+modesOffset,1);
    temperatureSolutions(1:end-modesOffset+2) = temperatureSolutionsGlobal(1:end);
    temperatureSolutions(end-modesOffset+2:end) = temperatureSolutionsEnriched(2:end);
    
%    temperatureSolutionsEnriched = reducedBasisCoefficients;
%     
%     %Update and merge temperature into global domain
% %     temperatureSolutions(end-numberOfModesSupports:end) = temperatureSolutions(end-numberOfModesSupports:end) + temperatureSolutionsEnriched;
%     temperatureSolutions(end-numberOfModesSupports:end) = temperatureSolutionsEnriched;

    refinedMesh = activeMesh;
    refinedTemperatureSolutions = temperatureSolutions;

    temperaturePostProcessing(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh,...
        poissonTransientProblemEnriched, temperatureSolutions, temperatureSolutionsGlobal, 0) ;
    heatFluxes(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh, poissonTransientProblemEnriched,...
        temperatureSolutions, temperatureSolutionsGlobal, 1);
    
%     [~, K, ~] = assembly(poissonTransientProblem);
%     internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
    
  end
  
end






