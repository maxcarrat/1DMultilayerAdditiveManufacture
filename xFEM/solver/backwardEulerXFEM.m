function [temperaturePostProcessing, heatFluxes, internalEnergy]= backwardEulerXFEM(coords, postProcessingCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refinementDepth, numberOfTrainingTimeSteps)
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

refinedTemperatureSolutions = zeros(2^refinementDepth+1, timeSteps);
localRefinedTemperatureSolutions = zeros(2^refinementDepth, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
  for t = 2:numberOfTrainingTimeSteps    
    formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Refinement
    refinedMesh = refineMesh(coords, refinementDepth, t, timeSteps);
    
    %Project old solution onto the new mesh
    refinedTemperatureSolutions = L2projection(refinedTemperatureSolutions, refinedMesh, coords);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [~, Kunconstrained, ~] = assembly(poissonTransientProblem);
    [M, K, f] = assemblyAndApplyStrongBCs(poissonTransientProblem);
    
    %Backward Euler Scheme
    RHS = timeStepSize * (f - K * refinedTemperatureSolutions);
    LHS = M + timeStepSize * K;
    
    temperatureIncrement = LHS\RHS;
    
    %Update temperature 
    refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
    localRefinedTemperatureSolutions(:, t) = refinedTemperatureSolutions(t-1:end);    

    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, refinedTemperatureSolutions, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, refinedTemperatureSolutions, 1);
    internalEnergy(t) = refinedTemperatureSolutions'*Kunconstrained*refinedTemperatureSolutions;
  end
  
 %% Generate the reduced basis
 
  [solutionReductionOperator] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,2:numberOfTrainingTimeSteps));
  
 %% Enriched mesh using RB

  enrichmentCoefficents = zeros(2+1, 1);

  for t = (numberOfTrainingTimeSteps+1):timeSteps
    
    formatSpec = 'Backward Euler Time Step (Reduced Basis): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Get active mesh
    activeMesh = getActiveCoordinates(coords, t, timeSteps);
    modes = size(solutionReductionOperator, 2);
    
    %Project old solution onto the new mesh
    temperatureSolutions = L2projectionEnriched(refinedTemperatureSolutions, activeMesh, refinedMesh, modes);
    
    %Solve Global/Coarse problem
    poissonTransientProblem = poissonProblemTransient(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
   temperatureSolutions = solveGlobalProblem(temperatureSolutions, poissonTransientProblem, timeStepSize);
    
    %Solve Local problem

    poissonTransientProblemEnriched = poissonProblemTransientEnriched(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
   temperatureSolutionsEnriched = solveLocalProblem(temperatureSolutions, enrichmentCoefficents, poissonTransientProblemEnriched, timeStepSize);
      
    
    %Update and merge temperature into global domain
    temperatureSolutions(end-2:end) = temperatureSolutions(end-2:end) + temperatureSolutionsEnriched;
    refinedMesh = activeMesh;
    refinedTemperatureSolutions = temperatureSolutions;
    enrichmentCoefficents = enrichmentCoefficents + temperatureSolutionsEnriched;

%     %Generate the refined Poisson problem at timeStep t and Assembly the
%     %linear system
%     poissonTransientProblemEnriched = poissonProblemTransientEnriched(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
%         rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
%     [~, Kunconstrained, ~] = assemblyEnrichedProblem(poissonTransientProblemEnriched);
%     [M, K, f] = applyBCs(poissonTransientProblemEnriched);
%     
%     %Backward Euler Scheme
%     RHS = timeStepSize * (f - K * refinedTemperatureSolutions);
%     LHS = M + timeStepSize * K;
%     
%     temperatureIncrement = LHS\RHS;
%     
%     %Update and merge temperature into global domain
%     refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
%     refinedMesh = activeMesh;
    
    %Post-Processing
%     temperaturePostProcessing(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh(end-1:end),...
%         poissonTransientProblemEnriched, temperatureSolutions, 0) ;
%     heatFluxes(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh(end-1:end), poissonTransientProblemEnriched,...
%         temperatureSolutions, 1);

    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords,...
        poissonTransientProblemEnriched, temperatureSolutions, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblemEnriched,...
        temperatureSolutions, 1);
%     internalEnergy(t) = temperatureSolutions'*Kunconstrained*temperatureSolutions;
    
  end
  
end






