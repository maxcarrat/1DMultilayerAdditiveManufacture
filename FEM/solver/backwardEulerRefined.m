function [temperaturePostProcessing, heatFluxes, refinedTemperatureSolutions, refinedHeatFluxes, internalEnergy]= backwardEulerRefined(coords, postProcessingCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refinementDepth)
% BACKWARDEULERREFINED computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );

temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

refinedTemperatureSolutions = zeros(2^refinementDepth, timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
refinedHeatFluxes = zeros(2^refinementDepth+1, timeSteps);

internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

  for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Refinement
    refinedMesh = refineMesh(coords, refinementDepth, t, timeSteps);
    
    %Project old solution onto the new mesh
    refinedTemperatureSolutions = L2projection(refinedTemperatureSolutions, refinedMesh, coords);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assemblyAndApplyStrongBCs(poissonTransientProblem);
    
    %Backward Euler Scheme
    RHS = timeStepSize * (f - K * refinedTemperatureSolutions);
    LHS = M + timeStepSize * K;
    
    temperatureIncrement = LHS\RHS;
    
    %Update and merge temperature into global domain
    mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions + temperatureIncrement, size(refinedTemperatureSolutions,1));
    refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
    
    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
    internalEnergy(t) = mergedTemperature'*K*mergedTemperature;
  end
  
end






