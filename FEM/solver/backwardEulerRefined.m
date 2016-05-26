function [temperaturePostProcessing, heatFluxes, refinedTemperatureSolutions, refinedHeatFluxes, internalEnergy]= backwardEulerRefined(coords, postProcessingCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refinementDepth)
% BACKWARDEULERREFINED computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );

temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

refinedTemperatureSolutions = zeros(2^refinementDepth, 1);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
refinedHeatFluxes = zeros(2^refinementDepth+1, timeSteps);
previousMesh = coords;
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

  for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Refinement
    refinedMesh = refineMesh(coords, refinementDepth, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assembly(poissonTransientProblem);
    
    %Project old solution onto the new mesh
    if (norm(refinedTemperatureSolutions))~=0.0
        refinedTemperatureSolutions = L2projection(poissonTransientProblem, refinedTemperatureSolutions, refinedMesh, previousMesh);
    end
    
    %Backward Euler Scheme
    [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, refinedTemperatureSolutions, timeStepSize);
    temperatureIncrement = LHS\RHS;
    
    %Update and merge temperature into global domain
    refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
    mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, size(coords, 2));
    previousMesh= refinedMesh;
    
    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
    internalEnergy(t) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
  end
  
end






