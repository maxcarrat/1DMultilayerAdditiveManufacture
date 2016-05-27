function [ temperaturePostProcessing, heatFluxes, internalEnergy ] = nonLinearBackwardEuler( coords, postProcessingCoords, rhs,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, maxIterations, tolerance )
% nonLinearBackwardEuler computes the 1D h-FEM numerical solution of a non-linear boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   rhs = external heat source
%   leftDirichletBoundaryConditionValue
%   rightDirichletBoundaryConditionValue
%   k = conductivity of the material
%   heatCapacity = heat capacity of the material
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );

temperatureSolutions = zeros(size(coords, 2), timeSteps);
temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Backward Euler Time Integration
  for t = 2:timeSteps    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Prepare the Analysis
    activeCoords = getActiveCoordinates(coords, t, timeSteps);
    activeTemperatureSolution = getActiveTemperatureSolution(temperatureSolutions, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonNonLinearProblemTransient(activeCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    
    %Solve the problem using Newton-Raphson scheme
    activeTemperatureSolution = newtonRaphsonGlobalProblem(activeTemperatureSolution, poissonTransientProblem, tolerance, maxIterations, timeStepSize);

    %Post-Processing
    mergedTemperature = mergeActiveSolutionInGlobalDomain( activeTemperatureSolution, size(temperatureSolutions,1));
    temperatureSolutions(:, t) = mergedTemperature;

    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
%     internalEnergy(t) = temperatureSolution'*K*temperatureSolution;
    
  end
end

