function [temperatureSolutionsPostProcessed, heatFluxesPostProcessed, internalEnergy]= backwardEulerMultilayers(coords, postProcessingCoords, rhs, initialTemperature,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector)
% BackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );
temperatureSolutions = zeros(size(coords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

temperatureSolutionsPostProcessed = zeros(size(postProcessingCoords, 2), timeSteps);
heatFluxesPostProcessed = zeros(size(postProcessingCoords, 2), timeSteps);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

  for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Prepare the Analysis
    activeCoords = getActiveCoordinates(coords, t, timeSteps);
    activeTemperatureSolution = getActiveTemperatureSolution(temperatureSolutions, initialTemperature, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(activeCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assembly(poissonTransientProblem);
    
    %Backward Euler Scheme
    [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, activeTemperatureSolution, timeStepSize);
    
    temperatureIncrement = LHS\RHS;
        
    %Update temperature and post-process
    activeTemperatureSolution = activeTemperatureSolution + temperatureIncrement;
    mergedTemperature = mergeActiveSolutionInGlobalDomain( activeTemperatureSolution, size(temperatureSolutions,1));
    temperatureSolutionsPostProcessed(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 0) ;
    heatFluxesPostProcessed(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, mergedTemperature, 1);
    internalEnergy(t) = activeTemperatureSolution'*K*activeTemperatureSolution;

    temperatureSolutions(:, t) = mergedTemperature;
  end
  
end






