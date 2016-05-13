function [temperatureSolutions, heatFluxes, refinedTemperatureSolutions, refinedHeatFluxes]= backwardEulerRefined(coords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refineDepth)
% BACKWARDEULERREFINED computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );

temperatureSolutions = zeros(size(coords, 2), timeSteps);
refinedTemperatureSolutions = zeros(2^refineDepth+1, timeSteps);

heatFluxes = zeros(size(coords, 2), timeSteps);
refinedHeatFluxes = zeros(2^refineDepth+1, timeSteps);


formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

  for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Prepare the Analysis
    activeCoords = getActiveCoordinates(coords, t, timeSteps);
    activeTemperatureSolution = getActiveTemperatureSolution(temperatureSolutions, t, timeSteps);
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(activeCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assemblyAndApplyStrongBCs(poissonTransientProblem);
    
    %Backward Euler Scheme
    RHS = timeStepSize * (f - K * activeTemperatureSolution);
    LHS = M + timeStepSize * K;
    
    temperatureIncrement = LHS\RHS;
    
    %Update temperature and post-process
    mergedTemperature = mergeActiveSolutionInGlobalDomain( activeTemperatureSolution + temperatureIncrement, size(temperatureSolutions,1));
    temperatureSolutions(:, t) = mergedTemperature;
    mergedHeatFluxes = mergeActiveSolutionInGlobalDomain( -evaluateHeatFlux( poissonTransientProblem, temperatureSolutions, t ), size(heatFluxes,1));
    heatFluxes(:, t) = mergedHeatFluxes;
    
    %Generate the local problem
    poissonTransientProblemRefined = poissonProblemTransientRefined(poissonTransientProblem, (activeTemperatureSolution + temperatureIncrement), refineDepth);
    [M_ref, K_ref, f_ref] = assemblyAndApplyStrongBCs(poissonTransientProblemRefined);
    
    %Backward Euler Scheme
    RHS = timeStepSize * (f_ref - K_ref * refinedTemperatureSolutions(:,t-1));
    LHS = M_ref + timeStepSize * K_ref;
    
    temperatureIncrement_ref = LHS\RHS;
    
    %Update temperature and post-process
    refinedTemperatureSolutions(:, t) = refinedTemperatureSolutions(:, t-1) + temperatureIncrement_ref;
    refinedHeatFluxes(:, t) = -evaluateHeatFlux( poissonTransientProblemRefined, refinedTemperatureSolutions, t );
    
  end
  
end






