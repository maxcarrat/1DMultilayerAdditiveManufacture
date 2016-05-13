function [temperatureSolutions, heatFluxes, internalEnergy]= backwardEulerSolver(coords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector)
% BackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );
temperatureSolutions = zeros(size(coords, 2), timeSteps);
heatFluxes = zeros(size(coords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);


formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

  for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    poissonTransientProblem = poissonProblemTransient(coords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assemblyAndApplyStrongBCs(poissonTransientProblem);
    
    RHS = timeStepSize*(f - K*temperatureSolutions(:,t-1));
    LHS = M + timeStepSize*K;
    
    temperatureIncrement = LHS\RHS;
    
  % update temperature
  
    temperatureSolutions(:, t) = temperatureSolutions(:,t-1) + temperatureIncrement;
    heatFluxes(:, t) = -evaluateHeatFlux( poissonTransientProblem, temperatureSolutions, t );
    
  end
  
end






