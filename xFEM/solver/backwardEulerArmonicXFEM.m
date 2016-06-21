function [temperaturePostProcessing, heatFluxes]= backwardEulerArmonicXFEM(coords, postProcessingCoords, rhs,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, numberOfModes)
% BackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem. 
% Moreover, the numerical solution for each element is also computed
%   coords = coordinates of the mesh points
%   problem = struct that defines the boundary value problem
%   timeVector = vector of time steps for Backward Euler implicit scheme

timeSteps=size(timeVector,2);
timeStepSize=max(timeVector)/( timeSteps );

temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);
heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);


formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

    poissonTransientProblem = poissonProblemArmonicXFEM(coords, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, 0, numberOfModes);
    temperatureSolutions = zeros(poissonTransientProblem.gdof, 1);


for t = 2:timeSteps
    
    formatSpec = 'Backward Euler Time Step: %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    poissonTransientProblem = poissonProblemArmonicXFEM(coords, rhs, leftDirichletBoundaryConditionValue,...
        rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, numberOfModes);

    [M, K, f] = assemblyXFEM(poissonTransientProblem);
    [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, temperatureSolutions(:), timeStepSize);
    
    temperatureIncrement = LHS\RHS;
    
    % update temperature
    
    temperatureSolutions(:) = temperatureSolutions(:) + temperatureIncrement;
    
    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateXFEMNumericalResults(postProcessingCoords, poissonTransientProblem, temperatureSolutions(:), 0) ;
    heatFluxes(:, t) = evaluateXFEMNumericalResults(postProcessingCoords, poissonTransientProblem, temperatureSolutions(:), 1);
end

end






