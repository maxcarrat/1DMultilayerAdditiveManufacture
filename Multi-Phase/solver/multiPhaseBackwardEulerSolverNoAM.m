function [ temperaturePostProcessing, heatFluxes, analyticalTemperaturePostProcessing, computationalTime, DOFs ]...
    = multiPhaseBackwardEulerSolverNoAM( coords, postProcessingCoords, rhs, T0,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, neumannBoundaryconditionValue,...
    k, heatCapacity, timeVector, tolerance,...
    refinementDepth, integrationOrder, maxIterations )
%MULTIPHASEBACKWARDEULERSOLVERNOAM Summary of this function goes here
%   Detailed explanation goes here

timeSteps=size(timeVector,2);

timeStepSize=max(timeVector)/( timeSteps );
temperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);
analyticalTemperaturePostProcessing = zeros(size(postProcessingCoords, 2), timeSteps);

refinedTemperatureSolutions = linspace(T0,...
    T0, (numel(coords)-1) * 2^refinementDepth)';

convergence = 0.0;
computationalTime = [];

heatFluxes = zeros(size(postProcessingCoords, 2), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%Refinement
refinedMesh = linspace(coords(1),coords(end),(numel(coords)-1) * 2^refinementDepth);

timeToGenerateAndSolveTheSystem = 0.0;

for t = 1:timeSteps
    tic
    
    formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
    fprintf(formatSpec, t)
    
    currentTime = timeStepSize * t;
    
    %Generate the Poisson problem at timeStep t
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs,...
        leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
        neumannBoundaryconditionValue, k, heatCapacity, currentTime);
    
    %Update and merge temperature into global domain
    [refinedTemperatureSolutions, convergenceFlag] = solveMultiPhaseProblem( poissonTransientProblem,...
        currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
    
    convergence = convergence + convergenceFlag;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    
    %Post-Processing
    temperaturePostProcessing(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, refinedTemperatureSolutions, 0) ;
    heatFluxes(:, t+1) = evaluateNumericalResults(postProcessingCoords, currentTime, poissonTransientProblem, refinedTemperatureSolutions, 1);
    analyticalTemperaturePostProcessing(:, t+1) = analyticalTemperature(postProcessingCoords, currentTime);
    
    DOFs = numel(refinedTemperatureSolutions);
    
end

computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];

if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
end
    
end


function solution = analyticalTemperature(x_PostProcessing, currentTime)

solution = zeros(size(x_PostProcessing, 2), 1);

lambda_anal = 0.246395135644045;
alpha = 6.8224e-06;
X_meltingFront = 2*lambda_anal*sqrt(alpha*currentTime);

for i=1:numel(x_PostProcessing)
    if x_PostProcessing(i) <= X_meltingFront
        solution(i) = 2000 - (330)*(erf(x_PostProcessing(i)/ 2 / sqrt(alpha*currentTime))/erf(lambda_anal));
    else
        solution(i) = 1000 + (670)*(erfc(x_PostProcessing(i)/ 2 / sqrt(alpha*currentTime))/erfc(lambda_anal));        
    end
end
end

