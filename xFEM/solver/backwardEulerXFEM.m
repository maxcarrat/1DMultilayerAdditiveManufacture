function [temperaturePostProcessing, heatFluxes, internalEnergy, modes]= backwardEulerXFEM(coords, postProcessingCoords, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, timeVector, refinementDepth, numberOfTrainingTimeSteps, maxIterations, tolerance)
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

refinedTemperatureSolutions = zeros(2^refinementDepth, 1);
localRefinedTemperatureSolutions = zeros(2^refinementDepth, timeSteps);

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
    poissonTransientProblem = poissonProblemTransient(refinedMesh, rhs, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
    [M, K, f] = assembly(poissonTransientProblem);
    
    %Project old solution onto the new mesh
    previousSolution = refinedTemperatureSolutions;
    if (norm(refinedTemperatureSolutions))~=0.0
        refinedTemperatureSolutions = L2projection(poissonTransientProblem, previousSolution, refinedMesh, previousMesh);
    end
    
    %Backward Euler Scheme
    [LHS, RHS] = applyBCs(M, K, f, poissonTransientProblem, refinedTemperatureSolutions, timeStepSize);
    temperatureIncrement = LHS\RHS;
    
    %Update and merge temperature into global domain
    refinedTemperatureSolutions = refinedTemperatureSolutions + temperatureIncrement;
    localRefinedTemperatureSolutions(:, t) = refinedTemperatureSolutions(t-1:end);    
    previousMesh = refinedMesh;

    %Post-Processing
    temperaturePostProcessing(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, refinedTemperatureSolutions, 0) ;
    heatFluxes(:, t) = evaluateNumericalResults(postProcessingCoords, poissonTransientProblem, refinedTemperatureSolutions, 1);
    internalEnergy(t) = refinedTemperatureSolutions'*K*refinedTemperatureSolutions;
  end
  
 %% Generate the reduced basis
 
  [solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,2:numberOfTrainingTimeSteps));
  
 %% Enriched mesh using RB

 numberOfModesSupports=1;

  for t = (numberOfTrainingTimeSteps+1):timeSteps
    
    formatSpec = 'Backward Euler Time Step (Reduced Basis): %1.1f \n' ;
    fprintf(formatSpec,t-1)
    
    currentTime = timeStepSize * (t-1);
    
    %Get active mesh
    activeMesh = getActiveCoordinates(coords, t, timeSteps);
    
    %Project old solution onto the new mesh
    temperatureSolutionsProjected = L2projectionEnriched(refinedTemperatureSolutions, activeMesh, refinedMesh, modes);
   
   %Solve Local problem iteratively
    for i=1:maxIterations
        
        %Solve Global/Coarse problem
        poissonTransientProblem = poissonProblemTransient(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
            rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime);
        temperatureSolutionsGlobal = solveGlobalProblem(temperatureSolutionsProjected, poissonTransientProblem, timeStepSize);
        
        %Update solutions
        temperatureSolutions = temperatureSolutionsGlobal;
        
        %Solve Local problem
        poissonTransientProblemEnriched = poissonProblemTransientEnriched(activeMesh, rhs, leftDirichletBoundaryConditionValue,...
            rightDirichletBoundaryConditionValue, k, heatCapacity, currentTime, refinementDepth, solutionReductionOperator);
        temperatureSolutionsLocal = solveLocalProblem(temperatureSolutionsProjected, poissonTransientProblemEnriched, timeStepSize, numberOfModesSupports);
      
        %Update solutions
        temperatureSolutions(end-numberOfModesSupports:end) = temperatureSolutionsLocal(end-numberOfModesSupports:end);
        temperatureSolutionsProjected(1:end-numberOfModesSupports) = temperatureSolutions(1:end-numberOfModesSupports);
        
        %Check convergence
        residualSolution = norm(temperatureSolutions(end-numberOfModesSupports:end)-temperatureSolutionsLocal(end-numberOfModesSupports:end));
        if  residualSolution<= tolerance
            break;
        end
        
    end

%    temperatureSolutionsEnriched = reducedBasisCoefficients;
%     
%     %Update and merge temperature into global domain
% %     temperatureSolutions(end-numberOfModesSupports:end) = temperatureSolutions(end-numberOfModesSupports:end) + temperatureSolutionsEnriched;
%     temperatureSolutions(end-numberOfModesSupports:end) = temperatureSolutionsEnriched;

    refinedMesh = activeMesh;
    refinedTemperatureSolutions = temperatureSolutions;

    temperaturePostProcessing(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh,...
        poissonTransientProblemEnriched, temperatureSolutions, 0) ;
    heatFluxes(:, t) = evaluateNumericalResultsEnriched(postProcessingCoords, activeMesh, poissonTransientProblemEnriched,...
        temperatureSolutions, 1);
    
    [~, K, ~] = assembly(poissonTransientProblem);
    internalEnergy(t) = temperatureSolutions'*K*temperatureSolutions;
    
  end
  
end






