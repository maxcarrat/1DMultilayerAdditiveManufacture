function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]=...
    XIGAMultiPhaseBackwardEulerSolver( p, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance,...
    maxIterations, numberOfRefinedElementsToBeKept,...
    TotalNumberOfControlPoints, refinementDepth, numberOfEnrichedControlPoints,...
    numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes,...
    integrationOrder, integrationModesOrder, layerThickness)
% IGAMultiPhaseBackwardEulerSolver computes the 1D h-FEM numerical solution of a boundary value problem.
% Moreover, the numerical solution for each element is also computed
%   CPs = control points
%   rhs = external heat source
%   leftDirichletBoundaryConditionValue
%   rightDirichletBoundaryConditionValue
%   k = conductivity of the material
%   heatCapacity = heat capacity of the material
%   timeVector = vector of time steps for Backward Euler implicit scheme
%   refinementDepth = depth of the refined mesh
%   numberOfTrainingTimeSteps = number of time steps in the training phase

%initialize variables
timeSteps=size(timeVector,2);

timeStepSize=max(timeVector)/( timeSteps );
temperaturePostProcessing = zeros(length(postProcessingCoords), timeSteps);

convergence = 0.0;
computationalTime = [];

heatFluxes = zeros(length(postProcessingCoords), timeSteps);
internalEnergy = zeros(timeSteps, 1);

formatSpec = 'Begin Time Integration Scheme \n' ;
fprintf(formatSpec)

%% Training phase
for layer = 1:numberOfTrainingLayers
    
    timeToGenerateAndSolveTheSystem = 0.0;
    tic 
    
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %% Refinement
    % In order to get POD modes corresponding to the last layer the
    % interface between the last ( refined ) layer and the other has to be
    % C0 continuous, i.e. we have to insert an additional knot at the layer
    % interface.

    newKnots = linspace( knotVector(end-(p+1)),  knotVector(end-p), 2^refinementDepth + p);
   
    if layer == 1 || p == 1
        [ CPs, refinedKnotVector] = refineKnotVector( p, CPs, knotVector, newKnots(2:end-1));
    else
        [ CPs, refinedKnotVector] = refineKnotVector( p, CPs, knotVector, newKnots(1:end-1));
    end

    %Project old solution onto the new mesh
    if layer > 1
        
        % create a fictitious problem for the projection
        poissonProblem = poissonProblemTransientIGA(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, 0.0,...
            refinedKnotVector, p, refinementDepth);
        
        %L2 projection onto new IGA mesh
        refinedTemperatureSolutions = L2projectionIGA( refinedTemperatureSolutions, poissonProblem, integrationOrder,...
            0, layerThickness, initialTemperature, 'false' );
    else
        refinedTemperatureSolutions = zeros(length(CPs), 1);
        localRefinedTemperatureSolutions = zeros(length(CPs), timeSteps);
    end
    
    % foor loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        poissonProblem = poissonProblemTransientIGA(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            refinedKnotVector, p, refinementDepth);
        
        %Update and merge temperature into global domain
        [refinedTemperatureSolutions, convergenceFlag] = solveMultiPhaseProblemIGA( poissonProblem,...
            currentTime, timeStepSize, integrationOrder, tolerance, maxIterations, refinedTemperatureSolutions );
        mergedTemperature = mergeActiveSolutionInGlobalDomain(refinedTemperatureSolutions, TotalNumberOfControlPoints);
        
        %Cache the local solution for POD
        localRefinedTemperatureSolutions(:, t+1) = getIGASolution(refinedTemperatureSolutions,...
            layer, poissonProblem, numberOfRefinedElementsToBeKept);
        
        %Cache knot vector
        previousKnotVector = poissonProblem.knotVector;
        
        %Update convergence flag and register time to solve the time step
        convergence = convergence + convergenceFlag;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 1);
        
        %Register number of Dofs
        DOFs = numel(refinedTemperatureSolutions);
        
    end
 
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];

end

%% Generate the reduced basis
% POD on the local/layer solution snapshots, omitt the first zero-solution
% vector
[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,3*numberOfLayersTimeSteps+2:numberOfTrainingLayers*numberOfLayersTimeSteps+1), numberOfPODModes);

%% Project refined solution onto enriched mesh

%next layer
layer = (numberOfTrainingLayers+1);

%Generate a coarse IGA mesh
knotVector = getOpenKnotVector( layer, p );
CPs = getControlPoints( layer, layerThickness, p );

%Enrich only the bar tip CP
activeNumberOfCPs = numberOfEnrichedControlPoints;

%Generate the XIGA Poisson problem at layer for projection
poissonXProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
    knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);

%Project global solution onto the enriched modal space
temperatureSolutions = L2projectionIGA( refinedTemperatureSolutions, poissonXProblem, integrationOrder,...
            integrationModesOrder, layerThickness, initialTemperature, 'transition', poissonProblem );

%% Enriched mesh using RB

% for loop ROM phase layers
for layer = (numberOfTrainingLayers+1):numberOfLayers
    timeToGenerateAndSolveTheSystem = 0.0;
    tic
    
    if layer > (numberOfTrainingLayers+1)
        %Generate a coarse IGA mesh
        knotVector = getOpenKnotVector( layer, p );
        CPs = getControlPoints( layer, layerThickness, p );
        
        %Enrich only the bar tip CP
        activeNumberOfCPs = numberOfEnrichedControlPoints;
        
        %Generate the XIGA Poisson problem at layer for projection
        poissonXProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        
        % L2 projection
        temperatureSolutions = L2projectionIGA( refinedTemperatureSolutions, poissonXProblem, integrationOrder,...
            integrationModesOrder, layerThickness, initialTemperature, 'true', poissonProblem );
        
    end
    

    % for loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate Local problem
        poissonXProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        [temperatureSolutions, convergenceFlag] = solveXIGAMultiPhase( poissonXProblem, currentTime,...
            timeStepSize,integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
        
        %Check convergence
        convergence = convergence + convergenceFlag;
        
        %Post-Processing temperatures and heat fluxes
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
        poissonXProblem, temperatureSolutions,...
            layer, numberOfLayers, 0);
        heatFluxes(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
            poissonXProblem,temperatureSolutions,...
            layer, numberOfLayers, 1);
        
        %register dofs
        DOFs = numel(temperatureSolutions);
        
    end
    
    poissonProblem = poissonXProblem;
    refinedTemperatureSolutions = temperatureSolutions;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;   

    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
    
end

if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
end
end
