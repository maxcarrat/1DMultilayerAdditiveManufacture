function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]=...
    multiscalePODhIGASolver( p, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance,...
    maxIterations, totalNumberOfControlPoints, refinementDepth, numberOfEnrichedControlPoints,...
    numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes,...
    integrationSplinesOrder, integrationModesOrder, layerThickness)
%%MULTISCALEPODHIGASOLVER the function solve a AM process of a 1D bar using
%multiscale POD h-IGA method. The process is devided in two sperate phases:
%1. Training phase: a base IGA mesh and an overlay h-FEM mesh are
%solved using multilevel hp-d approach (Rank, E. (1992b). A combination of 
%hp-version finite elements and a domain decomposition method. Proceedings
%of the First European Conference on Numerical Methods in Engineering). The
%overlay h-FEM mesh is created at each time step at the last layer of the
%bar and compries 2^refinementDepth linear elements.
%2. ROM phase: the base IGA mesh does not change but now the overlay mesh
%is a coarse(one single element) X-FEM mesh enriched by means of POD modes.
%The POD modes are obtained applying the POD to the set of solution
%vectores cahced at each time step during the training phase.
%
%Input:
%p = order of the BSplines
%postProcessingCoords = coordinates of post-processing grid
%rhs = external heat source,
%initialTemperature = initial temperature of the powder and of the basement
%leftDirichletBoundaryConditionValue = BCs on the bottom of the bar
%rightDirichletBoundaryConditionValue = BCs at the top of the bar
%neumannBoundaryconditionValue = external heat fluxes
%k = conductivity of the material
%heatCapacity = heat capacity of the material
%timeVector = time discretization vector
%tolerance = converegence tolerance for Gauss-Seidel-Newton nonlinear
%solver
%maxIterations = max iteration number in Gauss-Seidel_Newton solver
%totalNumberOfControlPoints = number of control pioints on the initial 
%coarse mesh
%refinementDepth = refinement level depth at the last layer
%numberOfEnrichedControlPoints = enriched control points on the last layer
%numberOfTrainingLayer = number of layer for the training phase
%numberOfLayer = total number of layers in the process
%numberOfPODModes = number of POD modes to be considered
%integrationSplinesOrder = number of Gauss points on the base mesh
%integrationModesOrder = number of Gauss points on the overlay eXtended mesh
%layerThickness = thickness of the single layer
%
%Output:
%temperaturePostProcessing = temperature solution at the post-processing
%grid nodes
%heatFluxes = heat fluxes solution at the post-processing grid nodes
%internalEnergy = internal energy of the system
%computationalTime = computational time of each time step of the process
%DOFs = degrees of freedom

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

%% Training phase ---------------------------------------------------------
for layer = 1:numberOfTrainingLayers
    
    timeToGenerateAndSolveTheSystem = 0.0;
    tic 
    
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %Generate a fine overlay FEM mesh
    nodalCoordinates = linspace( knotVector(end-(p+1)),  knotVector(end-p), 2^refinementDepth);
    
    %Project old solution onto the new mesh
    if layer > 1
        %Create a fictitious problem for the projection
            %Base mesh problem
            baseProblem = poissonProblemTransientIGA(CPs, rhs,...
                leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
                neumannBoundaryconditionValue, k, heatCapacity, 0.0,...
                knotVector, p, refinementDepth);

            %Overlay mesh problem
            overlayProblem = poissonProblemOverlayMesh(nodalCoordinates,...
                numberOfEnrichedControlPoints, rhs,...
                leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
                neumannBoundaryconditionValue, p, k, heatCapacity, 0.0);   
        
        %L2 projection onto new IGA mesh
        [baseTemperatureSolutions, overlayTemperatureSolutions] = L2MultiscaleProjection...
            ( baseTemperatureSolutions, overlayTemperatureSolutions,baseProblem,...
            overlayProblem, integrationSplinesOrder, 0, initialTemperature, 'false', previousBaseProblem );
    else
        baseTemperatureSolutions = zeros(length(CPs), 1);
        overlayTemperatureSolutions = zeros(length(nodalCoordinates), 1);
        localRefinedTemperatureSolutions = zeros(length(nodalCoordinates), timeSteps);
    end
    
    % foor loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(Training Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate the Poisson problem at timeStep t
        %Base mesh problem
        baseProblem = poissonProblemTransientIGA(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            knotVector, p, refinementDepth);
        
        %Overlay meah problem
        overlayProblem = poissonProblemOverlayMesh(nodalCoordinates, numberOfEnrichedControlPoints, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, p, k, heatCapacity, currentTime);
        
        %Update and merge temperature into global domain
        overlayIntegrationOrder = 2;
        [baseTemperatureSolutions, overlayTemperatureSolutions, convergenceFlag] = solveOneStepGaussSeidelNewton( baseProblem,...
            overlayProblem, currentTime, timeStepSize, integrationSplinesOrder, overlayIntegrationOrder,...
            tolerance, maxIterations, baseTemperatureSolutions, overlayTemperatureSolutions );
                
        %Cache the local solution for POD
        localRefinedTemperatureSolutions(:, t+1) = overlayTemperatureSolutions;
        
        %Update convergence flag and register time to solve the time step
        convergence = convergence + convergenceFlag;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsMultiscale(postProcessingCoords,...
            currentTime, baseProblem, overlayProblem, baseTemperatureSolutions,...
            overlayTemperatureSolutions, layer, numberOfLayers, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResultsMultiscale(postProcessingCoords, currentTime,...
            baseProblem, overlayProblem, baseTemperatureSolutions, overlayTemperatureSolutions,...
            layer, numberOfLayers, 1);
        
        %Register number of Dofs
        DOFs = numel(baseTemperatureSolutions) + numel(overlayTemperatureSolutions);
        
    end
 
    previousBaseProblem = baseProblem;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];

end

%% Generate the reduced basis ---------------------------------------------
% POD on the local/layer solution snapshots, omitt the first zero-solution
% vector
[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,3*numberOfLayersTimeSteps+2:numberOfTrainingLayers*numberOfLayersTimeSteps+1), numberOfPODModes);

% Project refined solution onto enriched mesh

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

%% ROM phase --------------------------------------------------------------

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

