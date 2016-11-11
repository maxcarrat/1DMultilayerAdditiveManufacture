function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]=...
    multiscalePODhIGASolver( p, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance,...
    maxIterations, totalNumberOfControlPoints, refinementDepth, numberOfEnrichedRefinementDepth,...
    numberOfTrainingLayers, numberOfLayersTimeSteps, numberOfLayers, numberOfPODModes,...
    integrationSplinesOrder, integrationModesOrder, layerThickness, steelThermalConductivityDerivative, heatCapacityDerivative)
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
%steelThermalConductivityDerivative = derivative of the conductivity of the
%material
%heatCapacityDerivative = derivative of the heat capacity of the
%material
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
for layer = 2:numberOfTrainingLayers
    
    timeToGenerateAndSolveTheSystem = 0.0;
    tic
    if layer==1
        ansatzOrder=p;
    else
        ansatzOrder=p;
    end
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, ansatzOrder );
    CPs = getControlPoints( layer, layerThickness, ansatzOrder );
    
    %Generate a fine overlay FEM mesh
    nodalCoordinates = linspace( (layer-1) * layerThickness,  CPs(end), 2^refinementDepth + 1);
    
    %Project old solution onto the new mesh
    if layer > 2
        %Create a fictitious problem for the projection
        %Base mesh problem
        baseProblem = poissonProblemBaseMesh(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, 0.0,...
            knotVector, p, refinementDepth);
        
        %Overlay mesh problem
        overlayProblem = poissonProblemOverlayMesh(nodalCoordinates,...
            numberOfEnrichedRefinementDepth, rhs,...
            neumannBoundaryconditionValue, p, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, 0.0);
        
        %L2 projection onto new IGA mesh
        [baseTemperatureSolutions, overlayTemperatureSolutions] = L2MultiscaleProjection...
            ( baseTemperatureSolutions, overlayTemperatureSolutions, baseProblem,...
             overlayProblem, integrationSplinesOrder, ansatzOrder + 1, layerThickness, initialTemperature,...
            'false', previousBaseProblem, previousOverlayProblem, layer );
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
        baseProblem = poissonProblemBaseMesh(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, currentTime,...
            knotVector, ansatzOrder, refinementDepth);
        
        %Overlay mesh problem
        overlayProblem = poissonProblemOverlayMesh(nodalCoordinates, numberOfEnrichedRefinementDepth, rhs,...
            neumannBoundaryconditionValue, ansatzOrder, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, currentTime);
        
        %Solve multiscale problem using one-step Gauss-Seidel-Newton method
        %set the eXtended index = 'false' since we are in the training
        %phase!
        overlayIntegrationOrder = integrationSplinesOrder - 1;
        [baseTemperatureSolutions, overlayTemperatureSolutions, convergenceFlag] = solveOneStepGaussSeidelNewton( ...
            overlayProblem, overlayProblem, baseProblem, currentTime, timeStepSize, integrationSplinesOrder, overlayIntegrationOrder,...
            tolerance, maxIterations, baseTemperatureSolutions, overlayTemperatureSolutions, initialTemperature, 'false' );
        
        %Cache the local solution for POD
        %in order to get the correct POD modes we need to decompose the
        %temperature evaluated at the nodes of the overlay mesh and not
        %only the overlay mesh nodal values!

%         temperatureAtOverlayNodes = evaluateNumericalResultsMultiscale(nodalCoordinates,...
%             currentTime, baseProblem, overlayProblem, baseTemperatureSolutions,...
%             overlayTemperatureSolutions, 1, numberOfLayers, 0);
%        localRefinedTemperatureSolutions(:, t+1) = temperatureAtOverlayNodes;

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
    previousOverlayProblem = overlayProblem;
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
    
end

%% Generate the reduced basis ---------------------------------------------
% POD on the local/layer solution snapshots, omitt the first zero-solution
% vector
[solutionReductionOperator, ~] = properOrthogonalDecomposition...
    (localRefinedTemperatureSolutions(:,3*numberOfLayersTimeSteps+2:numberOfTrainingLayers*numberOfLayersTimeSteps+1), numberOfPODModes);

% Project refined solution onto enriched mesh

%next layer
layer = (numberOfTrainingLayers+1);

%Generate a coarse IGA mesh
knotVector = getOpenKnotVector( layer, p );
CPs = getControlPoints( layer, layerThickness, p );

%Generate a one-element FEM mesh
nodalCoordinates = linspace( (layer-1) * layerThickness,  CPs(end), 2^(numberOfEnrichedRefinementDepth)+1);

%Generate the Poisson problem at timeStep t
%Cache previous problem
previousBaseProblem = baseProblem;

% %Base mesh problem
baseProblem = poissonProblemBaseMesh(CPs, rhs,...
    leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, steelThermalConductivityDerivative,...
    heatCapacity, heatCapacityDerivative, currentTime,...
    knotVector, p, refinementDepth);

%eXtended overlay mesh problem
overlayXProblem = poissonProblemXtendedOverlayMesh(nodalCoordinates,...
    numberOfEnrichedRefinementDepth, rhs, neumannBoundaryconditionValue, p, k,...
    steelThermalConductivityDerivative, heatCapacity, heatCapacityDerivative,...
    solutionReductionOperator, currentTime);

%         problemXtended = poissonProblemTransientXIGA(CPs, numberOfEnrichedControlPoints, rhs,...
%             leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
%             neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
%             knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        

% %Project global solution onto the enriched modal space
[ baseTemperatureSolutions, overlayTemperatureSolutions ] = L2MultiscaleProjection(...
    baseTemperatureSolutions, overlayTemperatureSolutions,...
    baseProblem, overlayXProblem, integrationSplinesOrder, integrationModesOrder,...
    layerThickness, initialTemperature, 'transition', previousBaseProblem,...
    overlayProblem, layer);

% temperatureSolutions = L2MultiscaleOntoEnrichedProjection(...
%     baseTemperatureSolutions, problemXtended, integrationSplinesOrder, integrationModesOrder,...
%     layerThickness, initialTemperature);

% %% TEST
% %Post-Processing temperatures and heat fluxes
% temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords,...
%     currentTime, baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions,...
%     overlayTemperatureSolutions, layer, numberOfLayers, 0) ;
% heatFluxes(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords, currentTime,...
%     baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions, overlayTemperatureSolutions,...
%     layer, numberOfLayers, 1);
% %% END TEST

%% ROM phase --------------------------------------------------------------

% for loop ROM phase layers
for layer = (numberOfTrainingLayers+1):numberOfLayers
    timeToGenerateAndSolveTheSystem = 0.0;
    tic
    
    if layer > (numberOfTrainingLayers+1)
        %Generate a coarse IGA mesh
        knotVector = getOpenKnotVector( layer, p );
        CPs = getControlPoints( layer, layerThickness, p );
        
        %Generate a one-element FEM mesh
        nodalCoordinates = linspace( (layer-1) * layerThickness,  CPs(end), 2^(numberOfEnrichedRefinementDepth)+1);
        
        %Base mesh problem
        baseProblem = poissonProblemBaseMesh(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, currentTime,...
            knotVector, p, refinementDepth);
        
        %eXtended overlay mesh problem
        overlayXProblem = poissonProblemXtendedOverlayMesh(nodalCoordinates,...
            numberOfEnrichedRefinementDepth, rhs, neumannBoundaryconditionValue, p, k,...
            steelThermalConductivityDerivative, heatCapacity, heatCapacityDerivative,...
            solutionReductionOperator, currentTime);

%         problemXtended = poissonProblemTransientXIGA(CPs, numberOfEnrichedControlPoints, rhs,...
%             leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
%             neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
%             knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        
        
        % L2 projection
                [ baseTemperatureSolutions, overlayTemperatureSolutions ]  = L2MultiscaleProjection(...
                    baseTemperatureSolutions, overlayTemperatureSolutions,...
                    baseProblem, overlayXProblem, integrationSplinesOrder, integrationModesOrder,...
                    layerThickness, initialTemperature, 'reduced', previousBaseProblem,...
                    overlayProblem, layer);
        
%         temperatureSolutions = L2projectionIGA( temperatureSolutions, problemXtended, integrationSplinesOrder,...
%             integrationModesOrder, layerThickness, initialTemperature, 'true', previousXtendedProblem );
        
        
        %         %% TEST
        %         %Post-Processing temperatures and heat fluxes
        %         temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords,...
        %             currentTime, baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions,...
        %             overlayTemperatureSolutions, layer, numberOfLayers, 0) ;
        %         heatFluxes(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords, currentTime,...
        %             baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions, overlayTemperatureSolutions,...
        %             layer, numberOfLayers, 1);
        %         %% END TEST

    end
    
    
    % for loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Base mesh problem
        baseProblem = poissonProblemBaseMesh(CPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, steelThermalConductivityDerivative,...
            heatCapacity, heatCapacityDerivative, currentTime,...
            knotVector, p, refinementDepth);
        
        %eXtended overlay mesh problem
        overlayXProblem = poissonProblemXtendedOverlayMesh(nodalCoordinates,...
            numberOfEnrichedRefinementDepth, rhs, neumannBoundaryconditionValue, p, k,...
            steelThermalConductivityDerivative, heatCapacity, heatCapacityDerivative,...
            solutionReductionOperator, currentTime);
        

%         problemXtended = poissonProblemTransientXIGA(CPs, numberOfEnrichedControlPoints, rhs,...
%             leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
%             neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
%             knotVector, p, refinementDepth, numberOfEnrichedControlPoints, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve using one-step Gauss-Seidel-Newton solver
        %%N.B. set the eXtended flag to true
        [baseTemperatureSolutions, overlayTemperatureSolutions, convergenceFlag] = solveOneStepGaussSeidelNewton( ...
            overlayProblem, overlayXProblem, baseProblem, currentTime, timeStepSize, integrationSplinesOrder, integrationModesOrder,...
            tolerance, maxIterations, baseTemperatureSolutions, overlayTemperatureSolutions, initialTemperature, 'true' );
  
%         [temperatureSolutions, convergenceFlag] = solveXIGAMultiPhase( problemXtended, currentTime,...
%             timeStepSize, integrationSplinesOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
        
        %Check convergence
        convergence = convergence + convergenceFlag;
        
        %Post-Processing temperatures and heat fluxes
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords,...
            currentTime, baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions,...
            overlayTemperatureSolutions, layer, numberOfLayers, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResultsXtendedMultiscale(postProcessingCoords, currentTime,...
            baseProblem, overlayProblem, overlayXProblem, baseTemperatureSolutions, overlayTemperatureSolutions,...
            layer, numberOfLayers, 1);
       
%         modes = problemXtended.modes;
%         
%         temperaturePostProcessing(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
%         problemXtended, temperatureSolutions, layer, numberOfLayers, 0);
%         heatFluxes(:, t+1) = evaluateNumericalResultsXIGA(postProcessingCoords, currentTime, modes,...
%             problemXtended,temperatureSolutions, layer, numberOfLayers, 1);
        
        %register dofs
        DOFs = numel(baseTemperatureSolutions)+numel(overlayTemperatureSolutions);
        
    end
    
    timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
%     previousXtendedProblem = problemXtended;
    previousBaseProblem = baseProblem;
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
    
end

if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
end
end

