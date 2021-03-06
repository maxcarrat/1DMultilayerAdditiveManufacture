function [temperaturePostProcessing, heatFluxes, internalEnergy, computationalTime, DOFs]=...
    IGAMultiPhaseBackwardEulerSolver( p, postProcessingCoords, rhs, ...
    initialTemperature, leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
    neumannBoundaryconditionValue, k, heatCapacity, timeVector, tolerance,...
    maxIterations, numberOfRefinedElementsToBeKept,...
    TotalNumberOfControlPoints, refinementDepth, PODRefinementDepth,...
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
    
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %Refinement
    newKnots = linspace( knotVector(end-(p+1)),  knotVector(end-p), 2^refinementDepth + 1);
    [ CPs, refinedKnotVector] = refineKnotVector( p, CPs, knotVector, newKnots(2:end-1));
    
    %Project old solution onto the new mesh
    if layer > 1
        refinedTemperatureSolutions = L2projectionIGA(poissonProblem,...
            refinedTemperatureSolutions, CPs, previousKnotVector, initialTemperature);
    else
        refinedTemperatureSolutions = zeros(length(CPs), 1);
        localRefinedTemperatureSolutions = zeros(length(CPs), timeSteps);
    end
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    % foor loop time integration @layer
    for iTime = 1:numberOfLayersTimeSteps
        tic
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
            layer, numberOfLayers, length(CPs), refinementDepth, numberOfRefinedElementsToBeKept);
        
        %Cache knot vector
        previousKnotVector = poissonProblem.knotVector;
        
        %Update convergence flag and register time to solve the time step
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 0) ;
        heatFluxes(:, t+1) = evaluateNumericalResultsIGA(postProcessingCoords, currentTime, poissonProblem, mergedTemperature, layer, numberOfLayers, 1);
        
        %Register number of Dofs
        DOFs = numel(refinedTemperatureSolutions);
        
    end
    
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

%% Generate the reduced basis
% POD on the local/layer solution snapshots, omitt the first zero-solution
% vector
[solutionReductionOperator, modes] = properOrthogonalDecomposition(localRefinedTemperatureSolutions(:,1+numberOfLayersTimeSteps:numberOfTrainingLayers*numberOfLayersTimeSteps), numberOfPODModes);

%% Enriched mesh using RB
% for loop ROM phase layers
for layer = (numberOfTrainingLayers+1):numberOfLayers
    
    %Generate a coarse IGA mesh
    knotVector = getOpenKnotVector( layer, p );
    CPs = getControlPoints( layer, layerThickness, p );
    
    %Enrich only the bar tip CP
    activeNumberOfCPs = 2^PODRefinementDepth;
    
    %Generate the XIGA Poisson problem at layer for projection
    poissonProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
        leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
        neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
        knotVector, p, refinementDepth, PODRefinementDepth, solutionReductionOperator);
    
    %Project global solution onto the enriched modal space
    temperatureSolutions = projectOntoXIGAMesh(poissonProblem, refinedTemperatureSolutions,...
    modes, knotVector, previousKnotVector, PODRefinementDepth, initialTemperature);
    
    timeToGenerateAndSolveTheSystem = 0.0;
    
    % foor loop time integration @layer    
    for iTime = 1:numberOfLayersTimeSteps
        tic
        t = (layer - 1) * numberOfLayersTimeSteps + iTime;
        
        formatSpec = 'Backward Euler Time Step(ROM Phase): %1.1f \n' ;
        fprintf(formatSpec, t)
        
        currentTime = timeStepSize * t;
        
        %Generate Local problem
        poissonProblem = poissonProblemTransientXIGA(CPs, activeNumberOfCPs, rhs,...
            leftDirichletBoundaryConditionValue, rightDirichletBoundaryConditionValue,...
            neumannBoundaryconditionValue, k, heatCapacity, currentTime,...
            knotVector, p, refinementDepth, PODRefinementDepth, solutionReductionOperator);
        
        disp(' Solve Local Enriched Problem ');
        
        %Solve Local problem enriched
        [temperatureSolutions, convergenceFlag] = solveMultiPhaseProblemIGA( poissonProblem, currentTime,...
            timeStepSize,integrationOrder, integrationModesOrder, tolerance, maxIterations, temperatureSolutions );
        
        convergence = convergence + convergenceFlag;
        timeToGenerateAndSolveTheSystem = timeToGenerateAndSolveTheSystem + toc;
        
        %Post-Processing
        temperaturePostProcessing(:, t+1) = postProcessingProjection(postProcessingCoords, currentTime, poissonTransientProblemEnriched, temperatureSolutions,...
            modes, 0.0);
        heatFluxes(:, t+1) = postProcessingProjection(postProcessingCoords, currentTime, poissonTransientProblemEnriched,temperatureSolutions,...
            modes, 1);
        
        DOFs = numel(temperatureSolutions);
        
    end
    
    previousKnotVector = knotVector;
    computationalTime = [computationalTime, timeToGenerateAndSolveTheSystem];
end

if convergence < 1
    disp('The analysis always converged')
else
    disp('The analysis did not always converged !!!')
    fprintf(num2str(convergence));
end
end


function [ projectedCoefficients ] = postProcessingProjection(x, t, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTION project the previous solution onto the post-processing mesh.
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives


projectedCoefficients=zeros(size(x));
X1 = problem.coords(1);
X2 = problem.coords(2);

projectedCoefficients(x>=X1 & x<=X2) = postProcessingProjectionElement(1, x(x>=X1 & x<=X2), t, problem,... % first element
    solutionCoefficients, modes, derivative);

% loop over elements
for e=2:problem.N
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    projectedCoefficients(x>X1 & x<=X2) = postProcessingProjectionElement(e, x(x>X1 & x<=X2), t, problem,...
        solutionCoefficients, modes, derivative);
    
end

end

function [ projectedCoefficients ] = postProcessingProjectionElement(e, x, t, problem, solutionCoefficients,...
    modes, derivative)
% POSTPROCESSINGPROJECTIONELEMENT project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);

X1 = problem.coords(e);
X2 = problem.coords(e+1);

if e > problem.N - problem.XN  % element is active
    
    elementEnrichedIndex = e - (problem.N - problem.XN);
    
    if elementEnrichedIndex == 1
        indexLocalEnrichedNodes = 2; %lhs node
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    % On active elements use the refined domain as integration domain
    refinedNodes = 2^problem.refinementDepth+1;
    integrationDomain = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    subDomainShapeFunctionCoefficients = linspace(-1, +1, ceil(refinedNodes/problem.XN));
    
    Xi1 = integrationDomain(1);
    Xi2 = integrationDomain(2);
    
    localCoordinates = mapGlobalToLocal( x, X1, X2);
    projectedCoefficients=zeros(size(localCoordinates));
    
    projectedCoefficients(localCoordinates>=Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, 1, elementEnrichedIndex,...
        localCoordinates(localCoordinates>=Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
        solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
    
    for integrationSubDomainIndex=2:ceil(refinedNodes/problem.XN)-1
        
        Xi1 = integrationDomain(integrationSubDomainIndex);
        Xi2 = integrationDomain(integrationSubDomainIndex+1);
        
        projectedCoefficients(localCoordinates>Xi1 & localCoordinates<=Xi2) = postProcessingProjectionSubElements(x, integrationSubDomainIndex, elementEnrichedIndex,...
            localCoordinates(localCoordinates>Xi1 & localCoordinates<=Xi2), t, problem, integrationDomain,...
            solutionCoefficients, modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes);
        
        
    end
    
    projectedCoefficients = projectedCoefficients .* (2/(X2-X1)) ^ derivative;
    
    
    
else     % Element not enriched
    
    projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));
    
    N = zeros(length(x), 2);
    B = zeros(length(x), 2);
    
    localCoords = mapGlobalToLocal( x, X1, X2);
    
    for k=1:length(x)
        [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    end
    
    if derivative == 0
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
    else
        for i=1:length(x)
            projectionOperator(i,1:size(B,2)) = B(i,:);
        end
    end
    
    projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:)) .* (2/(X2-X1)) ^ derivative;
    
    if derivative == 1
        for i=1:length(x)
            projectionOperator(i,1:size(N,2)) = N(i,:);
        end
        temperature = projectionOperator * solutionCoefficients(problem.LM(e,:));
        for i=1:numel(projectedCoefficients)
            projectedCoefficients(i) = projectedCoefficients(i) * ...
                problem.k(x, t, temperature(i));
        end
        
    end
end


end


function [ projectedCoefficients ] = postProcessingProjectionSubElements(globalCoords, subDomainIndex, e, x, t, problem, integrationDomain, solutionCoefficients,...
    modes, derivative, subDomainShapeFunctionCoefficients, indexLocalEnrichedNodes)
% POSTPROCESSINGPROJECTIONSUBELEMENTS project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh in local coordinates of the integration domain
%   problem
%   integrationDomain = integration domain of the enriched element in local
%   coords
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, 2 + modes * length(indexLocalEnrichedNodes) );

localCoords = x;
refinedNodes = 2^problem.refinementDepth + 1;
N = zeros(length(x), 2);
B = zeros(length(x), 2);

F = zeros(length(x), modes *length(indexLocalEnrichedNodes));
G = zeros(length(x), modes *length(indexLocalEnrichedNodes));

PODCoefficients = problem.reductionOperator(...
    (e-1)*(floor(refinedNodes/problem.XN))+1:(e-1)*(floor(refinedNodes/problem.XN)) +...
    ceil(refinedNodes/problem.XN),:);

for k=1:length(x)
    
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
    
    [F(k,:), G(k,:)] = PODModesAndDerivativesGaussIntegration(problem, localCoords(k), modes, PODCoefficients,...
        subDomainShapeFunctionCoefficients, subDomainIndex, indexLocalEnrichedNodes);
end

if derivative == 0
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
else
    for i=1:length(localCoords)
        projectionOperator(i,1:size(B,2)+size(G,2)) = [B(i,:), G(i,:)];
    end
end

if length(indexLocalEnrichedNodes) == 1 %lhs
    coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
        solutionCoefficients(problem.N  + (e - 1) * length(indexLocalEnrichedNodes) * modes + 2 : problem.N...
        + (e) * length(indexLocalEnrichedNodes) * modes + 1 )];
else
    coefficients = [solutionCoefficients((problem.N - problem.XN) + e : (problem.N - problem.XN) + e + 1); ...
        solutionCoefficients(problem.N  + (e - 2) * modes + 2 : problem.N...
        + (e - 2) * modes + length(indexLocalEnrichedNodes)  * modes + 1 )];
end

projectedCoefficients = projectionOperator * coefficients ;

if derivative == 1
    for i=1:length(localCoords)
        projectionOperator(i,1:size(N,2)+size(F,2)) = [N(i,:), F(i,:)];
    end
    
    temperature = projectionOperator * coefficients ;
    for i=1:numel(projectedCoefficients)
        projectedCoefficients(i) = projectedCoefficients(i) * ...
            problem.k(globalCoords, t, temperature(i));
    end
end

end