function [ M, f ] = assemblyL2ProjectionMultiscalehIGA( solutionBaseCoefficients,...
    solutionOverlayCoefficients, baseProblem, overlayProblem,...
    previousBaseProblem, previousOverlayProblem, integrationOrder,...
    initialTemperature, layer )
%ASSEMBLYL2PROJECTIONMULTISCALEPODXIGA assemble the matrix M and the vector f for
%the L2 projection system
%Input:
%solutionCoefficients = coefficients of the solutions
%baseProblem = poisson problem struct of the base mesh
%overlayProblem = poisson problem struct of the overlay mesh
%previousBaseProblem = problem of the previous layer base mesh
%layerLength = length of the layer
%initialtemperature = initial temperature of the powder
%layer = index of the actual layer
%Output:
%M = L2 projection matrix
%f = L2 projection rhs vector

%% Allocate matrices
%base  matrix
M_bb = zeros(baseProblem.gdof, baseProblem.gdof);
%overlay matrix
M_oo = zeros(overlayProblem.gdof, overlayProblem.gdof);
%coupling matrix
M_ob = zeros(overlayProblem.gdof, baseProblem.gdof);
%base  vector
f_b = zeros(baseProblem.gdof, 1);
%overlay vector
f_o = zeros(overlayProblem.gdof, 1);

%gauss points base elements
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

%loop over non-overlapped base elements of previous mesh
for e=1:previousBaseProblem.N-1
    
    ldof_b = previousBaseProblem.p+1;
    
    Xi1 = baseProblem.knotVector(e + baseProblem.p);
    Xi2 = baseProblem.knotVector(e + baseProblem.p + 1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        parametricGPCoord = mapParentToLocal(rGP(iGP),Xi1, Xi2);
        
        [Nsplines, ~] = BsplinesShapeFunctionsAndDerivatives(parametricGPCoord, baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        
        %% Integrate IGA block
        %Base vector
        f_b(baseProblem.LM(e,1:ldof_b)) = f_b(baseProblem.LM(e,1:ldof_b)) + Nsplines(baseProblem.LM(e, 1:ldof_b))' * evaluateTemperature(e, rGP(iGP), previousOverlayProblem,...
            previousBaseProblem, solutionBaseCoefficients) * wGP(iGP) * baseProblem.F_map(Xi1,Xi2) * detJacobianX_Xi;
        
        %Base matrix
        M_bb(previousBaseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = M_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            (Nsplines(baseProblem.LM(e, 1:ldof_b))' * Nsplines(baseProblem.LM(e, 1:ldof_b))) * wGP(iGP) * baseProblem.F_map(Xi1,Xi2) * detJacobianX_Xi;
    end
    
end

%Last layer of the previous mesh using composed integration
[rGP, wGP] = gaussPoints( 2 );
numberOfIntegrationPoints = length(rGP);

%loop over overlapping elements of the previous mesh
for e=1:previousOverlayProblem.N
    
    ldof_b = previousBaseProblem.p+1;
    
    %overlay mesh element boundaries
    X1 = previousOverlayProblem.coords(e);
    X2 = previousOverlayProblem.coords(e+1);
    
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);
        
        [Nsplines, ~] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        
        %% Integrate IGA block
        %Base vector
        f_b(previousBaseProblem.LM(end,1:ldof_b)) = f_b(previousBaseProblem.LM(end,1:ldof_b)) + Nsplines(layer-1:end-1)' * evaluateTemperatureComposed(e, rGP(iGP), previousOverlayProblem, previousBaseProblem,...
            solutionOverlayCoefficients, solutionBaseCoefficients) * wGP(iGP) * previousBaseProblem.F_map(X1,X2);
        
        %Base matrix
        M_bb(previousBaseProblem.LM(end, 1:ldof_b), previousBaseProblem.LM(end, 1:ldof_b)) = M_bb(previousBaseProblem.LM(end, 1:ldof_b), previousBaseProblem.LM(end, 1:ldof_b)) +...
            (Nsplines(layer-1:end-1)' * Nsplines(layer-1:end-1)) * wGP(iGP) * previousBaseProblem.F_map(X1,X2);
        
        
    end
end

%New layer at initial temperature
%Composed Integration------------------------------------------------------

%gauss points composed integration
%Composed integration allows to integrate discontinuous functions over a
%domain. The quadrature points are distributed over the smaller sub-domain
%where the function is still continuous. This technique is essential to
%correctly integrate the basis function of the overlapping mesh.

[rGPComposed, wGPComposed] = gaussPoints( 2 ); %linear elements in the overlay mesh
numberOfComposedIntegrationPoints = length(rGPComposed);

%loop over non-overlapped base elements
for e=1:overlayProblem.N
    
    ldof_b = baseProblem.p+1;
    
    %overlay mesh element boundaries
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    
    % Gauss integration
    for iGP = 1:numberOfComposedIntegrationPoints
        
        [N, ~] = shapeFunctionsAndDerivatives(rGPComposed(iGP));
        globalGPCoord = mapLocalToGlobal(rGPComposed(iGP),X1, X2);
        
        [Nsplines, ~] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        
        %% Integrate IGA block
        %Base vector
        f_b(baseProblem.LM(end,1:ldof_b)) = f_b(baseProblem.LM(end,1:ldof_b)) + Nsplines(layer:end)' * initialTemperature...
            * wGPComposed(iGP) * detJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %Base matrix
        M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) = M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
             (Nsplines(layer:end)' * Nsplines(layer:end)) * wGPComposed(iGP) * detJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %% Integrate Coupling block
        %Coupling matrix
        M_ob(overlayProblem.LM(e,:), baseProblem.LM(end,:)) = M_ob(overlayProblem.LM(e,:), baseProblem.LM(end,:)) +...
             (N' * Nsplines(layer:end)) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %% Integrate FEM block
        %Overlay vector
        f_o(overlayProblem.LM(e,:)) = f_o(overlayProblem.LM(e,:)) + N' * initialTemperature * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Overlay matrix
        M_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) = M_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
             (N' * N) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
    end
end

[M_oo, f_o] = constrainOverlayMeshL2Projection( M_oo, f_o, overlayProblem );

M = [M_bb, M_ob'; M_ob, M_oo];
f = [f_b; f_o];

end

%% Evaluate temperature in the non-overlapped base mesh
function [ projectedCoefficients ] = evaluateTemperature(e, x, overlayProblem,...
    baseProblem, solutionBaseCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   overlayProblem = problem struct of the overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh


numberOfProjectionPoints = length(x);
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
Nspline = zeros(length(x), size(baseProblem.LM, 2));

Xi1 = baseProblem.knotVector(e + baseProblem.p);
Xi2 = baseProblem.knotVector(e + baseProblem.p + 1);

localCoord = mapParentToLocal(x,Xi1, Xi2);

for k=1:length(x)
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(localCoord, baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(Nspline,2)) = NIga(baseProblem.LM(e,:));
end

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) *...
    solutionBaseCoefficients(baseProblem.LM(e,:));
projectedCoefficients = projectedBaseCoefficients;

end

%% Evaluate temperature for composed inetgartion
function [ projectedCoefficients ] = evaluateTemperatureComposed(e, x, overlayProblem,...
    baseProblem, solutionOverlayCoefficients, solutionBaseCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   overlayProblem = problem struct of the overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionOverlayCoefficients = temeprature distribution of the previous overlay mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh


numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2));
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));

N = zeros(length(x), 2);
Nspline = zeros(length(x), size(baseProblem.LM, 2));

localCoords = x;
globalCoord = mapLocalToGlobal(x,overlayProblem.coords(e), overlayProblem.coords(e+1));

for k=1:length(x)
    [projectionOperator(k,1:size(N,2)), ~] = shapeFunctionsAndDerivatives(localCoords(k));    
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric...
        (globalCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(Nspline,2)) = NIga(baseProblem.LM(end,:));
end

%evaluate the overlay solution
projectedOverlayCoefficients = projectionOperator * solutionOverlayCoefficients(overlayProblem.LM(e,:));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(baseProblem.LM(end,:));

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end