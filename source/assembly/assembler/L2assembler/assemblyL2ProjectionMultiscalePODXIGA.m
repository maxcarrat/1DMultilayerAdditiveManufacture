function [ M, f ] = assemblyL2ProjectionMultiscalePODXIGA( previousBaseTemperature,...
    previousOverlayTemperature, baseProblem, overlayProblem, integrationOrder,...
    integrationModalOrder, layer, layerLength,...
    initialTemperature, previousBaseProblem, previousOverlayProblem )
%ASSEMBLYL2PROJECTIONMULTISCALEPODXIGA assemble the matrix M and the vector f for
%the L2 projection system
%Input:
%previousBaseTemperature = temperature coefficients of the base mesh
%solution in the previous time step
%previousOverlayTemperature = temperature coefficients of the overlay mesh
%solution in the previous time step
%baseProblem = struct which defines the initial boundary value problem on
%the base mesh
%overlayProblem = problem struct which defines the IBVP on the overlay mesh
%integrationOrder = number of quadrature points
%integrationModalOrder = numeber of quadratur points in enriched elements
%layerLength = length of bar layers
%initialTemperature = initial temperature of the system
%previousBaseProblem = struct of the poisson problem in the previous layer
%Output:
%M = L2 projection matrix integral(N'N detJx detJxi)
%f = L2 rhs integral(N' functionToProject detJx detJxi)

%% Allocate matrices

%Overlay FEM  matrix
M_FEM = zeros(overlayProblem.overDofs,overlayProblem.overDofs);
%Overlay FEM  vector
f_FEM = zeros(overlayProblem.overDofs, 1);
%Overlay XFEM  matrix
M_XFEM = zeros(overlayProblem.Xdof,overlayProblem.Xdof);
%Overlay XFEM vector
f_XFEM = zeros(overlayProblem.Xdof, 1);
%Coupling overlay matrix
M_Coupling = zeros(overlayProblem.Xdof,overlayProblem.overDofs);
%base  matrix
M_bb = zeros(baseProblem.gdof, baseProblem.gdof);
%coupling matrix
M_ob = zeros(overlayProblem.gdof, baseProblem.gdof);
%base  vector
f_b = zeros(baseProblem.gdof, 1);

%Gauss integration
[rGP, wGP] = gaussPoints( integrationOrder );
[rGPXFEM, wGPXFEM] = gaussPoints( integrationModalOrder );

numberOfIntegrationPoints = length(rGP);
numberOfModalIntegrationPoints = length(rGPXFEM);

modes = overlayProblem.modes;
incrementModes = 0;

%On active elements use the refined domain as integration domain
refinedDofs = previousOverlayProblem.N + 1;
integrationDomain = linspace(-1, 1, ceil(refinedDofs/overlayProblem.XN));

%loop over non-overlapped base elements of previous mesh
for e=1:previousBaseProblem.N-1
    
    ldof_b = previousBaseProblem.p+1;
    
    Xi1 = baseProblem.knotVector(e + baseProblem.p);
    Xi2 = baseProblem.knotVector(e + baseProblem.p + 1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        parametricGPCoord = mapParentToLocal(rGP(iGP),Xi1, Xi2);
        
        [Nsplines, ~] = BsplinesShapeFunctionsAndDerivatives(parametricGPCoord, baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = previousBaseProblem.coords(end) - previousBaseProblem.coords(1);
        
        %% Integrate IGA block
        %Base vector
        f_b(baseProblem.LM(e,1:ldof_b)) = f_b(baseProblem.LM(e,1:ldof_b)) + Nsplines(baseProblem.LM(e, 1:ldof_b))' * evaluateTemperature(e, rGP(iGP),...
            previousBaseProblem, previousBaseTemperature) * wGP(iGP) * baseProblem.F_map(Xi1,Xi2) * detJacobianX_Xi;
        
        %Base matrix
        M_bb(previousBaseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = M_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            (Nsplines(baseProblem.LM(e, 1:ldof_b))' * Nsplines(baseProblem.LM(e, 1:ldof_b))) * wGP(iGP) * baseProblem.F_map(Xi1,Xi2) * detJacobianX_Xi;
    end
    
end

%Last layer of the previous mesh using composed integration
[rGP, wGP] = gaussPoints( 2 );
numberOfIntegrationPoints = length(rGP);

%% ------------------------------------------------------------------------
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
        f_b(previousBaseProblem.LM(end,1:ldof_b)) = f_b(previousBaseProblem.LM(end,1:ldof_b)) + Nsplines(layer-1:end-1)'...
            * evaluateTemperatureComposed(e, rGP(iGP), previousOverlayProblem, overlayProblem, previousBaseProblem,...
            previousOverlayTemperature, previousBaseTemperature) * wGP(iGP) * previousBaseProblem.F_map(X1,X2);
        
        %Base matrix
        M_bb(previousBaseProblem.LM(end, 1:ldof_b), previousBaseProblem.LM(end, 1:ldof_b)) = ...
            M_bb(previousBaseProblem.LM(end, 1:ldof_b), previousBaseProblem.LM(end, 1:ldof_b)) +...
            (Nsplines(layer-1:end-1)' * Nsplines(layer-1:end-1)) * wGP(iGP) * previousBaseProblem.F_map(X1,X2);
        
        
    end
end

%% ------------------------------------------------------------------------
%New layer at initial temperature
%Composed Integration------------------------------------------------------

%gauss points composed integration
%Composed integration allows to integrate discontinuous functions over a
%domain. The quadrature points are distributed over the smaller sub-domain
%where the function is still continuous. This technique is essential to
%correctly integrate the basis function of the overlapping mesh.

[rGPComposed, wGPComposed] = gaussPoints( integrationModalOrder );

for e=1:overlayProblem.N
    
    ldof = 2;
    ldof_b = baseProblem.p+1;

    %overlay mesh element boundaries
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    elementGloabalCoords = [X1, X2];
    
    %% overlay element mesh
    
    %In this element the temperature values are constantly equal to the
    %initial temperature, but they are projected onto an enriched
    %element.
    
    % Gauss integration
    for iGP = 1:numberOfModalIntegrationPoints
        
        %evaluate linear shape functions at integration points
        [N, ~] = shapeFunctionsAndDerivatives(rGPXFEM(iGP));
        
        %Jacobian from local to global
        detJacobianLocalToGlobal = overlayProblem.F_map(X1, X2);
        
        %global coords
        globalGPCoord = mapLocalToGlobal(rGPComposed(iGP),X1, X2);
        
        %BSplines
        [Nspline, ~] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        
        %% Integrate IGA block
        %Base vector
        f_b(baseProblem.LM(end,1:ldof_b)) = f_b(baseProblem.LM(end,1:ldof_b)) + Nspline(layer:end)' * initialTemperature...
            * wGPComposed(iGP) *  detJacobianLocalToGlobal;
        
        %Base matrix
        M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) = M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
            (Nspline(layer:end)' * Nspline(layer:end)) * wGPComposed(iGP) * detJacobianLocalToGlobal;
        
        %% Integrate FEM block
        
        % vector
        f_FEM(overlayProblem.LM(e,:)) = f_FEM(overlayProblem.LM(e,:)) + N' * initialTemperature...
            * wGPXFEM(iGP) * detJacobianLocalToGlobal;
        
        % matrix
        M_FEM(overlayProblem.LM(e, :), overlayProblem.LM(e, :)) = M_FEM(overlayProblem.LM(e, :), overlayProblem.LM(e, :)) +...
            (N' * N) * wGPXFEM(iGP) * detJacobianLocalToGlobal;
        
        elementEnrichedIndex = e - (overlayProblem.N - overlayProblem.XN);
        
        if elementEnrichedIndex == 1
            indexLocalEnrichedNodes = 2;
        else
            indexLocalEnrichedNodes = [1, 2];
        end
        
        modalDofs = length(indexLocalEnrichedNodes)*modes;
        incrementModes = (length(indexLocalEnrichedNodes) - 1) * modes;
        
        %loop over sub-domains
        for integrationSubDomainIndex =...
                1 : ceil(refinedDofs/overlayProblem.XN)-1
            
            %if the GP is inside the sub-domain
            if rGPXFEM(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                    rGPXFEM(iGP) > integrationDomain(integrationSubDomainIndex)
                
                PODCoefficients = overlayProblem.reductionOperator(...
                    (elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)+1:(elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)...
                    + ceil(refinedDofs/overlayProblem.XN),:);
                
                %evaluate BSplines and POD-modal enrichment functions
                %@GP
%                 [ N, ~] = shapeFunctionsAndDerivativesSubElements( rGPXFEM(iGP), ...
%                     integrationSubDomainIndex, subDomainShapeFunctionCoefficients, integrationDomain );
                [ N, ~] = shapeFunctionsAndDerivatives( rGPXFEM(iGP));                
                [ F, ~ ] = PODModesAndDerivativesMultiscale( rGPXFEM(iGP), elementGloabalCoords, modes,...
                    PODCoefficients, integrationDomain,...
                    integrationSubDomainIndex, indexLocalEnrichedNodes );
                                
                %% Integrate Coupling block
                % matrix
                M_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) =...
                    M_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) +...
                    F' * N * wGPXFEM(iGP) *  detJacobianLocalToGlobal;
                
                %% Integrate XFEM block
                % vector
                f_XFEM(overlayProblem.LME(elementEnrichedIndex,1:modalDofs)) = f_XFEM(overlayProblem.LME(elementEnrichedIndex,1:modalDofs)) +  F' *...
                    initialTemperature * wGPXFEM(iGP) * detJacobianLocalToGlobal;
                
                % matrix
                M_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                    M_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                    ( F' * F ) * wGPXFEM(iGP) * detJacobianLocalToGlobal;
                
             %% Base/Overlay Coupling block
                N = [N, F];

                %Coupling matrix
                M_ob(overlayProblem.LMBC(e,1:modalDofs+2), baseProblem.LM(end,:)) = M_ob(overlayProblem.LMBC(e,1:modalDofs+2), baseProblem.LM(end,:)) +...
                    (N' * Nspline(baseProblem.LM(end,1:ldof_b))) * wGPComposed(iGP) * detJacobianLocalToGlobal;    
                
                break;
            end
        end
        
    end
end

%% Assembly
M_oo = [M_FEM, M_Coupling'; M_Coupling, M_XFEM];
f_o = [f_FEM; f_XFEM];

%Constrain overlay mesh
[M_oo, f_o] = constrainOverlayMeshL2Projection( M_oo, f_o, overlayProblem );

M = [M_bb, M_ob'; M_ob, M_oo];
f = [f_b; f_o];

end


%% Evaluate temperature in the non-overlapped base mesh
function [ projectedCoefficients ] = evaluateTemperature(e, x,...
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
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(baseProblem.LM(e,:));

projectedCoefficients = projectedBaseCoefficients;

end

%% Evaluate temperature for composed integration
function [ projectedCoefficients ] = evaluateTemperatureComposed(e, x, overlayProblem,...
   overlayXProblem, baseProblem, solutionOverlayCoefficients, solutionBaseCoefficients)
% EVALUATETEMPERATURECOMPOSED project the previous solution onto the element
%using composed integration.
%   e = element index
%   x = post-processing mesh
%   overlayProblem = problem struct of the initial overlay mesh
%   overlayXProblem = problem struct of the POD-enriched overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionOverlayCoefficients = temeprature distribution of the previous overlay mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh


numberOfProjectionPoints = length(x);
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
% integrationDomain = linspace(-1, +1, numel(overlayXProblem.reductionOperator));
% 
% localCoords = x;
globalCoord = mapLocalToGlobal(x,overlayXProblem.coords(1), overlayXProblem.coords(end));

% elementEnrichedIndex = e - (overlayProblem.N - overlayXProblem.XN);
% 
% if elementEnrichedIndex == 1
%     indexLocalEnrichedNodes = 2;
% else
%     indexLocalEnrichedNodes = [1, 2];
% end
% 
% modalDofs = length(indexLocalEnrichedNodes)*overlayXProblem.modes;
% 
% projectionOperator = zeros(numberOfProjectionPoints, size(overlayXProblem.LM, 2) + modalDofs);

for k=1:length(x)
%     [ N, ~] = shapeFunctionsAndDerivatives( localCoords(k) );
%     
%     [F, ~] = PODModesAndDerivativesMultiscale(localCoords(k), overlayXProblem.modes,...
%         overlayXProblem.reductionOperator, integrationDomain, e, indexLocalEnrichedNodes);
%     
%     projectionOperator(k,1:size(N, 2) + size(F, 2)) =...
%         [N, F];    
    
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric...
        (globalCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
    
    projectionOperatorSpline(k,1:size(NIga,2)) = NIga(:);
end

% %evaluate the overlay solution
% projectedOverlayCoefficients = projectionOperator * solutionOverlayCoefficients(overlayXProblem.LMBC(e,1:modalDofs+2));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline * solutionBaseCoefficients;

% projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;
projectedCoefficients = projectedBaseCoefficients;

end




