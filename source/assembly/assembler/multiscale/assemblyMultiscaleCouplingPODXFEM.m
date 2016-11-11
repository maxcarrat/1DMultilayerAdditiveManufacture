function [ M_ob, K_ob ] = assemblyMultiscaleCouplingPODXFEM( initialOverlayProblem, overlayProblem,...
    baseProblem, time, integrationOrder, overlaySolutionCoefficients,...
    baseSolutionCoefficients, lastConvergentBaseSolution, lastConvergentOverlaySolution )
%ASSEMBLYMULTISCALECOUPLINGPODXFEM  assembles the mass, the conductivity matrix
%and load vector
%Input:
%overlayProblem = problem struct of the overlay mesh
%baseProblem = problem struct of the base mesh
%time = current time
%integrationOrder = number of quadrature points for the integration
%overlaySolutionCoefficients = coefficients of the overlay solution
%baseSolutionCoefficients = coefficients of the base solution
%lastConvergentBaseSolution = last convergent base solution
%lastConvergentOverlaySolution = last convergent overlay solution
%Output:
%K_ob = coupling conductivity matrix
%M_ob = overlay capacity matrix

%% Allocate matrices
%coupling K matrix and Kprime matrix
K_ob = zeros(overlayProblem.gdof,baseProblem.gdof);
%coupling M matrix
M_ob = zeros(overlayProblem.gdof,baseProblem.gdof);

%On active elements use the refined domain as integration domain
refinedDofs = initialOverlayProblem.N + 1;
integrationDomain = linspace(-1, 1, ceil(refinedDofs/overlayProblem.XN));

%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
invDetJacobianX_Xi = 1 / detJacobianX_Xi;

for e=1:overlayProblem.N
    
    ldof_b = baseProblem.p + 1;
    
    %overlay mesh element boundaries
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    elementGloabalCoords = [X1, X2];

    detJacobianLocalToGlobal = overlayProblem.F_map(X1, X2);
    inverseJacobianLocalToGlobal = overlayProblem.B_map(X1,X2);
    
    elementEnrichedIndex = e - (overlayProblem.N - overlayProblem.XN);
    
    if elementEnrichedIndex == 1
        indexLocalEnrichedNodes = 2;
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    modalDofs = length(indexLocalEnrichedNodes) * overlayProblem.modes;
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);
        
        [Nspline, Bspline] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        Bspline = Bspline * invDetJacobianX_Xi;
        
        %% Integrate Coupling block using composed integration
        %Composed integration allows to integrate discontinuous functions over a
        %domain. The quadrature points are distributed over the smaller sub-domain
        %where the function is still continuous. This technique is essential to
        %correctly integrate the basis function of the overlapping mesh.
        
        %loop over sub-domains
        for integrationSubDomainIndex =...
                1 : ceil(refinedDofs/overlayProblem.XN) - 1
            
            %if the GP is inside the sub-domain
            if rGP(iGP) <= integrationDomain(integrationSubDomainIndex + 1) &&...
                    rGP(iGP) > integrationDomain(integrationSubDomainIndex)
                
                PODCoefficients = overlayProblem.reductionOperator(...
                    (elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)+1:(elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)...
                    + ceil(refinedDofs/overlayProblem.XN),:);
                
                %evaluate BSplines and POD-modal enrichment functions
                %@GP
                [ N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
                B = B * inverseJacobianLocalToGlobal;
                
                [ F, G ] = PODModesAndDerivativesMultiscale( rGP(iGP), elementGloabalCoords, overlayProblem.modes,...
                    PODCoefficients, integrationDomain,...
                    integrationSubDomainIndex, indexLocalEnrichedNodes );
                
                N = [N, F];
                B = [B, G];
                
                %Capacity matrix
                M_ob(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2), baseProblem.LM(end,:)) = M_ob(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2), baseProblem.LM(end,:)) +...
                    overlayProblem.heatCapacity( mapLocalToGlobal(rGP(iGP), X1, X2),...
                    evaluateTemperature(integrationSubDomainIndex, e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
                    evaluateTemperature(integrationSubDomainIndex, e, rGP(iGP), overlayProblem, baseProblem, lastConvergentOverlaySolution, lastConvergentBaseSolution))...
                    * (N' * Nspline(baseProblem.LM(end,1:ldof_b))) * wGP(iGP) * detJacobianLocalToGlobal;

                %Conductivity matrix
                K_ob(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2), baseProblem.LM(end,:)) = K_ob(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2), baseProblem.LM(end,:)) +...
                    overlayProblem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
                    time, evaluateTemperature(integrationSubDomainIndex, e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
                    * (B' * Bspline(baseProblem.LM(end,1:ldof_b))) * wGP(iGP) * detJacobianLocalToGlobal;
                break;
            end
        end
    end
    
end

end

function [ projectedCoefficients ] = evaluateTemperature(subDomainIndex, element, x, overlayProblem,...
    baseProblem, solutionOverlayCoefficients, solutionBaseCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = subDomain index
%   x = post-processing mesh
%   overlayProblem = problem struct of the overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionOverlayCoefficients = temeprature distribution of the previous overlay mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh

numberOfProjectionPoints = length(x);
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
refinedDofs = size(overlayProblem.reductionOperator,1);
integrationDomain = linspace(-1, +1, ceil(refinedDofs/overlayProblem.XN));

X1 = overlayProblem.coords(element);
X2 = overlayProblem.coords(element +1);
elementGloabalCoords = [X1, X2];

localCoords = x;
globalCoord = mapLocalToGlobal(x, X1, X2);

parametricGPCoord = mapGlobalToParametric(...
    globalCoord, baseProblem.coords(1), baseProblem.coords(end));

elementEnrichedIndex = element - (overlayProblem.N - overlayProblem.XN);

if elementEnrichedIndex == 1
    indexLocalEnrichedNodes = 2;
else
    indexLocalEnrichedNodes = [1, 2];
end

PODCoefficients = overlayProblem.reductionOperator(...
    (elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)+1:(elementEnrichedIndex-1)*floor(refinedDofs/overlayProblem.XN)...
    + ceil(refinedDofs/overlayProblem.XN),:);

modalDofs = length(indexLocalEnrichedNodes)*overlayProblem.modes;

projectionOperator = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2) + modalDofs);

for k=1:length(x)
    [ N, ~] = shapeFunctionsAndDerivatives(localCoords(k));
    
    [F, ~] = PODModesAndDerivativesMultiscale(localCoords(k), elementGloabalCoords, overlayProblem.modes,...
        PODCoefficients, integrationDomain, subDomainIndex, indexLocalEnrichedNodes);
    
    projectionOperator(k,1:size(N, 2) + size(F, 2)) =...
        [N, F];
    
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(parametricGPCoord(k),...
        baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(NIga,2)) = NIga(:);
end

%evaluate the overlay solution
projectedOverlayCoefficients = projectionOperator(:,:) * solutionOverlayCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(:);

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end
