function [ M_oo, Mprime_oo, K_oo, Kprime_oo, f_o ] = assemblyMultiscaleOverlayPODXFEM( initialOverlayProblem, overlayProblem, baseProblem, time, integrationOrder,...
    overlayIntegrationOrder, overlaySolutionCoefficients, baseSolutionCoefficients, lastConvergedBaseSolution, lastConvergentOverlaySolution )
%ASSEMBLYMULTISCALEOVERLAYPODXFEM assembles the mass, the conductivity matrix
%and load vector
%Input:
%overlayProblem = problem struct of the overlay mesh
%baseProblem = problem struct of the base mesh
%time = current time
%integrationOrder = number of quadrature points for the integration
%overlayIntegrationOrder = number of quadrature points for the integration
%of POD modes
%overlaySolutionCoefficients = coefficients of the overlay solution
%baseSolutionCoefficients = coefficients of the base solution
%lastConvergedBaseSolution = last convergent base solution
%lastConvergentOverlaySolution = last convergent overlay solution
%Output:
%K_oo = overlay conductivity matrix
%Kprime_oo = derivative of the diffussion matrix
%M_oo = overlay capacity matrix
%Mprime_oo = derivative of the capacity matrix
%f_oo = rhs overlay vector

%% Allocate matrices
%overlay K matrix and Kprime matrix
K_oo = zeros(overlayProblem.gdof,overlayProblem.gdof);
Kprime_oo = zeros(overlayProblem.gdof,overlayProblem.gdof);
%overlay M matrix Mprime matrix
M_oo = zeros(overlayProblem.gdof,overlayProblem.gdof);
Mprime_oo = zeros(overlayProblem.gdof,overlayProblem.gdof);
%overlay rhs vector
f_o = zeros(overlayProblem.gdof, 1);
%Overlay FEM  matrix
M_FEM = zeros(overlayProblem.overDofs,overlayProblem.overDofs);
K_FEM = zeros(overlayProblem.overDofs,overlayProblem.overDofs);
MPrime_FEM = zeros(overlayProblem.overDofs,overlayProblem.overDofs);
KPrime_FEM = zeros(overlayProblem.overDofs,overlayProblem.overDofs);
%Overlay FEM  vector
f_FEM = zeros(overlayProblem.overDofs, 1);
%Overlay XFEM  matrix
M_XFEM = zeros(overlayProblem.Xdof,overlayProblem.Xdof);
K_XFEM = zeros(overlayProblem.Xdof,overlayProblem.Xdof);
MPrime_XFEM = zeros(overlayProblem.Xdof,overlayProblem.Xdof);
KPrime_XFEM = zeros(overlayProblem.Xdof,overlayProblem.Xdof);
%Overlay XFEM vector
f_XFEM = zeros(overlayProblem.Xdof, 1);
%Coupling overlay matrix
M_Coupling = zeros(overlayProblem.Xdof,overlayProblem.overDofs);
K_Coupling = zeros(overlayProblem.Xdof,overlayProblem.overDofs);
MPrime_Coupling = zeros(overlayProblem.Xdof,overlayProblem.overDofs);
KPrime_Coupling = zeros(overlayProblem.Xdof,overlayProblem.overDofs);


%On active elements use the refined domain as integration domain
refinedDofs = initialOverlayProblem.N + 1;
integrationDomain = linspace(-1, 1, ceil(refinedDofs/overlayProblem.XN));

%gauss points
[rGP, wGP] = gaussPoints( overlayIntegrationOrder );
numberOfIntegrationPoints = length(rGP);

modes = overlayProblem.modes;

%loop over the overlay elements
for e=1:overlayProblem.N
    
    ldof_b = baseProblem.p + 1;

    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    elementGloabalCoords = [X1, X2];
    
    elementEnrichedIndex = e - (overlayProblem.N - overlayProblem.XN);
    
    if elementEnrichedIndex == 1
        indexLocalEnrichedNodes = 2;
    else
        indexLocalEnrichedNodes = [1, 2];
    end
    
    modalDofs = length(indexLocalEnrichedNodes)*modes;
    
    %Jacobian from local to global
    detJacobianLocalToGlobal = overlayProblem.F_map(X1, X2);
    inverseJacobianLocalToGlobal = overlayProblem.B_map(X1,X2);
    
    %Jacobian and Inverse from parametric to global
    detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
    inverseJacobianX_Xi = 1 / detJacobianX_Xi;
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);
        
        parametricGPCoord = mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end));
        
        [Nsplines, Bspline] = BsplinesShapeFunctionsAndDerivatives(parametricGPCoord,...
            baseProblem.p, baseProblem.knotVector);
        Bspline = Bspline * inverseJacobianX_Xi;
        
        deltaBaseSolution = -(lastConvergedBaseSolution - baseSolutionCoefficients);
        deltaOverlaySolution = (lastConvergentOverlaySolution - overlaySolutionCoefficients);
        
        
        
        %% Composed Integration -------------------------------------------
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
                
                %Evaluate POD modal basis
                [ N, B ] = shapeFunctionsAndDerivatives(rGP(iGP));
                B = B * inverseJacobianLocalToGlobal;
                
                [ F, G ] = PODModesAndDerivativesMultiscale( rGP(iGP), elementGloabalCoords, modes,...
                    PODCoefficients, integrationDomain,...
                    integrationSubDomainIndex, indexLocalEnrichedNodes );
                
                NTotal = [N, F];
                BTotal = [B, G];
                
                %% Integrate FEM block --------------------------------------------
                
                temperature = evaluateTemperature(integrationSubDomainIndex, e, rGP(iGP), overlayProblem, baseProblem,...
                    overlaySolutionCoefficients, baseSolutionCoefficients);
                
                lastTemperature = evaluateTemperature(integrationSubDomainIndex, e, rGP(iGP), overlayProblem, baseProblem,...
                    lastConvergentOverlaySolution, lastConvergedBaseSolution);

                
                % vector
                f_FEM(overlayProblem.LM(e,:)) = f_FEM(overlayProblem.LM(e,:))...
                    + N' * overlayProblem.rhs(globalGPCoord, time)...
                    * wGP(iGP) * detJacobianLocalToGlobal;
      
                % matrix
                M_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) =...
                    M_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
                    overlayProblem.heatCapacity(globalGPCoord, temperature, lastTemperature)* ...
                    (N' * N) * wGP(iGP) *  detJacobianLocalToGlobal;
             
                K_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) =...
                    K_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
                    overlayProblem.k(globalGPCoord, time, temperature) *...
                    (B' * B) * wGP(iGP) *  detJacobianLocalToGlobal;
            
                MPrime_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) =...
                    MPrime_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
                    overlayProblem.heatCapacityDerivative(globalGPCoord, temperature, lastTemperature)*...
                    (N' * N) * (Nsplines(baseProblem.LM(end, 1:ldof_b)) * deltaBaseSolution(baseProblem.LM(end, 1:ldof_b))...
                    + NTotal * deltaOverlaySolution(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) * detJacobianLocalToGlobal;
            
                KPrime_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) =...
                    KPrime_FEM(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
                    overlayProblem.kDerivative(globalGPCoord, time, temperature) *...
                    (B' * B) * (Bspline(baseProblem.LM(end, 1:ldof_b)) * baseSolutionCoefficients(baseProblem.LM(end, 1:ldof_b))...
                    + BTotal * overlaySolutionCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
                


                
                %% Integrate Coupling block -------------------------------
                % matrix
                M_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) =...
                    M_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) +...
                    overlayProblem.heatCapacity( globalGPCoord, temperature, lastTemperature) *...
                    F' * N * wGP(iGP) *  detJacobianLocalToGlobal;
            
                K_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) =...
                    K_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) +...
                    overlayProblem.k(globalGPCoord, time, temperature) *...
                    G' * B * wGP(iGP) *  detJacobianLocalToGlobal;
               
                MPrime_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) =...
                    MPrime_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) +...
                    overlayProblem.heatCapacityDerivative( globalGPCoord, temperature, lastTemperature) *...
                    F' * N * (Nsplines(baseProblem.LM(end, 1:ldof_b)) * deltaBaseSolution(baseProblem.LM(end, 1:ldof_b))...
                    + NTotal * deltaOverlaySolution(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
             
                KPrime_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) =...
                    KPrime_Coupling(overlayProblem.LME(elementEnrichedIndex,1:modalDofs), overlayProblem.LMC(elementEnrichedIndex, :)) +...
                    overlayProblem.kDerivative(globalGPCoord, time, temperature) *...
                    G' * B * (Bspline(baseProblem.LM(end, 1:ldof_b)) * baseSolutionCoefficients(baseProblem.LM(end, 1:ldof_b))...
                    + BTotal * overlaySolutionCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
                
                %% Integrate XFEM block -----------------------------------
                % vector
                f_XFEM(overlayProblem.LME(elementEnrichedIndex,1:modalDofs)) = f_XFEM(overlayProblem.LME(elementEnrichedIndex,1:modalDofs))...
                    + F' * overlayProblem.rhs(globalGPCoord, time)...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
                % matrix
                M_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                    M_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                    overlayProblem.heatCapacity( globalGPCoord, temperature, lastTemperature) *...
                    (F' * F) * wGP(iGP) *  detJacobianLocalToGlobal;
               
                K_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                    K_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                    overlayProblem.k(globalGPCoord, time, temperature) *...
                    (G' * G) * wGP(iGP) *  detJacobianLocalToGlobal;
              
                MPrime_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                    MPrime_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                    overlayProblem.heatCapacityDerivative( globalGPCoord, temperature, lastTemperature) *...
                    (F' * F) * (Nsplines(baseProblem.LM(end, 1:ldof_b)) * deltaBaseSolution(baseProblem.LM(end, 1:ldof_b))...
                    + NTotal * deltaOverlaySolution(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
         
                KPrime_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) =...
                    KPrime_XFEM(overlayProblem.LME(elementEnrichedIndex, 1:modalDofs), overlayProblem.LME(elementEnrichedIndex, 1:modalDofs)) +...
                    overlayProblem.kDerivative(globalGPCoord,time, temperature) *...
                    (G' * G) * (Bspline(baseProblem.LM(end, 1:ldof_b)) * baseSolutionCoefficients(baseProblem.LM(end, 1:ldof_b))...
                    + BTotal * overlaySolutionCoefficients(overlayProblem.LMBC(elementEnrichedIndex,1:modalDofs+2)))...
                    * wGP(iGP) *  detJacobianLocalToGlobal;
                
                break;
            end
        end %end loop over subDomains
        
        %Constrain Xtended Dofs
%         [K_XFEM,KPrime_XFEM, f_XFEM] = constrainXtendedMesh(K_XFEM, KPrime_XFEM, f_XFEM, overlayProblem);
        
        %external heat source
        f_o = [f_FEM; f_XFEM];
        
        %Capacity matrix
        M_oo = [M_FEM, M_Coupling'; M_Coupling, M_XFEM];
        Mprime_oo = [MPrime_FEM, MPrime_Coupling'; MPrime_Coupling, MPrime_XFEM];
        
        %Conductivity matrix
        K_oo = [K_FEM, K_Coupling'; K_Coupling, K_XFEM];
        Kprime_oo = [KPrime_FEM, KPrime_Coupling'; KPrime_Coupling, KPrime_XFEM];
        
    end %end loop over integration points
end %end loop over elements
end

function [ projectedCoefficients ] = evaluateTemperature(subDomainIndex, element, x, overlayXProblem,...
    baseProblem, solutionOverlayCoefficients, solutionBaseCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   subDomainIndex = subDomain index
%   element = element index
%   x = evaluation local coordinates
%   overlayXProblem = problem struct of the eXtended overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionOverlayCoefficients = temeprature distribution of the previous overlay mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh

numberOfProjectionPoints = length(x);
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));
refinedDofs = size(overlayXProblem.reductionOperator,1);
integrationDomain = linspace(-1, +1, ceil(refinedDofs/overlayXProblem.XN));

X1 = overlayXProblem.coords(element);
X2 = overlayXProblem.coords(element + 1);
elementGlobalCoords = [X1, X2];

localCoords = x;
globalCoord = mapLocalToGlobal(localCoords, X1, X2);

parametricGPCoord = mapGlobalToParametric(...
    globalCoord, baseProblem.coords(1), baseProblem.coords(end));

elementEnrichedIndex = element - (overlayXProblem.N - overlayXProblem.XN);

if elementEnrichedIndex == 1
    indexLocalEnrichedNodes = 2;
else
    indexLocalEnrichedNodes = [1, 2];
end

modalDofs = length(indexLocalEnrichedNodes)*overlayXProblem.modes;

projectionOperator = zeros(numberOfProjectionPoints, size(overlayXProblem.LM, 2) + modalDofs);

PODCoefficients = overlayXProblem.reductionOperator(...
    (elementEnrichedIndex-1)*floor(refinedDofs/overlayXProblem.XN)+1:(elementEnrichedIndex-1)*floor(refinedDofs/overlayXProblem.XN)...
    + ceil(refinedDofs/overlayXProblem.XN),:);

for k=1:length(x)
    [ N, ~] = shapeFunctionsAndDerivatives(localCoords(k));
    
    
    [F, ~] = PODModesAndDerivativesMultiscale(localCoords(k), elementGlobalCoords, overlayXProblem.modes,...
        PODCoefficients, integrationDomain, subDomainIndex, indexLocalEnrichedNodes);
    
    projectionOperator(k,1:size(N, 2) + size(F, 2)) =...
        [N, F];
    
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(parametricGPCoord,...
        baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(NIga,2)) = NIga(:);
end

%evaluate the overlay solution
projectedOverlayCoefficients = projectionOperator(:,:) * solutionOverlayCoefficients(overlayXProblem.LMBC(elementEnrichedIndex,1:modalDofs+2));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(:);

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end

