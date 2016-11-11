function [ M_bo, K_bo] = assemblyMultiscaleCouplingBaseFEM(overlayProblem, baseProblem, time,...
    integrationOrder, overlaySolutionCoefficients, baseSolutionCoefficients,...
    lastConvergedBaseSolution, lastConvergentOverlaySolution)
%ASSEMBLYMULTISCALEOVERLAYFEM  assembles the mass, the conductivity matrix
%and load vector
%Input:
%overlayProblem = problem struct of the overlay mesh
%baseProblem = problem struct of the base mesh
%time = current time
%integrationOrder = number of quadrature points for the integration
%overlaySolutionCoefficients = coefficients of the overlay solution
%baseSolutionCoefficients = coefficients of the base solution
%lastConvergedBaseSolution = last convergent base solution
%lastConvergentOverlaySolution = last convergent overlay solution
%Output:
%K_ob = coupling conductivity matrix
%M_ob = overlay capacity matrix

%% Allocate matrices
%coupling K matrix and Kprime matrix
K_bo = zeros(baseProblem.gdof, overlayProblem.gdof);
%coupling M matrix
M_bo = zeros(baseProblem.gdof, overlayProblem.gdof);

%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

for e=1:overlayProblem.N
    
    ldof_b = baseProblem.p + 1;
    
    %overlay mesh element boundaries
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);
        
        [Nspline, Bspline] = BsplinesShapeFunctionsAndDerivatives( mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)),...
            baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        invDetJacobianX_Xi = 1 / detJacobianX_Xi;
        
        %% Integrate Coupling block using composed integration
        %Composed integration allows to integrate discontinuous functions over a
        %domain. The quadrature points are distributed over the smaller sub-domain
        %where the function is still continuous. This technique is essential to
        %correctly integrate the basis function of the overlapping mesh.
        
        %Capacity matrix
        M_bo(baseProblem.LM(end,:), overlayProblem.LM(e,:)) = M_bo(baseProblem.LM(end,:), overlayProblem.LM(e,:)) +...
            overlayProblem.heatCapacity(globalGPCoord,...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, lastConvergentOverlaySolution, lastConvergedBaseSolution))...
                    * (Nspline(baseProblem.LM(end,1:ldof_b))' * N) * wGP(iGP) * overlayProblem.F_map(X1,X2);

        %Conductivity matrix
        K_bo(baseProblem.LM(end,:), overlayProblem.LM(e,:)) = K_bo(baseProblem.LM(end,:), overlayProblem.LM(e,:)) +...
            overlayProblem.k(globalGPCoord,...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
                    * (Bspline(baseProblem.LM(end,1:ldof_b))' * B) * wGP(iGP) * invDetJacobianX_Xi;

    end
    
end

end

function [ projectedCoefficients ] = evaluateTemperature(e, x, overlayProblem,...
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

localCoords = x;

X1 = overlayProblem.coords(e);
X2 = overlayProblem.coords(e+1);

globalGPCoord = mapLocalToGlobal(x,X1, X2);

parametricGPCoord = mapGlobalToParametric(...
    globalGPCoord, baseProblem.coords(1), baseProblem.coords(end));

for k=1:length(x)
    [projectionOperator(k,1:size(N,2)), ~] = shapeFunctionsAndDerivatives(localCoords(k));
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(parametricGPCoord(k), baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(NIga,2)) = NIga(:);
end

%evaluate the overlay solution
projectedOverlayCoefficients = projectionOperator(:,:) * solutionOverlayCoefficients(overlayProblem.LM(e,:));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(:);

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end