function [ M_oo, Mprime_oo, K_oo, Kprime_oo, f_o ] = assemblyMultiscaleOverlayFEM(overlayProblem, baseProblem, time,...
    integrationOrder, overlaySolutionCoefficients, baseSolutionCoefficients, lastConvergedBaseSolution, lastConvergentOverlaySolution)
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

%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

for e=1:overlayProblem.N
    
    ldof_b = baseProblem.p + 1;
    
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);
        

        [Nsplines, Bsplines] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);

        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        inverseJacobianX_Xi = 1 / detJacobianX_Xi;

        
        %% Integrate FEM block
        %external heat source
        f_o(overlayProblem.LM(e,:)) = f_o(overlayProblem.LM(e,:)) + N' * overlayProblem.rhs(globalGPCoord,...
            time) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Capacity matrix
        M_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) = M_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
            overlayProblem.heatCapacity( globalGPCoord, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, lastConvergentOverlaySolution, lastConvergedBaseSolution))...
            * (N' * N) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Capacity matrix derivative
        deltaBaseSolution = -(lastConvergedBaseSolution - baseSolutionCoefficients);
        deltaOverlaySolution = (lastConvergentOverlaySolution - overlaySolutionCoefficients);
        
        Mprime_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) = Mprime_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
            overlayProblem.heatCapacityDerivative( globalGPCoord, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem,overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem,lastConvergentOverlaySolution, lastConvergedBaseSolution))...
            * (N' * N) * (Nsplines(baseProblem.LM(end, 1:ldof_b)) * deltaBaseSolution(baseProblem.LM(end, 1:ldof_b))...
            + N * deltaOverlaySolution(overlayProblem.LM(e,:))) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Conductivity matrix
        K_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) = K_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
            overlayProblem.k(globalGPCoord,...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (B' * B) * wGP(iGP) * overlayProblem.B_map(X1,X2);
        
        %Conductivity matrix derivative
        Kprime_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) = Kprime_oo(overlayProblem.LM(e,:), overlayProblem.LM(e,:)) +...
            overlayProblem.kDerivative(globalGPCoord,...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (B' * N) * (Bsplines(baseProblem.LM(end, 1:ldof_b)) * baseSolutionCoefficients(baseProblem.LM(end, 1:ldof_b)) * inverseJacobianX_Xi...
            + B * overlaySolutionCoefficients(overlayProblem.LM(e,:)) * overlayProblem.B_map(X1,X2)) * wGP(iGP);
        
    end
end

end

function [ projectedCoefficients ] = evaluateTemperature(e, x, overlayProblem,...
    baseProblem, solutionOverlayCoefficients, solutionBaseCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = evaluation local coordinates
%   overlayProblem = problem struct of the overlay mesh
%   baseProblem = problem struct of the base mesh
%   solutionOverlayCoefficients = temeprature distribution of the previous overlay mesh
%   solutionBaseCoefficients = temeprature distribution of the previous base mesh


numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, size(overlayProblem.LM, 2));
projectionOperatorSpline = zeros(numberOfProjectionPoints, size(baseProblem.LM, 2));

N = zeros(length(x), 2);
X1 = overlayProblem.coords(e);
X2 = overlayProblem.coords(e+1);

localCoords = x;
globalGPCoord = mapLocalToGlobal(localCoords,X1, X2);

parametricGPCoord = mapGlobalToParametric(...
    globalGPCoord, baseProblem.coords(1), baseProblem.coords(end));

for k=1:length(x)
    [projectionOperator(k,1:size(N,2)), ~] = shapeFunctionsAndDerivatives(localCoords(k));
    

    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(parametricGPCoord, baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(NIga,2)) = NIga(:);
end

%evaluate the overlay solution
projectedOverlayCoefficients = projectionOperator(:,:) * solutionOverlayCoefficients(overlayProblem.LM(e,:));

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(:);

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end