function [ K_bb, Kprime_bb, M_bb, Mprime_bb, f_ ] = assemblyMultiscaleBaseIGA(overlayProblem, baseProblem, time,...
    integrationOrder, overlaySolutionCoefficients, baseSolutionCoefficients, lastConvergedBaseSolution, lastConvergentOverlaySolution)
%ASSEMBLYMULTISCALEBASEIGA  assembles the mass, the conductivity matrix
%and load vector of the base IGA mesh.
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
%K_bb = base conductivity matrix
%Kprime_bb = derivative of the diffussion matrix
%M_bb = base capacity matrix
%Mprime_bb = derivative of the capacity matrix
%f_b = rhs base vector

%% Allocate matrices
%overlay K matrix and Kprime matrix
K_bb = zeros(baseProblem.gdof,baseProblem.gdof);
Kprime_bb = zeros(baseProblem.gdof,baseProblem.gdof);
%overlay M matrix Mprime matrix
M_bb = zeros(baseProblem.gdof,baseProblem.gdof);
Mprime_bb = zeros(baseProblem.gdof,baseProblem.gdof);
%overlay rhs vector
f_ = zeros(baseProblem.gdof, 1);

%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

for e=1:overlayProblem.N
    
    ldof = 2;
    
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
        globalGPCoord = mapLocalToGlobal(rGP(iGP),X1, X2);

        [Nsplines, Bspline] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, overlayProblem.coords(1), overlayProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        
        %% Integrate FEM block
        %extrnal heat source
        f_(overlayProblem.LM(e,1:ldof)) = f_(overlayProblem.LM(e,1:ldof)) + N' * overlayProblem.rhs(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Capacity matrix
        M_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) = M_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) +...
            overlayProblem.heatCapacity( rGP(iGP), evaluateTemperature(e, rGP(iGP), overlayProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), lastConvergentOverlaySolution, lastConvergedBaseSolution))...
            * (N' * N) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Capacity matrix derivative
        deltaBaseSolution = lastConvergedBaseSolution - baseSolutionCoefficients;
        deltaOverlaySolution = lastConvergentOverlaySolution - overlaySolutionCoefficients;
        
        Mprime_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) = Mprime_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) +...
            overlayProblem.heatCapacityDerivative( rGP(iGP), evaluateTemperature(e, rGP(iGP), overlayProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, lastConvergentOverlaySolution, lastConvergedBaseSolution))...
            * (N' * N) * (Nsplines * deltaBaseSolution + N * deltaOverlaySolution) * wGP(iGP) * overlayProblem.F_map(X1,X2);
        
        %Diffusion matrix
        K_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) = K_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) +...
            overlayProblem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (B' * B) * wGP(iGP) * overlayProblem.B_map(X1,X2);
        
        %Diffusion matrix derivative
        Kprime_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) = Kprime_bb(overlayProblem.LM(e, 1:ldof), overlayProblem.LM(e, 1:ldof)) +...
            overlayProblem.kderivative(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (B' * N) * (Bspline * baseSolutionCoefficients + B * overlaySolutionCoefficients) * wGP(iGP) * overlayProblem.B_map(X1,X2);
        
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
Nspline = zeros(length(x), 2);

localCoords = x;
globalCoord = mapLocalToGlobal(x,overlayProblem.coords(e), overlayProblem.coords(e+1));

for k=1:length(x)
    [N(k,:), ~] = shapeFunctionsAndDerivatives(localCoords(k));
    [Nspline(k,:), ~] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric...
        (globalCoord, overlayProblem.coords(1), overlayProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
end

%evaluate the overlay solution
for i=1:length(x)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end
projectedOverlayCoefficients = projectionOperator * solutionOverlayCoefficients(overlayProblem.LM(e,:));

%evaluate the base solution
for i=1:length(x)
    projectionOperatorSpline(i,1:size(N,2)) = N(i,:);
end
projectedBaseCoefficients = projectionOperatorSpline * solutionBaseCoefficients(baseProblem.LM(e,:));

projectedCoefficients = projectedOverlayCoefficients + projectedBaseCoefficients;

end

end

