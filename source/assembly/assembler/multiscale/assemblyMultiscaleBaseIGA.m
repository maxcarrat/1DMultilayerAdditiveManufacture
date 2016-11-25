function [ M_bb, Mprime_bb, K_bb, Kprime_bb, f_b ] = assemblyMultiscaleBaseIGA(overlayProblem, baseProblem, time,...
    integrationOrder, overlaySolutionCoefficients, baseSolutionCoefficients,...
    lastConvergentBaseSolution, lastConvergentOverlaySolution)
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

%% Allocate matrices ======================================================
%overlay K matrix and Kprime matrix
K_bb = zeros(baseProblem.gdof,baseProblem.gdof);
Kprime_bb = zeros(baseProblem.gdof,baseProblem.gdof);
%overlay M matrix Mprime matrix
M_bb = zeros(baseProblem.gdof,baseProblem.gdof);
Mprime_bb = zeros(baseProblem.gdof,baseProblem.gdof);
%overlay rhs vector
f_b = zeros(baseProblem.gdof, 1);

%gauss points base elements
[rGP, wGP] = gaussPoints( integrationOrder );
numberOfIntegrationPoints = length(rGP);

%loop over non-overlapped base elements
for e=1:baseProblem.N-1
    
    ldof_b = baseProblem.p+1;
    
    X1 = baseProblem.knotVector(e + baseProblem.p);
    X2 = baseProblem.knotVector(e + baseProblem.p + 1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        parametricGPCoord = mapParentToLocal(rGP(iGP),X1, X2);
        globalGPCoord = mapParametricToGlobal(parametricGPCoord, baseProblem);
        
        [Nsplines, Bsplines] = BsplinesShapeFunctionsAndDerivatives(parametricGPCoord, baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        inverseJacobianX_Xi = 1 / detJacobianX_Xi;
        
        %% Integrate IGA block
        %extrnal heat source
        f_b(baseProblem.LM(e,1:ldof_b)) = f_b(baseProblem.LM(e,1:ldof_b)) + Nsplines(baseProblem.LM(e,1:ldof_b))' * baseProblem.rhs(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time) * wGP(iGP) * detJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %Capacity matrix
        M_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = M_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            baseProblem.heatCapacity( globalGPCoord, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, lastConvergentBaseSolution))...
            * (Nsplines(baseProblem.LM(e,1:ldof_b))' * Nsplines(baseProblem.LM(e,1:ldof_b)))...
            * wGP(iGP) * detJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %Capacity matrix derivative
        deltaBaseSolution = ( lastConvergentBaseSolution - baseSolutionCoefficients );
        
        Mprime_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = Mprime_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            baseProblem.heatCapacityDerivative( globalGPCoord, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, baseSolutionCoefficients), ...
            evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, lastConvergentBaseSolution))...
            * (Nsplines(baseProblem.LM(e,1:ldof_b))' * Nsplines(baseProblem.LM(e,1:ldof_b))) * (Nsplines(baseProblem.LM(e,1:ldof_b)) * deltaBaseSolution(baseProblem.LM(e,1:ldof_b)))...
            * wGP(iGP) * detJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %Conductivity matrix
        K_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = K_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            baseProblem.k(globalGPCoord,...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, baseSolutionCoefficients))...
            * (Bsplines(baseProblem.LM(e,1:ldof_b))' * Bsplines(baseProblem.LM(e,1:ldof_b)))...
            * wGP(iGP) * inverseJacobianX_Xi * baseProblem.F_map(X1,X2);
        
        %Conductivity matrix derivative
        Kprime_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) = Kprime_bb(baseProblem.LM(e, 1:ldof_b), baseProblem.LM(e, 1:ldof_b)) +...
            baseProblem.kDerivative(globalGPCoord,...
            time, evaluateTemperature(e, rGP(iGP), overlayProblem, baseProblem, baseSolutionCoefficients))...
            * (Bsplines(baseProblem.LM(e,1:ldof_b))' * Nsplines(baseProblem.LM(e,1:ldof_b)))...
            * (inverseJacobianX_Xi *Bsplines(baseProblem.LM(e,1:ldof_b)) * baseSolutionCoefficients(baseProblem.LM(e,1:ldof_b))) * wGP(iGP) * baseProblem.F_map(X1,X2);
        
    end
    
end

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
    ldof_o = 2;
    
    %overlay mesh element boundaries
    X1 = overlayProblem.coords(e);
    X2 = overlayProblem.coords(e+1);

    % Gauss integration
    for iGP = 1:numberOfComposedIntegrationPoints
        
        [N, B] = shapeFunctionsAndDerivatives(rGPComposed(iGP));
        globalGPCoord = mapLocalToGlobal(rGPComposed(iGP),X1, X2);
        
        [Nsplines, Bsplines] = BsplinesShapeFunctionsAndDerivatives(mapGlobalToParametric(...
            globalGPCoord, baseProblem.coords(1), baseProblem.coords(end)), baseProblem.p, baseProblem.knotVector);
        
        detJacobianX_Xi = baseProblem.coords(end) - baseProblem.coords(1);
        inverseJacobianX_Xi = 1 / detJacobianX_Xi;
        
        Bsplines = Bsplines * inverseJacobianX_Xi;
        B = B *  overlayProblem.B_map(X1,X2);
        
        %% Integrate IGA block
        %external heat source
        f_b(baseProblem.LM(end,1:ldof_b)) = f_b(baseProblem.LM(end,1:ldof_b)) + Nsplines(baseProblem.LM(end,1:ldof_b))' * baseProblem.rhs(mapLocalToGlobal(rGPComposed(iGP), X1, X2),...
            time) * wGPComposed(iGP) * baseProblem.F_map(X1,X2);
        
        %Capacity matrix
        M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) =...
            M_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
            baseProblem.heatCapacity( globalGPCoord, evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
            evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, lastConvergentOverlaySolution, lastConvergentBaseSolution))...
            * (Nsplines(baseProblem.LM(end,1:ldof_b))' * Nsplines(baseProblem.LM(end,1:ldof_b)))...
            * wGPComposed(iGP) * baseProblem.F_map(X1,X2);
        
        %Capacity matrix derivative
        deltaBaseSolution = lastConvergentBaseSolution - baseSolutionCoefficients;
        deltaOverlaySolution = -(lastConvergentOverlaySolution - overlaySolutionCoefficients);
        
        Mprime_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) =...
            Mprime_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
            baseProblem.heatCapacityDerivative( globalGPCoord, evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients), ...
             evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, lastConvergentOverlaySolution, lastConvergentBaseSolution))...
            * (Nsplines(baseProblem.LM(end,1:ldof_b))' * Nsplines(baseProblem.LM(end,1:ldof_b))) * (Nsplines * deltaBaseSolution +...
            N * deltaOverlaySolution(overlayProblem.LM(e, 1:ldof_o))) * wGPComposed(iGP) * baseProblem.F_map(X1,X2);
        
        %Diffusion matrix
        K_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) =...
            K_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
            baseProblem.k(globalGPCoord,...
            time, evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (Bsplines(baseProblem.LM(end,1:ldof_b))' * Bsplines(baseProblem.LM(end,1:ldof_b))) * wGPComposed(iGP) * baseProblem.F_map(X1,X2);
        
        %Diffusion matrix derivative
        Kprime_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) =...
            Kprime_bb(baseProblem.LM(end, 1:ldof_b), baseProblem.LM(end, 1:ldof_b)) +...
            baseProblem.kDerivative(globalGPCoord,...
            time, evaluateTemperatureComposed(e, rGPComposed(iGP), overlayProblem, baseProblem, overlaySolutionCoefficients, baseSolutionCoefficients))...
            * (Bsplines(baseProblem.LM(end,1:ldof_b))' * Nsplines(baseProblem.LM(end,1:ldof_b))) * (Bsplines * baseSolutionCoefficients  +...
            B * overlaySolutionCoefficients(overlayProblem.LM(e, 1:ldof_o))) * wGPComposed(iGP) * baseProblem.F_map(X1,X2);       
    end    
end

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
    [NIga, ~] =  BsplinesShapeFunctionsAndDerivatives(localCoord(k), baseProblem.p, baseProblem.knotVector);
    projectionOperatorSpline(k,1:size(Nspline,2)) = NIga(baseProblem.LM(e,:));
end

%evaluate the base solution
projectedBaseCoefficients = projectionOperatorSpline(:,:) * solutionBaseCoefficients(baseProblem.LM(e,:));

projectedCoefficients = projectedBaseCoefficients;

end

%% Evaluate temperature for composed integration
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

localCoords = x;

globalCoord = mapLocalToGlobal(x,overlayProblem.coords(e), overlayProblem.coords(e + 1));

parametricGPCoord = mapGlobalToParametric(...
    globalCoord, baseProblem.coords(1), baseProblem.coords(end));

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