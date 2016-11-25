function [M, K, f] = assemblyMultiPhaseSystemIGA(problem, time, integrationOrder,...
    solutionCoefficients, oldSolutionCoefficients)
%   ASSEMBLYMULTIPHASESYSTEMX assembles the mass and the conductivity matrix and load vector
%   problem = definition of the boundary value problem
%   time = current time
%   integrationOrder = number of integration points

%% Allocate matrices
%global conductivity matrix
K = zeros(problem.gdof,problem.gdof);
%global capacity matrix
M = zeros(problem.gdof,problem.gdof);
%global load vector
f = zeros(problem.gdof, 1);

%gauss points
[rGP, wGP] = gaussPoints( integrationOrder );

numberOfIntegrationPoints = length(rGP);

for e=1:problem.N
    
    ldof = problem.p+1;
    
    Xi1 = problem.knotVector(e + problem.p);
    Xi2 = problem.knotVector(e + problem.p + 1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        localCoords = mapParentToLocal(rGP(iGP), Xi1, Xi2);
        globalCoords = mapParametricToGlobal(localCoords, problem);
        [N, B] = BsplinesShapeFunctionsAndDerivatives(localCoords,problem.p, problem.knotVector);
        
        JacobianX_Xi = B(problem.LM(e, 1:ldof)) * problem.coords(problem.LM(e, 1:ldof))';

%         detJacobianX_Xi = norm(JacobianX_Xi);
        detJacobianX_Xi = problem.coords(end) - problem.coords(1);

        inverseJacobianX_Xi = 1 / detJacobianX_Xi;
        
        %% Integrate FEM block
        %external heat source
        f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + N(problem.LM(e, 1:ldof))' * problem.rhs( globalCoords,...
            time) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
        
        %Capacity matrix
        M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
            problem.heatCapacity( rGP(iGP), evaluateTemperature(e, localCoords, problem, solutionCoefficients), ...
            evaluateTemperature(e, localCoords, problem, oldSolutionCoefficients))...
            * (N(problem.LM(e, 1:ldof))' * N(problem.LM(e, 1:ldof))) * wGP(iGP) * detJacobianX_Xi * problem.F_map(Xi1, Xi2);
        
        %Diffusion matrix
        K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
            problem.k(globalCoords, time, evaluateTemperature(e, localCoords, problem, solutionCoefficients))...
            * (B(problem.LM(e, 1:ldof))' * B(problem.LM(e, 1:ldof))) * wGP(iGP) * inverseJacobianX_Xi * problem.F_map(Xi1, Xi2);
        
    end
    
end

end

function [ projectedCoefficients ] = evaluateTemperature(e, x, problem, solutionCoefficients)
% EVALUATETEMPERATURE project the previous solution onto the element.
%   e = element index
%   x = post-processing mesh
%   problem
%   solutionCoefficients = temeprature distribution of the previous mesh
%   modes = number of enrichment modes
%   derivative = order of derivatives

numberOfProjectionPoints = length(x);
projectionOperator = zeros(numberOfProjectionPoints, size(problem.LM, 2));

m=length(problem.knotVector);

N = zeros(length(x), m - 1 - problem.p);
B = zeros(length(x), m - 1 - problem.p);

localCoords = x;

for k=1:length(x)
    [N(k,:), B(k,:)] = BsplinesShapeFunctionsAndDerivatives(localCoords, problem.p, problem.knotVector);
end

for i=1:length(x)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end


projectedCoefficients = projectionOperator(problem.LM(e,:)) * solutionCoefficients(problem.LM(e,:));

end