function [M, K, f] = assemblyNonLinearSystem(problem, time, integrationOrder, solutionCoefficients)
%   ASSEMBLYNONLINEARSYSTEMXFEM assembles the mass and the conductivity matrix and load vector
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
    
    ldof = 2;
    
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    % Gauss integration
    for iGP = 1:numberOfIntegrationPoints
        
        [N, B] = shapeFunctionsAndDerivatives(rGP(iGP));
        
        %% Integrate FEM block
        %extrnal heat source
        f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + N' * problem.rhs(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time) * wGP(iGP) * problem.F_map(X1,X2);
        
        %Capacity matrix
        M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
            problem.heatCapacity * (N' * N) * wGP(iGP) * problem.F_map(X1,X2);
        
        %Diffusion matrix
        K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) +...
            problem.B_map(X1,X2) * problem.k(mapLocalToGlobal(rGP(iGP), X1, X2),...
            time, evaluateTemperature(e, rGP(iGP), problem, solutionCoefficients))...
            * (B' * B) * wGP(iGP);
        
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

N = zeros(length(x), 2);
B = zeros(length(x), 2);

localCoords = x;

for k=1:length(x)
    [N(k,:), B(k,:)] = shapeFunctionsAndDerivatives(localCoords(k));
end

for i=1:length(x)
    projectionOperator(i,1:size(N,2)) = N(i,:);
end


projectedCoefficients = projectionOperator * solutionCoefficients(problem.LM(e,:));

end