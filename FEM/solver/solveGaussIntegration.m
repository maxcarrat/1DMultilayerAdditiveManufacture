function [ localTemperatureCoefficients ] = solveGaussIntegration( localTemperatureCoefficients,...
    problem, time, timeStepSize, integrationOrder )
%SOLVEGAUSSINTEGRATION returns the nodal temperature values of the
%coarse/global problem using Backward Euler implictit scheme
%   problemCoarse = problem struct on the coarse mesh
%   activeMesh = active elements at teh current time step



%assembly and apply BCs
[M, K, f] = assemblyFastIntegration(problem, time, integrationOrder);
% [M, K, f] = assemblyXFEMGaussIntegrationLinearSystem(problem, time, integrationOrder);
% [M, K, f] = assemblyXFEMLinearSystem(problem, time, integrationOrder);

[M, K, f] = applyGlobalBCs(problem, M, K, f);

% Apply Neumann BCs
f(problem.neumann_bc(1, 1)) = f(problem.neumann_bc(1, 1)) + problem.neumann_bc(1, 2);

%Solve
RHS = timeStepSize * (f - K * localTemperatureCoefficients);
LHS = M + timeStepSize * K;
increment = LHS\RHS;

localTemperatureCoefficients = localTemperatureCoefficients + increment;
   
end

