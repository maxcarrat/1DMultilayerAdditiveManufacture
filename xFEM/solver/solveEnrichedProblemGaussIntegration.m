function [ localTemperatureCoefficients ] = solveEnrichedProblemGaussIntegration( localTemperatureCoefficients, problem, time,...
    timeStepSize, integrationOrder )
%SOLVEGLOBALPROBLEMGAUSSINTEGRATION returns the nodal temperature values of the
%local/enriched problem.
%   problem = problem struct on the coarse mesh


%Assembly the local reduced basis
[M, K, f] = assemblyXFEMLinearSystem(problem, time, integrationOrder);

% Apply Weak Dirichlet BCs

for i=1:size(problem.dirichlet_bc,1)
    K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
    f(problem.dirichlet_bc(i, 1)) =  f(problem.dirichlet_bc(i, 1)) + problem.penalty * problem.dirichlet_bc(i, 2);
end

% Apply Neumann BCs
f(problem.neumann_bc(1, 1)) = f(problem.neumann_bc(1, 1)) + problem.neumann_bc(1, 2);

% constrainedNode = problem.coords(end);
% [K, f] = applyWeakDirichletBCs(problem, constrainedNode, constrainedTemperature, K, f);

%Solve
RHS = timeStepSize * (f - K * localTemperatureCoefficients);
LHS = M + timeStepSize * K;
increment = LHS\RHS;

localTemperatureCoefficients = localTemperatureCoefficients + increment;

end

