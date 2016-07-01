function [ localTemperatureCoefficients ] = solveEnrichedProblemPreIntegration( localTemperatureCoefficients, problem, time,...
    timeStepSize, integrationOrder, K_Coupling, K_XFEM, M_Coupling, M_XFEM, f_XFEM, modalDofs )
%SOLVEGLOBALPROBLEMGAUSSINTEGRATION returns the nodal temperature values of the
%local/enriched problem.
%   problem = problem struct on the coarse mesh


%Assembly the local reduced basis
[M_FEM, K_FEM, f_FEM] = assemblyLinearSystem(problem, time, integrationOrder);

K = [K_FEM, K_Coupling(1:modalDofs,:)'; K_Coupling(1:modalDofs,:), K_XFEM(1:modalDofs,1:modalDofs)];
M = [M_FEM, M_Coupling(1:modalDofs,:)'; M_Coupling(1:modalDofs,:), M_XFEM(1:modalDofs,1:modalDofs)];
f = [f_FEM; f_XFEM(1:modalDofs)];

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

