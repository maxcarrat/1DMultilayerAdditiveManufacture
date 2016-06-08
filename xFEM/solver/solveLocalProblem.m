function [ localTemperatureCoefficients ] = solveLocalProblem( localTemperatureCoefficients, problem, timeStepSize, numberOfModes )
%SOLVEGLOBALPROBLEM returns the nodal temperature values of the
%local/enriched problem.
%   problemCoarse = problem struct on the coarse mesh

%% Loop over POD modes 
%The enrichment solution of each mode is added to the final enrichment
%vector

%Assembly the local reduced basis
[M, K, f] = assemblyLocalProblem(problem);

%Apply Weak Dirichlet BCs at the end node of the enriched element
constrainedTemperature = problem.dirichlet_bc(2, 2);
constrainedNode = problem.coords(end);
[K, f] = applyWeakDirichletBCs(problem, constrainedNode, constrainedTemperature, K, f);

%Solve
RHS = timeStepSize * (f - K * localTemperatureCoefficients);
LHS = M + timeStepSize * K;
increment = LHS\RHS;

localTemperatureCoefficients = localTemperatureCoefficients + increment;

end

