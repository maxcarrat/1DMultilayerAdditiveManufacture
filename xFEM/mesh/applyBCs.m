function [ LHS, RHS ] = applyBCs( M, K, f, problem, previousSolution, timeStepSize )
%APPLYBCS apply Dirichlet BCs
%   left Dirichlet BC
%   right Dirichlet BC

for i=1:size(problem.dirichlet_bc,1)
    K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
    f(problem.dirichlet_bc(i, 1)) =  f(problem.dirichlet_bc(i, 1)) + problem.penalty * problem.dirichlet_bc(i, 2);
end

    RHS = timeStepSize*(f - K*previousSolution);
    LHS = M + timeStepSize*K;
end