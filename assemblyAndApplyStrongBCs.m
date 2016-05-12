function [ M, K, f ] = assemblyAndApplyStrongBCs( problem )
%APPLYBCS apply Dirichlet BCs
%   left Dirichlet BC
%   right Dirichlet BC

[M, K, f] = assembly(problem);

    for i=1:size(problem.dirichlet_bc,1)
        M(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = M(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
        K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
        f(problem.dirichlet_bc(i, 1)) =  f(problem.dirichlet_bc(i, 1)) + problem.penalty * problem.dirichlet_bc(i, 2);
    end

end

