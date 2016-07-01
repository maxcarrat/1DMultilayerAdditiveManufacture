function [ M, K, f ] = applyGlobalBCs( problem, M, K, f )
%assemblyGlobalBCs apply Dirichlet BCs
%   left Dirichlet BC
%   right Dirichlet BC

    for i=1:size(problem.dirichlet_bc,1)
        K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
        f(problem.dirichlet_bc(i, 1)) =  f(problem.dirichlet_bc(i, 1)) + problem.penalty * problem.dirichlet_bc(i, 2);
    end
    
% apply Neumann BCs
        f(problem.neumann_bc(1, 1)) =  f(problem.neumann_bc(1, 1)) + problem.neumann_bc(1, 2);

end
