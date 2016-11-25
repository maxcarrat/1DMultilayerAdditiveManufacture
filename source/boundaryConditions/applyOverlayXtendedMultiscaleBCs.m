function [ K, Kprime, KCoupling, M, MCoupling, Mprime, f ] = applyOverlayXtendedMultiscaleBCs( problem, K, Kprime, KCoupling, M, MCoupling, Mprime, f )
%%Xtended MultiscaleBCs 
% apply Dirichlet BCs
%   left Dirichlet BC
%   right Dirichlet BC

for i=1:size(problem.dirichlet_bc,1)
    K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = K(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
    Kprime(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = Kprime(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
%     KCoupling(problem.dirichlet_bc(i, 1), :) = zeros( 1, size(KCoupling, 2));
    
%     M(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = M(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
%     Mprime(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = Mprime(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penalty;
%     MCoupling(problem.dirichlet_bc(i, 1), :) = zeros( 1, size(MCoupling, 2));
    
    f(problem.dirichlet_bc(i, 1)) =  f(problem.dirichlet_bc(i, 1)) + problem.penalty * problem.dirichlet_bc(i, 2);
end

% apply Neumann BCs
f(problem.neumann_bc(1, 1)) =  f(problem.neumann_bc(1, 1)) + problem.neumann_bc(1, 2);

end

