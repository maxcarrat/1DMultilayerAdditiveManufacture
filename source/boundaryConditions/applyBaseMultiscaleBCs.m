function [ K, Kprime, KCoupling, MCoupling, f ] = applyBaseMultiscaleBCs( baseProblem, KCoupling, K, Kprime, MCoupling, f )
%assemblyBaseMultiscaleBCs apply Dirichlet BCs
%   left Dirichlet BC
%   right Dirichlet BC

for i=1:size(baseProblem.dirichlet_bc,1)
    K(baseProblem.dirichlet_bc(i, 1), baseProblem.dirichlet_bc(i, 1)) = K(baseProblem.dirichlet_bc(i, 1), baseProblem.dirichlet_bc(i, 1)) + baseProblem.penalty;
    Kprime(baseProblem.dirichlet_bc(i, 1), baseProblem.dirichlet_bc(i, 1)) = Kprime(baseProblem.dirichlet_bc(i, 1), baseProblem.dirichlet_bc(i, 1)) + baseProblem.penalty;
    KCoupling(baseProblem.dirichlet_bc(i, 1), :) = zeros( 1, size(KCoupling, 2));

%     MCoupling(baseProblem.dirichlet_bc(i, 1), :) = zeros( 1, size(MCoupling, 2));
    
    f(baseProblem.dirichlet_bc(i, 1)) =  f(baseProblem.dirichlet_bc(i, 1)) + baseProblem.penalty * baseProblem.dirichlet_bc(i, 2);
end

% apply Neumann BCs
f(baseProblem.neumann_bc(1, 1)) =  f(baseProblem.neumann_bc(1, 1)) + baseProblem.neumann_bc(1, 2);

end
