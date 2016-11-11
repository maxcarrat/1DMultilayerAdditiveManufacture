function [ K, Kprime, f ] = constrainXtendedMesh( K, Kprime, f,...
    problem )
%CONSTRAINXTENDEDMESH Constrain the Xtended mesh 

for i=1:max(max(problem.LME))
    K(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) =...
        K(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) + problem.penalty;
    
    Kprime(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) =...
        Kprime(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) + problem.penalty;
    
    
    f(problem.Xdirichlet_bc(i, 1)) =  f(problem.Xdirichlet_bc(i, 1)) * ...
        problem.Xdirichlet_bc(i, 2);
end

end


