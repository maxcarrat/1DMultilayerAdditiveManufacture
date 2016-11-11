function [ M, f ] = constrainXtendedMeshL2Projection( M, f, problem )
%CONSTRAINXTENDEDMESHL2PROJECTION  Constrain the Xtended mesh in L2
%projection

for i=1:size(problem.Xdirichlet_bc, 1)
    M(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) = M(problem.Xdirichlet_bc(i, 1), problem.Xdirichlet_bc(i, 1)) + problem.penalty;
    f(problem.Xdirichlet_bc(i, 1)) =  f(problem.Xdirichlet_bc(i, 1))*problem.Xdirichlet_bc(i, 2);
end

end

