function [ M_oo, f_o ] = constrainOverlayMeshL2Projection( M_oo, f_o, problem )
%CONSTRAINOVERLAYMESHL2PROJECTION Constrain the overlay mesh in L2
%projection

for i=1:size(problem.dirichlet_bc,1)
    M_oo(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) = M_oo(problem.dirichlet_bc(i, 1), problem.dirichlet_bc(i, 1)) + problem.penaltyL2;
    f_o(problem.dirichlet_bc(i, 1)) = f_o(problem.dirichlet_bc(i, 1)) + problem.penaltyL2*problem.dirichlet_bc(i, 2);
end

end
