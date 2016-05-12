function ME = localCapacityMatrix(problem)
%   ME = LOCALCAPACITYMATRIX(problem) computes the load vector local to one element
%   problem = definition of the boundary value problem. see also HP_FEM
 
    ldof = 2;
    
    ME = zeros(ldof, ldof);
    
    for i=1:ldof
        for j=1:ldof
            ME(i,j) = problem.M(@(x, d)problem.basis_fun(x, i, d), @(x, d)problem.basis_fun(x, j, d));
        end
    end

end
