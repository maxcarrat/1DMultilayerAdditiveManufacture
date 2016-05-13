function KE = localConductivityMatrix(problem)
%   KE = LOCALCONDUCTIVITYMATRIX(problem, coefficients) computes the conductivity local to one element
%   problem = definition of the boundary value problem.

    ldof = 2;
    
    KE = zeros(ldof, ldof);
    
    for i=1:ldof
        for j=1:ldof
            KE(i,j) = problem.B(@(x, d)problem.basis_fun(x, i, d), @(x, d)problem.basis_fun(x, j, d));
        end
    end

end
