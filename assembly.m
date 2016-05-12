function [M, K, f] = assembly(problem)
%   [M, K, f] = ASSEMBLY(problem) assembles the mass and the conductivity matrix and load vector 
%   problem = definition of the boundary value problem

    %global conductivity matrix
    K = zeros(problem.gdof,problem.gdof);
    %global capacity matrix
    M = zeros(problem.gdof,problem.gdof);
    %load vector
    f = zeros(problem.gdof, 1);

    KE = localConductivityMatrix(problem);
    ME = localCapacityMatrix(problem);
    
    for e=1:problem.N
        ldof = 2;
        X1 = problem.coords(e);
        X2 = problem.coords(e+1);
        
        fe = localLoadVector(e, problem);

        f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + problem.F_map(X1,X2) * fe;
        M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
        K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
    end
end