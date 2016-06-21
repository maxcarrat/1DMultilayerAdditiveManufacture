function [M, K, f] = assemblyXFEM(problem)
%   [M, K, f] = ASSEMBLYXFEM(problem) assembles the mass and the conductivity matrix and load vector 
%   problem = definition of the boundary value problem

    %global conductivity matrix
    K = zeros(problem.gdof,problem.gdof);
    %global capacity matrix
    M = zeros(problem.gdof,problem.gdof);
    %load vector
    f = zeros(problem.gdof, 1);

    KE = rbLocalConductivityMatrix(problem, problem.coords);
    ME = rbLocalCapacityMatrix(problem, problem.coords);
    
    for e=1:problem.N
        
        ldof = 2 + 2 * problem.modes;
        X1 = problem.coords(e);
        X2 = problem.coords(e+1);
        
        fe = rbLocalLoadVector(problem, problem.coords, e);

        f(problem.rbLM(e,1:ldof)) = f(problem.rbLM(e,1:ldof)) + problem.F_map(X1,X2) * fe;
        M(problem.rbLM(e, 1:ldof), problem.rbLM(e, 1:ldof)) = M(problem.rbLM(e, 1:ldof), problem.rbLM(e, 1:ldof)) + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
        K(problem.rbLM(e, 1:ldof), problem.rbLM(e, 1:ldof)) = K(problem.rbLM(e, 1:ldof), problem.rbLM(e, 1:ldof)) + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
   
    end
end