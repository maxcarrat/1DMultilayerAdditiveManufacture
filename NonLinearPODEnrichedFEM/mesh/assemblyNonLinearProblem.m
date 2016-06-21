function    [M, K, f] = assemblyNonLinearProblem(problem, solution)
%ASSEMBLYNONLINEARPROBLEM returns the linear system of the non-linear
%problem
%   problem = struct of the non-linear transient heat problem
%   solution = actula temperature values

%global conductivity matrix
K = zeros(problem.gdof,problem.gdof);
%global capacity matrix
M = zeros(problem.gdof,problem.gdof);
%load vector
f = zeros(problem.gdof, 1);

ME = localCapacityMatrix(problem);

for e=1:problem.N
    
    ldof = 2;
    X1 = problem.coords(e);
    X2 = problem.coords(e+1);
    
    KE = nonLinearConductivityMatrix(e, problem, solution);
    fe = localLoadVector(e, problem);
    
    f(problem.LM(e,1:ldof)) = f(problem.LM(e,1:ldof)) + problem.F_map(X1,X2) * fe;
    M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = M(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
    K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) = K(problem.LM(e, 1:ldof), problem.LM(e, 1:ldof)) + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
    
end

end

