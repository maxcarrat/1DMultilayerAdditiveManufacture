function [ M, K, f ] = assemblyLocalNonLinearProblem(problem, solution, iMode, numberOfModesSupports)
%ASSEMBLYLOCALNONLINEARPROBLEM  returns the linear system of the non-linear
%problem
%   problem = struct of the non-linear transient heat problem
%   solution = actula temperature values
%   iMode = ith POD mode index
%   numberOfModesSupport = number of supports of the modes

%conductivity matrix of the enriched problem
K = zeros(numberOfModesSupports+1,numberOfModesSupports+1);
%capacity matrix of the enriched problem
M = zeros(numberOfModesSupports+1,numberOfModesSupports+1);
%load vector of the enriched problem
f = zeros(numberOfModesSupports+1, 1);

enrichedElementCoords = linspace(problem.N-numberOfModesSupports, problem.N, 2^problem.refinementDepth);

ME = rbLocalCapacityMatrix(problem, iMode, enrichedElementCoords);
lastActiveElement = problem.N;

for modesSupport=1:numberOfModesSupports
    
    ldof = 2;
    X1 = problem.coords(lastActiveElement-2);
    X2 = problem.coords(lastActiveElement-1);
    
    KE = rbNonLinearConductivityMatrix(problem, solution, iMode, enrichedElementCoords, lastActiveElement);
    fe = rbLocalLoadVector(problem, iMode, enrichedElementCoords, lastActiveElement);
    
    f(problem.rbLM(modesSupport,1:ldof)) = f(problem.rbLM(modesSupport,1:ldof)) + problem.F_map(X1,X2) * fe;
    M(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof)) = M(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof))...
        + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
    K(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof)) = K(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof))...
        + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
end

end

