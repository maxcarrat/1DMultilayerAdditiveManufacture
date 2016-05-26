function [M, K, f] = assemblyLocalProblem(problem, iMode, numberOfModesSupports)
%   [M, K, f] = ASSEMBLYLOCALPROBLEM(problem) assembles the mass and the conductivity matrix and load vector 
%   problem = definition of the boundary value problem
%   iMode = ith POD mode


%conductivity matrix of the enriched problem
K = zeros(numberOfModesSupports+1,numberOfModesSupports+1);
%capacity matrix of the enriched problem
M = zeros(numberOfModesSupports+1,numberOfModesSupports+1);
%load vector of the enriched problem
f = zeros(numberOfModesSupports+1, 1);

enrichedElementCoords = linspace(problem.N-numberOfModesSupports, problem.N, 2^problem.refinementDepth);

KE = rbLocalConductivityMatrix(problem, iMode, enrichedElementCoords);
ME = rbLocalCapacityMatrix(problem, iMode, enrichedElementCoords);

for modesSupport=1:numberOfModesSupports
    
    ldof = 2;
    X1 = problem.coords(problem.N-2);
    X2 = problem.coords(problem.N-1);
    
    fel = rbLocalLoadVector(problem, iMode, enrichedElementCoords, problem.N);
        
    f(problem.rbLM(modesSupport,1:ldof)) = f(problem.rbLM(modesSupport,1:ldof)) + problem.F_map(X1,X2) * fel;
    M(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof)) = M(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof))...
        + problem.F_map(X1,X2) * ME(1:ldof, 1:ldof);
    K(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof)) = K(problem.rbLM(modesSupport, 1:ldof), problem.rbLM(modesSupport, 1:ldof))...
        + problem.B_map(X1,X2) * KE(1:ldof, 1:ldof);
    
end
  
end