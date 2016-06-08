function [Menr, Kenr, fenr] = assemblyLocalProblem(problem)
%   [Menr, Kenr, fenr] = ASSEMBLYLOCALPROBLEM(problem) assembles the mass and the conductivity matrix and load vector 
%   of the enriched problem at the last active element of the mesh.
%   problem = definition of the boundary value problem
%   iMode = ith POD mode
%   numberOfModes = number of POD basis 


enrichedElementCoords = linspace(problem.coords(end-1), problem.coords(end),...
    2^problem.refinementDepth);

%conductivity matrix of the locally enriched problem
Kenr = rbLocalConductivityMatrix(problem, enrichedElementCoords)*problem.B_map(problem.coords(end-1), problem.coords(end));
%capacity matrix of the locally enriched problem
Menr = rbLocalCapacityMatrix(problem, enrichedElementCoords)*problem.F_map(problem.coords(end-1), problem.coords(end));
%load vector of the locally enriched problem
fenr = rbLocalLoadVector(problem, enrichedElementCoords, problem.N)*problem.F_map(problem.coords(end-1), problem.coords(end));
  
end