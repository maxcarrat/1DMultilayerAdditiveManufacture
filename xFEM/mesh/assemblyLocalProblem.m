function [Menr, Kenr, fenr] = assemblyLocalProblem(problem)
%   [Menr, Kenr, fenr] = ASSEMBLYLOCALPROBLEM(problem) assembles the mass and the conductivity matrix and load vector 
%   of the enriched problem at the last active element of the mesh.
%   problem = definition of the boundary value problem
%   iMode = ith POD mode
%   numberOfModes = number of POD basis 


enrichedElementCoords = linspace(problem.coords(end-1), problem.coords(end),...
    2^problem.refinementDepth+1);

detJenr = 2/(enrichedElementCoords(end) - enrichedElementCoords(end-1));

%conductivity matrix of the locally enriched problem
Kenr = rbLocalConductivityMatrix(problem, enrichedElementCoords)*problem.B_map(problem.coords(end-1), problem.coords(end))*detJenr;
%capacity matrix of the locally enriched problem
Menr = rbLocalCapacityMatrix(problem, enrichedElementCoords)*problem.F_map(problem.coords(end-1), problem.coords(end))*1/detJenr;
%load vector of the locally enriched problem
fenr = rbLocalLoadVector(problem, enrichedElementCoords, problem.N)*problem.F_map(problem.coords(end-1), problem.coords(end))*1/detJenr;
  
end