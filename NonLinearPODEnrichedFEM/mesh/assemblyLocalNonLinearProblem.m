function [ Menr, Kenr, fenr ] = assemblyLocalNonLinearProblem(problem, solution)
%ASSEMBLYLOCALNONLINEARPROBLEM  returns the linear system of the non-linear
%problem
%   problem = struct of the non-linear transient heat problem
%   solution = actula temperature values

lastActiveElement = problem.N;

enrichedElementCoords = linspace(problem.coords(end-1), problem.coords(end),...
    2^problem.refinementDepth + 1);

detJenr = 2/(enrichedElementCoords(end) - enrichedElementCoords(end-1));

%conductivity matrix of the locally enriched problem
Kenr = rbNonLinearConductivityMatrix(problem, solution, enrichedElementCoords, lastActiveElement)*...
    problem.B_map(problem.coords(end-1), problem.coords(end)) * detJenr;
%capacity matrix of the locally enriched problem
Menr = rbLocalCapacityMatrix(problem, enrichedElementCoords)*...
    problem.F_map(problem.coords(end-1), problem.coords(end)) * 1/detJenr;
%load vector of the locally enriched problem
fenr = rbLocalLoadVector(problem, enrichedElementCoords, lastActiveElement)*...
    problem.F_map(problem.coords(end-1), problem.coords(end)) * 1/detJenr;
  
end

