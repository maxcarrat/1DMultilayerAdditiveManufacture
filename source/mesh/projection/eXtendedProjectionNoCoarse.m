function [ projectedCoefficients ] = eXtendedProjectionNoCoarse(problem, previousTemperature,...
    modes, updatedMesh, previousMesh, PODRefinementDepth, initialTemperature)
% EXTENDEDPROJECTIONNOCOARSE project the previous solution onto the updated mesh at the
% new time step considering the enrichement modes from POD.
%   previousTemperature = temeprature distribution of the previous mesh
%   updatedMesh = actual mesh configuration
%   coarseMesh = initial/coarse mesh configuration
%   modes = number of enrichment modes
%   initialTemperature = initial temperature of the powder

projectedCoefficients = zeros(problem.N + 1 +...
    modes * (2*problem.XN - 1), 1);

projectedCoefficients = projectOntoEnrichedMesh(problem, previousTemperature, modes, updatedMesh,...
    previousMesh, PODRefinementDepth, initialTemperature);

end