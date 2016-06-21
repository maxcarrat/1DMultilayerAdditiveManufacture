function [ refinedMesh ] = refineMesh( mesh, refinementDepth, time, timeSteps )
%REFINEMESH refine the last element of the initial mesh up to a given depth
%   mesh = original coarse mesh
%   refinementDepth = depth of the refinement
%   time = actual time step
%   timeSteps = number of time steps

%coarsening of the previous active configuration
mesh = linspace(mesh(1), mesh(end), size(mesh, 2));

%generate refined mesh for the new active configuration
activeMesh = getActiveCoordinates(mesh, time, timeSteps);
activeMesh = getLayerActiveCoords
refinedActiveElement = linspace(activeMesh(end-1), activeMesh(end), 2^(refinementDepth) + 1);
refinedMesh = [activeMesh(1:end-2) refinedActiveElement];

end

