function [ refinedMesh ] = refineLayerEnriched( activeMesh, numberOfElementsPerLayer, refinementDepth,...
    layer, numberOfTrainingLayers, numberOfRefinedElementsToBeKept )
%REFINELAYERENRICHED refine the elements layer of the initial mesh up to a given depth
%   activeMesh = original coarse mesh
%   refinementDepth = depth of the refinement


%generate refined mesh for the new active configuration

switch numberOfRefinedElementsToBeKept
    case(0)
        refinedMesh = activeMesh;
    case(4)
        refinedActiveElement_1 = linspace(activeMesh(end-4*numberOfElementsPerLayer), activeMesh(end-3*numberOfElementsPerLayer),...
            numberOfElementsPerLayer*2^(refinementDepth)+ 1);
        refinedActiveElement_2 = linspace(activeMesh(end-3*numberOfElementsPerLayer), activeMesh(end-2*numberOfElementsPerLayer),...
            numberOfElementsPerLayer*2^(refinementDepth)+ 1);
        refinedActiveElement_3 = linspace(activeMesh(end-2*numberOfElementsPerLayer), activeMesh(end-numberOfElementsPerLayer),...
            numberOfElementsPerLayer*2^(refinementDepth)+ 1);
        refinedActiveElement_4 = linspace(activeMesh(end-5*numberOfElementsPerLayer), activeMesh(end-4*numberOfElementsPerLayer),...
            numberOfElementsPerLayer*2^(refinementDepth)+ 1);
        
        refinedMesh = [activeMesh(1:end-6) refinedActiveElement_4(1:end-1) refinedActiveElement_1(1:end-1) refinedActiveElement_2(1:end-1)...
            refinedActiveElement_3(1:end) activeMesh(end)];
    otherwise
        disp('Case not implemented yet!!!');
end
end

