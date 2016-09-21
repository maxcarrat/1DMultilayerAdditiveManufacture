function [ refinedMesh, numberOfElementsPerLayer ] = refineLayerEnriched( activeMesh, numberOfElementsPerLayer, refinementDepth,...
    layer, numberOfTrainingLayers, PODRefinementDepth, numberOfRefinedElementsToBeKept )
%REFINELAYERENRICHED refine the elements layer of the initial mesh up to a given depth
%   activeMesh = original coarse mesh
%   refinementDepth = depth of the refinement


% generate refined mesh for the new active configuration

switch numberOfRefinedElementsToBeKept
    case(0)
        refinedMesh = activeMesh;
    case(1)
        % refine with a given depth (PODRefinementDepth) the POD enriched
        % element.
        numberOfElementsPerLayer = numberOfElementsPerLayer*2^(PODRefinementDepth);
        
        refinedActiveElement = linspace(activeMesh(end-1), activeMesh(end),...
            numberOfElementsPerLayer + 1);
        refinedMesh = [activeMesh(1:end-2) refinedActiveElement];
        
    case(5)
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

