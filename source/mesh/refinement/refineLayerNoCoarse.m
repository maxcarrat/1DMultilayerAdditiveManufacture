function [ refinedMesh ] = refineLayerNoCoarse( mesh, refinementDepth, layer, numberOfLayers, numberOfRefinedElementsToBeKept )
%REFINELAYER refine the elements layer of the initial mesh up to a given depth
%   mesh = original coarse mesh
%   refinementDepth = depth of the refinement
%   layer = actual layer
%   numberOfLayers = number of layers

%coarsening of the previous active configuration
mesh = linspace(mesh(1), mesh(end), size(mesh, 2));

%generate refined mesh for the new active configuration
[activeMesh, numberOfElementsPerLayer] = getLayerActiveCoords(mesh, layer, numberOfLayers);
refinedActiveElement = linspace(activeMesh(end-numberOfElementsPerLayer), activeMesh(end), numberOfElementsPerLayer*2^(refinementDepth)+ 1);
switch numberOfRefinedElementsToBeKept
    case(1)
        refinedMesh = [mesh(1:((layer-1) * numberOfElementsPerLayer)) refinedActiveElement];
        
    case(5)
        if layer == 1
            refinedMesh = [mesh(1:((layer-1) * numberOfElementsPerLayer)) refinedActiveElement];
        elseif layer == 2
            refinedActiveElement_2 = linspace(activeMesh(end-2*numberOfElementsPerLayer), activeMesh(end-numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedMesh = [mesh(1:((layer-2) * numberOfElementsPerLayer)) refinedActiveElement_2(1:end-1) refinedActiveElement];
            
        elseif layer == 3
            refinedActiveElement_2 = linspace(activeMesh(end-2*numberOfElementsPerLayer), activeMesh(end-numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_3 = linspace(activeMesh(end-3*numberOfElementsPerLayer), activeMesh(end-2*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedMesh = [mesh(1:((layer-3) * numberOfElementsPerLayer)) refinedActiveElement_3(1:end-1) refinedActiveElement_2(1:end-1) refinedActiveElement];
            
        elseif layer == 4
            refinedActiveElement_2 = linspace(activeMesh(end-2*numberOfElementsPerLayer), activeMesh(end-numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_3 = linspace(activeMesh(end-3*numberOfElementsPerLayer), activeMesh(end-2*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_4 = linspace(activeMesh(end-4*numberOfElementsPerLayer), activeMesh(end-3*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedMesh = [mesh(1:((layer-4) * numberOfElementsPerLayer)) refinedActiveElement_4(1:end-1) refinedActiveElement_3(1:end-1) refinedActiveElement_2(1:end-1) refinedActiveElement];
            
        else
            refinedActiveElement_2 = linspace(activeMesh(end-2*numberOfElementsPerLayer), activeMesh(end-numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_3 = linspace(activeMesh(end-3*numberOfElementsPerLayer), activeMesh(end-2*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_4 = linspace(activeMesh(end-4*numberOfElementsPerLayer), activeMesh(end-3*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedActiveElement_5 = linspace(activeMesh(end-5*numberOfElementsPerLayer), activeMesh(end-4*numberOfElementsPerLayer),...
                numberOfElementsPerLayer*2^(refinementDepth)+ 1);
            
            refinedMesh = [mesh(1:((layer-5) * numberOfElementsPerLayer)) refinedActiveElement_5(1:end-1) refinedActiveElement_4(1:end-1) refinedActiveElement_3(1:end-1) refinedActiveElement_2(1:end-1) refinedActiveElement];
            
        end
    otherwise
        disp('case not implemented!!!');
end
        
end

