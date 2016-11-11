function [ LM ] = locationMapXtendedBaseOverlayCoupling(...
    numberOfOverlayElements, modes, numberOfEnrichedDofs, numberOfCoupledBaseDofs )
%LOCATIONMAPXTENDEDBASEOVERLAYCOUPLING  generate the location map of the coupling terms
%in multiscale X-POD IGA.
%Input:
%numberOfOverlayElements = number of elements in the overlay mesh
%modes = number of POD modes
%numberOfEnrichedDofs = number od DOFs in the overlay mesh enriched using
%POD modes
%numberOfCoupledBaseDofs = number of coupled degrees of freedom in the base
%mesh
%Output:
%LM = location map of copling terms

LM = zeros(numberOfOverlayElements, 2 + modes);

%% Construct the location map of the Coupling-block

% Fill LM with locations of the modes (Coupling X-POD and IGA meshes)
for i=1:numberOfOverlayElements
    
    LM(i,1) =  i;
    LM(i,2) =  i+1;
    
    if i <= 2
        for j=1:2*modes
            LM(i,2 + j) = (numberOfOverlayElements + 1) + (i-1)*modes + j;
        end
    else
        for j=1:2*modes
            LM(i,2 + j) = (numberOfOverlayElements + 1) + (i-2)*2*modes + j + 1;
        end
    end
end

end


