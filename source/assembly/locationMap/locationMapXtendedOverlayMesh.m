function [ LM ] = locationMapXtendedOverlayMesh( modes, numberOfXElements )
%LOCATIONMAPXTENDEDOVERLAYMESH  it maps the shape functions to each
% element to the corresponding unknown index
%   modes = number of POD-modes
%   numberOfXDofs = number of enriched knot spans

LM = zeros( numberOfXElements, 2*modes );

%% Construct the location map of the XFEM-block
% Fill LM with locations of the modes

for i=1:numberOfXElements
    if i <= 2
        for j=1:2*modes
            LM(i, j) = (i-1) * modes + j;
        end
    else
        for j=1:2*modes
            LM(i, j) = modes + (i-2) * 2 * modes + j;
        end
    end
end

end

