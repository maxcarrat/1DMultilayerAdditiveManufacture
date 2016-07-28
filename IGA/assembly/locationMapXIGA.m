function LM = locationMapXIGA(modes, numberOfXelements, p)

% LOCATIONMAPXIGA it maps the shape functions to each
% element to the corresponding unknown index
%   modes = number of POD-modes
%   numberOfXelements = number of enriched knot spans
% p = BSpline order

LM = zeros( (p+1) * modes, numberOfXelements);

%% Construct the location map of the XFEM-block
% Fill LM with locations of the modes

for j=1:numberOfXelements
    if j == 1
        for i=1:(p+1)*modes
            LM(i, j) = (j-1) * modes + i;
        end
    else
        for i=1:(p+1)*modes
            LM(i, j) = (j-2) * modes + i;
        end
    end
end

end