function LM = locationMapXIGA(modes, numberOfXcontrolPoints, p)

% LOCATIONMAPXIGA it maps the shape functions to each
% element to the corresponding unknown index
%   modes = number of POD-modes
%   numberOfXelements = number of enriched knot spans
% p = BSpline order

LM = zeros( modes, numberOfXcontrolPoints);

%% Construct the location map of the XFEM-block
% Fill LM with locations of the modes

for j=1:numberOfXcontrolPoints
    if j == 1
        for i=1:modes
            LM(i, j) = (j-1) * modes + i;
        end
    else
        for i=1:modes
            LM(i, j) = (j-2) * modes + i;
        end
    end
end

end