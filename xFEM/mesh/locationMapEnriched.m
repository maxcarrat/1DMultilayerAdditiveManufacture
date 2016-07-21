function LM = locationMapEnriched(modes, numberOfXelements)

% LOCATIONMAPENRICHED it maps the shape functions to each
% element to the corresponding unknown index
%   modes = number of POD-modes
%   numberOfXelements = number of enriched elements

LM = zeros(2 * modes, numberOfXelements);

%% Construct the location map of the XFEM-block

% Fill LM with locations of the modes

for j=1:numberOfXelements
    if j == 1
        for i=1:2*modes
            
            LM(i, j) = (j-1) * modes + i;
            
        end
    else
        for i=1:2*modes
            
            LM(i, j) = (j-2) * modes + i;
            
        end
    end
end

end