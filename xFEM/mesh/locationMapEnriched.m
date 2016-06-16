function LM = locationMapEnriched(modes, numberOfXelements)

% LOCATIONMAPENRICHED it maps the shape functions to each
% element to the corresponding unknown index
%   modes = number of POD-modes
%   numberOfXelements = number of enriched elements

LM = zeros(numberOfXelements, 2*modes);

%% Construct the location map of the XFEM-block

% Fill LM with locations of the modes 
for i=1:numberOfXelements
    for j=1:2*modes
        LM(i, j) = (i-1)*2*modes + j;
    end
end

end