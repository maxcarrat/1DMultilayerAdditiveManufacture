function LM = locationMapCoupling(numberOfElements, numberOfXelements)

% LOCATIONMAPCOUPLING it maps the shape functions to each
% element to the corresponding unknown index
%   N = number of coarse mesh elements
%   modes = number of POD-modes

LM = zeros(numberOfXelements, 2);

%% Construct the location map of the Coupling-block

% Fill LM with locations of the modes (Coupling FEM-XFEM)
for i=1:numberOfXelements
        LM(i, 1) = (numberOfElements-numberOfXelements) + i;
        LM(i, 2) = (numberOfElements-numberOfXelements) + i + 1;

end

end