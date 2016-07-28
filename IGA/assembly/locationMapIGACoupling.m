function LM = locationMapIGACoupling(numberOfElements, numberOfXelements, p)

% LOCATIONMAPIGACOUPLING it maps the shape functions to each
% element to the corresponding unknown index
%   numberOfElements = number of knot spans
%   numberOfXelements = number of enriched knot spans
%   p = BSpline functions order

LM = zeros(numberOfXelements, p+1);

%% Construct the location map of the Coupling-block

% Fill LM with locations of the modes (Coupling IGA-XIGA)
% for i=1:numberOfXelements
%       for j=1:p+1
%           LM(i, j) = (numberOfElements-numberOfXelements) + i + (j-1);
%       end
% end

for i=1:numberOfXelements
      for j=1:p+1
          LM(i, j) = (numberOfElements-numberOfXelements) + i + (j-1);
      end
end

end