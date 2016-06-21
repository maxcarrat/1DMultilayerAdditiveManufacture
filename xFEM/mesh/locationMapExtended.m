function LM = locationMapExtended(numberOfElements, numberOfXelements, modes)

% LOCATIONMAPEXTENDED it maps the shape functions to each
% element to the corresponding unknown index
%   N = number of coarse mesh elements
%   modes = number of POD-modes

LM = zeros(numberOfElements + 2*modes, 2 + 2*modes);

%% Construct the location map of the FEM

% LM of FEM dofs...
for i=1:numberOfElements
    LM(i, 1) = i;
    LM(i, 2) = i+1;
end

% ...fill LM with locations of the modes (Coupling FEM-XFEM)
for i=(numberOfElements-numberOfXelements+1):numberOfElements
    for j=1:2*modes
        LM(i, 2 + j) = numberOfElements + 1 + j +...
            (i-(numberOfElements-numberOfXelements+1))*2*modes;
    end
end


%% Extend LM with the XFEM dofs

%XFEM LM...
iCounter = numberOfElements + 1;
for i=numberOfElements+1:2*modes:numberOfElements+2*modes*numberOfXelements
    for k=1:2*modes
        for j=1:2*modes
            
            LM(iCounter, j) = numberOfElements + 1 + j + (i - (numberOfElements+1));

        end
        iCounter = iCounter+1;
    end
end

%... LM of XFEM-FEM Coupling
% for i=(numberOfElements+1):(numberOfElements+2*modes*numberOfXelements)
%     iCounter = iCounter+1;
%     
% end

iCounter = numberOfElements + 1;
for i=1:numberOfXelements
    for k=1:2*modes
        LM(iCounter, 2*modes + 1) = numberOfElements - numberOfXelements - 1 ...
            + 1 + (i-1);
        LM(iCounter, 2*modes + 2) = numberOfElements - numberOfXelements - 1 ...
            + 2 + (i-1);
        
        iCounter = iCounter+1;
    end
end


end