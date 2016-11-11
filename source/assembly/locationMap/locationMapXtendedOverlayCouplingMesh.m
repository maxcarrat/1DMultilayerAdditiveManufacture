function [ LM ] = locationMapXtendedOverlayCouplingMesh(  numberOfOverlayElements, numberOfCoupledDofs )
%LOCATIONMAPXTENDEDOVERLAYMESH generate the location map of the coupling terms
%in multiscale X-POD IGA.
%Input:
%numberOfOverlayElements = number of elements in the overlay mesh
%numberOfCoupledBaseDofs = number of coupled degrees of freedom in the
%one-element overlay mesh
%Output:
%LM = location map of copling terms

LM = zeros(numberOfCoupledDofs, 2);

%% Construct the location map of the Coupling-block

% Fill LM with locations of the modes (Coupling POD modes and FEM)
for i=1:numberOfCoupledDofs
      for j=1:2
          LM(i, j) = (numberOfOverlayElements-numberOfCoupledDofs) + i + (j-1);
      end
end

end

