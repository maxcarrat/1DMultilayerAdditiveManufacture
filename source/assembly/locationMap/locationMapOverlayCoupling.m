function [ LM ] = locationMapOverlayCoupling( numberOfOverlayElements, numberOfCoupledBaseDofs, p )
%LOCATIONMAPOVERLAYCOUPLING generate the location map of the coupling terms
%in multiscale h-IGA.
%Input:
%numberOfOverlayElements = number of elements in the overlay mesh
%numberOfCoupledBaseDofs = number of coupled degrees of freedom in the base
%mesh
%p = BSpline functions order
%Output:
%LM = location map of copling terms

LM = zeros(numberOfCoupledBaseDofs, p+1);

%% Construct the location map of the Coupling-block

% Fill LM with locations of the modes (Coupling IGA-XIGA)
for i=1:numberOfCoupledBaseDofs
      for j=1:p+1
          LM(i, j) = (numberOfOverlayElements-numberOfCoupledBaseDofs) + i + (j-1);
      end
end

end

