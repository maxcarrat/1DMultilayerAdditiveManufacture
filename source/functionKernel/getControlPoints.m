function [ CPs ] = getControlPoints( layer, layerThickness, p )
%GETCONTROLPOINTS generate the control points vector for a given layer

CPs = linspace(0.0, layer*layerThickness, p+layer);

end

