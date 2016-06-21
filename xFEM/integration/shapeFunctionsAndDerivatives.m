function [ N, B ] = shapeFunctionsAndDerivatives( localCoordinate )
%SHAPEFUNCTIONSANDDERIVATIVES evaluate the shape functions matrix N and its
%derivative B in local coordinates
%   localCoordinate = local coordinates of the integration point

N = 0.5 * [(1-localCoordinate), (1+localCoordinate)];
B = 0.5 * [-1, 1];

end

