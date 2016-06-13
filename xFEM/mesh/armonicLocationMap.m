function [ armonicBasisFun, LM ] = armonicLocationMap( N, numberOfModes )
%ARMONICLOCATIONMAP returns the location mp and the armonic basis function
%   N = number of elements
%   numberOfModes = number of armonic modal function


LM = zeros(N, 2 + numberOfModes*2);
for i=1:N
    LM(i, 1) = i;
    LM(i, 2) = i+1;
    for j=1:numberOfModes
        for k=1:2
            LM(i, 2*(j-1)+2+k) = (N+1) + j + (i-1)*2*numberOfModes + (j+ k-2);
        end
    end
end
    
 armonicBasisFun = @armonicBasis;

end

