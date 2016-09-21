function  LM = locationMapIGA(N, p)
% location matrix
% it maps the shape functions local to each element to a global unknown
% index

LM = zeros(N, p+1);
for i=1:N
    for j=1:p+1
        LM(i, j) = i + (j-1);
    end
end

end
