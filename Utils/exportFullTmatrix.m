function [ F ] = exportFullTmatrix( stT, out, format)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 2)
    out = [];
end
if(nargin < 3) 
    format = '%d %d %.15e %.15e\n';
end


[T] = sparseTmatrix(stT);
F = full(T);

[r,c] = size(F);

im = repmat(1:r,c,1);
jm = im';

i=im(:); j = jm(:);

values = [i j real(F(:)) imag(F(:))];

% write to a file
if(~isempty(out))
    fileID = fopen(out, 'w');
    fprintf(fileID, '%% %d elements of T-matrix\n', r);
    fprintf(fileID, '%% i j Tr Ti \n');
    fprintf(fileID, format, values.');
    fclose(fileID);
end



end

