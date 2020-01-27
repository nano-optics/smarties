function stRes1 = pstGetResStructOneLambda(stResAll, lambda0)
  %% pstGetResStructOneLambda
% Extracts a single lambda from the result structure to calculate fields
%
% Input:
%       stResAll: [struct] Result structure with all lambda's
%       lambda0:  [1 x 1] Value of lambda to extract. If not found,
%                  lambda(1) is used
%
% Output:
%       stRes1: A struct with all the same fields as stAbcdnm, as well as
%       all of the other input parameters as fields of that struct, all for
%       one lambda
%
% Dependency: 
% none

indLambda = find(stResAll.lambda == lambda0,1);
if isempty(indLambda)
    indLambda = 1;
end

stRes1=stResAll;
stRes1.lambda=stResAll.lambda(indLambda);
if ~isscalar(stResAll.epsilon2)
    stRes1.epsilon2=stResAll.epsilon2(indLambda);
end
if ~isscalar(stResAll.epsilon1)
    stRes1.epsilon1=stResAll.epsilon1(indLambda);
end
stRes1.cnm = stResAll.cnm(indLambda,:);
stRes1.dnm = stResAll.dnm(indLambda,:);
stRes1.pnm = stResAll.pnm(indLambda,:);
stRes1.qnm = stResAll.qnm(indLambda,:);

if size(stResAll.anm,1)>1
    stRes1.anm = stResAll.anm(indLambda,:);
end
if size(stResAll.bnm,1)>1
    stRes1.bnm = stResAll.bnm(indLambda,:);
end
