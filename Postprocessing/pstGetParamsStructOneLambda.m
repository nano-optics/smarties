function stParams1 = pstGetParamsStructOneLambda(stParamsAll, lambda0)
  %% pstGetParamsStructOneLambda
% Extracts a single lambda from the result structure to calculate fields
%
% Input:
%       stParamsAll: [struct] Params structure with all lambda's
%       lambda0:  [1 x 1] Value of lambda to extract. If not found,
%                  lambda(1) is used
%
% Output:
%       stParams1: A struct with all the same fields, all for
%       one lambda
%
% Dependency: 
% none

indLambda = find(stParamsAll.lambda == lambda0,1);
if isempty(indLambda)
    indLambda = 1;
end

stParams1=stParamsAll;
stParams1.lambda=stParamsAll.lambda(indLambda);
if ~isscalar(stParamsAll.epsilon2)
    stParams1.epsilon2=stParamsAll.epsilon2(indLambda);
end
if ~isscalar(stParamsAll.epsilon1)
    stParams1.epsilon1=stParamsAll.epsilon1(indLambda);
end
if ~isscalar(stParamsAll.s)
    stParams1.s=stParamsAll.s(indLambda);
end
if ~isscalar(stParamsAll.k1)
    stParams1.k1=stParamsAll.k1(indLambda);
end
