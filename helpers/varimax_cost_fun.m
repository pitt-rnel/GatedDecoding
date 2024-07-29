function vcost = varimax_cost_fun(A)
    % standardize
    h = sqrt(sum(A.^2, 2)); 
    A = bsxfun(@rdivide, A, h); 
    p = size(A,1); 
    % varimax cost function
    vcost = 1/p*sum(sum(A.^4,1)- 1/p*sum(A.^2).^2); 
end