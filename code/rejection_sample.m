function [idx,L] = rejection_sample(A,w)
n = size(A,1);

if n <= 1000
    [idx,L] = rejection_helper(A,w);
    idx = idx(idx~=0);
    L = tril(L(idx,idx));
    return
end

k = floor(n/2);
[idx1,L11] = rejection_sample(A(1:k,1:k),w(1:k));

B = L11 \ A(idx1,k+1:end);
S = A(k+1:end,k+1:end) - B'*B;
[idx2,L22] = rejection_sample(S,w(k+1:end));
idx = [idx1;k+idx2];
L = [L11 B(:,idx2);zeros(length(idx2),length(idx1)) L22];
end