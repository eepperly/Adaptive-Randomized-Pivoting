function [idx,W] = sketchy_pivoting(A,k)
    Om = sparsestack(2*k,size(A,2),4)';
    [~,~,idx] = qr(full(A*Om)',"vector");
    idx = idx(1:k);
    S = sparsestack(2*k,size(A,2),4)';
    AS = full(A*S);
    AselS = full(A(idx,:)*S);
    [Q,R] = qr(AselS',"econ");
    W = (AS * Q) / R';
end