function A = random_sparse(m,n,s)
    rows = zeros(1,n*s);
    for i = 1:n
        rows((i-1)*s+1:i*s) = datasample(1:m,s,"Replace",false);
    end
    cols = kron(1:n,ones(1,s));
    vals = randn(1,n*s);
    A = sparse(rows,cols,vals,m,n);
end