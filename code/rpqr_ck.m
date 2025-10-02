function [idx,Q,R] = rpqr_ck(B)
    [k,n] = size(B);
    idx = zeros(k,1);
    Q = eye(k);
    for i = 1:k
        idx(i) = datasample(1:n,1,"Weights",sqcolnorms(B(i:end,:)));
        u = B(i:end,idx(i));
        e1 = zeros(k-i+1,1); e1(1) = 1;
        u = u + sign(sum(u))*norm(u)*e1; u = u / norm(u);
        B(i:end,:) = B(i:end,:) - 2*u*(u'*B(i:end,:));
        Q(:,i:end) = Q(:,i:end) - 2*(Q(:,i:end)*u)*u';
    end
    R = B(:,idx);
end