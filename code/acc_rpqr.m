function [idx,W] = acc_rpqr(A,k)
    A = full(A');
    [m,n] = size(A);
    idx = zeros(k,1);
    num_picked = 0;
    Q = zeros(m,k);
    F = zeros(n,k);
    while num_picked < k
        prop = datasample(1:n,k,"Replace",true,"Weights",sqcolnorms(A));
        H = A(:,prop); H = H'*H;
        [sel,~] = rejection_sample(H,diag(H));
        sel = prop(sel);
        sel = sel(1:min(length(sel),k-num_picked));
        idx(num_picked+1:num_picked+length(sel)) = sel;
        [Qnew,~] = qr(A(:,sel),"econ");
        Q(:,num_picked+1:num_picked+length(sel)) = Qnew;
        F(:,num_picked+1:num_picked+length(sel)) = A'*Qnew;
        A = A - Qnew*F(:,num_picked+1:num_picked+length(sel))';
        num_picked = num_picked + length(sel);
    end
    W = F/F(idx,:);
end