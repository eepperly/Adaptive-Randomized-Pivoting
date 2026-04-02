function [idx,L] = rejection_helper_m(A,w)
    n = size(A,1);
    idx = zeros(n,1);
    L = A;
    out_pos = 0;

    for i = 1:n
        if w(i) * rand() >= L(i,i)
            continue;
        end

        out_pos = out_pos + 1;
        idx(out_pos) = i;

        % Scale L(i:n,i) by 1/sqrt(L(i,i)).
        L(i:n,i) = L(i:n,i) / sqrt(L(i,i));

        % Rank-1 update: L(i+1:n,i+1:n) -= v*v', where v = L(i+1:n,i).
        if i < n
            v = L(i+1:n,i);
            L(i+1:n,i+1:n) = L(i+1:n,i+1:n) - v * v';
        end
    end
end
