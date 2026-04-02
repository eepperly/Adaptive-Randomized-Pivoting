function [idx,Q,R] = rejection_rpqr(B,varargin)
    [m,n] = size(B);
    if isempty(varargin) || isempty(varargin{1})
        k = m;
    else
        k = varargin{1};
    end
    if length(varargin) < 2 || isempty(varargin{2})
        l = k;
    else
        l = varargin{2};
    end
    if length(varargin) < 3 || isempty(varargin{3})
        use_cpp = false;
    else
        use_cpp = true;
    end
    scn = sqcolnorms(B);
    idx = zeros(k,1);
    prop = datasample(1:n,l,"Replace",true,"Weights",scn);
    H = B(:,prop); H = H'*H;
    sel = rejection_sample(H,scn(prop),use_cpp);
    sel = prop(sel);
    idx(1:length(sel)) = sel; num_picked = length(sel);
    iqr = incremental_qr(B(:,sel));
    while num_picked < k
        prop = datasample(1:n,l,"Replace",true,"Weights",scn);
        C = iqr.projectOut(B(:,prop));
        H = C'*C;
        sel = rejection_sample(H,scn(prop));
        sel = prop(sel);
        idx(num_picked+1:num_picked+length(sel)) = sel;
        num_picked = num_picked + length(sel);
        iqr.addcols(B(:,sel));
    end
    Q = iqr.getQ(); R = iqr.getR();
end