function [idx,W] = arp(A,k,varargin)
    if ~isempty(varargin) && ~isempty(varargin{1})
        method = varargin{1};
    else
        method = "rejection";
    end

    if length(varargin)>1 && ~isempty(varargin{2})
        interp_method = varargin{2};
    else
        interp_method = "default";
    end

    Om = sparsestack(k,size(A,2),4)';
    [Q,~] = qr(full(A*Om),"econ");
    
    if strcmp(method,"rejection")
        [idx,U,R] = rejection_rpqr(Q');
    elseif strcmp(method, "rejection-cpp") || strcmp(method, "cpp")
        [idx,U,R] = rejection_rpqr(Q',[],[],true);
    elseif strcmp(method,"ck") || strcmp(method,"original")
        [idx,U,R] = rpqr_ck(Q');
    else
        error("Method %s not recognized",method)
    end

    if strcmp(interp_method,"default")
        W = (Q * U) / R';
    elseif strcmp(interp_method,"optimal")
        [Q,R] = qr(A(idx,:)',"econ");
        W = (A * Q) / R';
    elseif strcmp(interp_method,"sketchy")
        S = sparsestack(2*k,size(A,2),4)';
        AS = full(A*S);
        AselS = full(A(idx,:)*S);
        [Q,R] = qr(AselS',"econ");
        W = (AS * Q) / R';
    else
        error("Interpolation method %s not recognized",interp_method)
    end
end