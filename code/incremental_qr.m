classdef incremental_qr < handle
% Implementation of the incremental QR decomposition of a tall matrix
% with an increasing number of columns

    properties
        H   % Matrix to store evolving QR decomposition
        tau % Vector to store scalars for Householder reflections
        k   % Number of columns in current matrix
        n   % Height of matrix
    end
    
    methods
        function obj = incremental_qr(Y)
            % Initialize the QR decomposition with matrix Y
            [obj.n,obj.k] = size(Y);
            obj.H = zeros(obj.n);
            obj.tau = zeros(obj.n,1);
            [obj.H(:,1:obj.k), obj.tau(1:obj.k)] =...
                matlab.internal.decomposition.compactQR(Y);
        end
        
        function obj = addcols(obj, Ynew)
            % Add new columns to the existing QR decomposition
            l = size(Ynew,2);
            
            Ynew = obj.applyQtfull(Ynew);   % Apply Qt to Ynew
            Ybot = Ynew(obj.k+1:end,:);     % Extract bottom of Ynew

            % Form QR decomposition of bottom
            [Ynew(obj.k+1:end,:),obj.tau(obj.k+1:obj.k+l)] =...
                matlab.internal.decomposition.compactQR(Ybot);

            obj.H(:,obj.k+1:obj.k+l) = Ynew; % Write to buffer
            obj.k = obj.k + l; % Increase size of k
        end
        
        function Qx = applyQ(obj, x)
            % Apply the Q matrix to a vector x
            Qx = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, [x;zeros(obj.n-obj.k,...
                size(x,2))], false, obj.k);
        end
        
        function Qtx = applyQt(obj, x)
            % Apply the transpose of the Q matrix to a vector x
            Qtx = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, x, true, obj.k);
            Qtx = Qtx(1:obj.k,:);
        end
        
        function Qtx = applyQtfull(obj, x)
            % Apply the transpose of the Q matrix to a vector x
            Qtx = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, x, true, obj.k);
        end
        
        function y = projectOut(obj, x)
            % Compute y = (I-Q*Q')*x
            y = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, x, true, obj.k);
            y(1:obj.k,:) = 0;
            y = matlab.internal.decomposition.applyHouseholder(...
                obj.H, obj.tau, y, false, obj.k);
        end
        
        function Q = getQ(obj)
            % Return the Q matrix from the current QR decomposition
            Q = obj.applyQ(eye(obj.k));
        end
        
        function R = getR(obj)
            % Return the R matrix from the current QR decomposition
            R = triu(obj.H);  % Extract the upper triangular part
            R = R(1:obj.k,1:obj.k);
        end
    end
end