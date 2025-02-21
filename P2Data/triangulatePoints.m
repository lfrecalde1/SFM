function X_4xN = triangulatePoints(x1_h, x2_h, P1, P2)
%TRIANGULATEPOINTS A simple linear triangulation of matched points in two cameras.
%
%   x1_h : 3xN homogeneous image coords in camera 1
%   x2_h : 3xN homogeneous image coords in camera 2
%   P1   : 3x4 projection matrix for camera 1
%   P2   : 3x4 projection matrix for camera 2
%
%   X_4xN: 4xN homogeneous 3D coordinates

    N = size(x1_h, 2);
    X_4xN = zeros(4, N);
    
    for i = 1:N
        A = [ ...
            (x1_h(1,i)*P1(3,:) - P1(1,:));
            (x1_h(2,i)*P1(3,:) - P1(2,:));
            (x2_h(1,i)*P2(3,:) - P2(1,:));
            (x2_h(2,i)*P2(3,:) - P2(2,:));
        ];
        % Solve A * X = 0 in least-squares sense => SVD or 'pinv'
        [~,~,V] = svd(A);
        X = V(:,end);
        % Store
        X_4xN(:,i) = X;
    end
end