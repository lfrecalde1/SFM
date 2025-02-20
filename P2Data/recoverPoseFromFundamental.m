function [R, t, inliers] = recoverPoseFromFundamental(F, K, pts1, pts2)

    E = F;
    
    
    [U, S, V] = svd(E);
    % Force the singular values to be [1, 1, 0]
    S(1,1) = 1; 
    S(2,2) = 1; 
    S(3,3) = 0;
    E = U * S * V';
    
    % Re-decompose after enforcing the singular values 
    [U, ~, V] = svd(E);
    
    % 3) Define W
    W = [0 -1 0; 
         1  0 0; 
         0  0 1];
    
    % 4) Two candidate rotations
    R1 = U * W  * V';
    R2 = U * W' * V';
    
    % 5) Candidate translation vectors (the 3rd column of U, +/- sign)
    u3 = U(:,3);
    t_candidates = [ u3, -u3 ];
    R_candidates = cat(3, R1, R2);
    
    % Ensure rotations are proper (det(R) = +1). If det < 0, flip sign.
    for i = 1:size(R_candidates, 3)
        if det(R_candidates(:,:,i)) < 0
            R_candidates(:,:,i) = -R_candidates(:,:,i);
        end
    end
    
    % Convert input points to homogeneous coords for triangulation
    pts1_h = [pts1, ones(size(pts1,1),1)]';
    pts2_h = [pts2, ones(size(pts2,1),1)]';
    
    % Camera 1 matrix: [I | 0]
    P1 = K * [eye(3), zeros(3,1)];
    
    bestCount = -inf;
    R = eye(3); 
    t = [0;0;0];
    inliers = true(size(pts1,1),1); % trivial init
    
    % 6) Evaluate all 4 combinations via cheirality check
    for iRot = 1:2
        for iT = 1:2
            R_test = R_candidates(:,:,iRot);
            t_test = t_candidates(:, iT);
            
            % Camera 2 matrix: [R | t]
            P2 = K * [R_test, t_test];
            
            % Triangulate points
            pts3D_4xN = triangulatePoints(pts1_h, pts2_h, P1, P2);  
            % 3D points are in homogeneous form [X; Y; Z; W].
            % Convert to Cartesian:
            pts3D = pts3D_4xN(1:3,:) ./ pts3D_4xN(4,:);
            
            % Check depth in Cam1: Z1 > 0
            Z1 = pts3D(3, :);
            % Check depth in Cam2: Z2 > 0
            % Transform pts3D into camera 2 frame: R_test*X + t_test
            pts3D_cam2 = R_test * pts3D + t_test;
            Z2 = pts3D_cam2(3, :);
            
            valid = (Z1 > 0) & (Z2 > 0);
            numInFront = sum(valid);
            
            if numInFront > bestCount
                bestCount = numInFront;
                R = R_test;
                t = t_test;
                inliers = valid(:);  % store which points are in front
            end
        end
    end
end


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