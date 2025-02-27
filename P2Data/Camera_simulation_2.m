%% Code to simulate a camera with internal and external parameters
clc; clear; close all;

% Create camera class (same as your original)
% ---------------------------------------------------------------
% Larger radial and tangential distortion
k1 = 0.1;    % bigger negative => stronger barrel distortion
k2 =  -0.1; 
k3 =  0.0;      
p1 =  0.00;   % non-zero tangential distortion
p2 = -0.00;   % negative sign for demonstration

cam = CentralCamera( ...
    'focal', 0.015, ...             % 15 mm focal length
    'pixel', 10e-6, ...             % 10 micrometers per pixel
    'resolution', [1280 1024], ...  % image size: (nu=1280, nv=1024)
    'centre', [640 512], ...        % principal point
    'k', [k1 k2 k3], ...           % radial distortion
    'p', [p1 p2],...
    'distortion', [k1 k2 k3 p1 p2]);                % tangential distortion

K = cam.K;  % camera intrinsics

% Define chessboard in the XY-plane (Z=0)
% ---------------------------------------------------------------
% Define the number of steps along each axis
nx = 4; 
ny = 4;
nz = 10;

% Define the step size along each axis
dx = 0.1; 
dy = 0.1;
dz = 0.1;

% Generate 3D coordinates for x, y, z
[x_vals, y_vals, z_vals] = ndgrid( ...
    0 : dx : (nx-1)*dx, ...
    0 : dy : (ny-1)*dy, ...
    0 : dz : (nz-1)*dz );

% Flatten them and pack into [x; y; z] format
points3D = [x_vals(:), y_vals(:), z_vals(:)]';

% Convert to homogeneous coordinates [x; y; z; 1]
points3D_H = [points3D; ones(1, size(points3D, 2))];


% Storage arrays
samples = 2;
data_uv = ones(3, size(points3D, 2), samples);  % (2, #points, #samples)
data_uv_noise = ones(3, size(points3D, 2), samples);  % (2, #points, #samples)
data_xy = ones(3, size(points3D, 2), samples); % (2, #points, #samples)

% Save rotations
R_plane = zeros(3, 3, samples);
t_plane = zeros(3, samples);

R_camera = zeros(3, 3, samples);
t_camera = zeros(3, samples);

% Camera positions and angles
roll_values = [0, -10, -10, 10];
yaw_values = [0, 20, 20, -20];
pitch_values = [0, 0, 0, 0];
z_values = [0, 0.2, 0.2, 0.4];
y_values = [0, 0, 0.1, 0];
x_values = [0, 0.1, 0.0, 0.1];


%% Aux values for normalized
F_identity = [1, 0, 0;...
              0, 1, 0;...
              0, 0, 1];
Identity = [1, 0, 0, 0;...
            0, 1, 0, 0;...
            0, 0, 1, 0];

% Create a figure for live plotting
figure('Name','Image Plane Animation','NumberTitle','off');
axis([0 cam.nu 0 cam.nv]);   % x from [0..1280], y from [0..1024]
axis ij;                     % flip y-axis to match image coords
hold on; grid on;
xlabel('u (pixels)'); ylabel('v (pixels)');
title('Projected Points in the Image Plane');
colors = lines(samples);  % Alternatively, try jet(numFrames) or hsv(numFrames)
for k = 1:samples
    
    % ------- 1) Generate random rotation angles (in degrees) -------
    roll_deg  = 0;  % from -5 to +5
    pitch_deg = 0;  % from -5 to +5
    yaw_deg   = 0;  % from -5 to +5
    
    % Convert chosen angles to radians (or fix some to zero as you did)
    roll  = deg2rad(roll_deg);
    pitch = deg2rad(pitch_deg);
    yaw   = deg2rad(yaw_deg);
    
    % ------- 2) Build the rotation matrix -------
    Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
    R_plane(:, :, k) = Rx * Ry * Rz;
    
    % ------- 3) Random translation in Z between [1.5, 2.0] -------
    translationZ = 2.5;
    translationy = -0.25;
    translationx = 0;
    t_plane(:, k) = [translationx; translationy; translationZ];
    

    %% Camera translation and rotation
    roll_camera  = roll_values(k);  % from -5 to +5
    pitch_camera = pitch_values(k);  % from -5 to +5
    yaw_camera   = yaw_values(k);  % from -5 to +5
    
    
    roll_camera  = deg2rad(roll_camera);
    pitch_camera = deg2rad(pitch_camera);
    yaw_camera   = deg2rad(yaw_camera);
    
    
    Rx = [1 0 0; 0 cos(roll_camera) -sin(roll_camera); 0 sin(roll_camera) cos(roll_camera)];
    Ry = [cos(pitch_camera) 0 sin(pitch_camera); 0 1 0; -sin(pitch_camera) 0 cos(pitch_camera)];
    Rz = [cos(yaw_camera) -sin(yaw_camera) 0; sin(yaw_camera) cos(yaw_camera) 0; 0 0 1];
    R_camera(:, :, k) = Rx * Ry * Rz;
    % ------- 3) Random translation in Z between [1.5, 2.0] -------
    translation_camera_Z = z_values(k);
    translation_camera_y = y_values(k);
    translation_camera_x = x_values(k);
    t_camera(:, k) = [translation_camera_x; translation_camera_y; translation_camera_Z];
    
    % ------- 4) Construct full 4x4 transform H_plane -------
    H_plane = [R_plane(:, :, k), t_plane(:, k); 0 0 0 1];
    H_camera = [R_camera(:, :, k), t_camera(:, k); 0 0 0 1];

    
    % ------- 5) Transform the chessboard points into world coords -------
    points3D_I_h = H_plane * points3D_H;  
    points3D_I   = points3D_I_h(1:3, :);
    Points3d_C = inv(H_camera)*points3D_I_h;

    %% Location origin respect to camara
    x_from_camara(:, k) = pinv(R_camera(:, :, k))*(-t_camera(:, k));
    %x_from_camara(:, k) = (t_camera(:, k));
    % ------- 6) Project onto image plane using the camera intrinsics -------
    %  or directly with 'cam.C'
    values_normalized = F_identity*Identity*Points3d_C;
    values_normalized = values_normalized(1:2, :)./values_normalized(3, :);
    radius = vecnorm(values_normalized);
    D = 1 + k1*radius.^2 + k2*radius.^4;
    x_warp = values_normalized.*D;
    x_warp = [x_warp; ones(1, length(x_warp))];


    pixels_aux = cam.K * x_warp;         % 3 x N in homogeneous
    pixels_aux = pixels_aux(1:2,:) ./ pixels_aux(3,:);  % 2 x N
    
    % Store the 2D points if you need them later:
    data_uv(1:2, :,k) = pixels_aux;
    data_uv_noise(1:2, :,k) = pixels_aux;
    % Increase this to get more outliers
    num_outliers = 0;
    outlier_indices = randperm(size(data_uv,2), num_outliers);

    % Increase or decrease these points by a large random amount
    % Adjust the multiplier (e.g., 25, 50, etc.) to control severity of outliers
    multiplier = 5;
    data_uv_noise(1, outlier_indices, k) = data_uv(1, outlier_indices, k) + multiplier * randn(size(outlier_indices));
    data_uv_noise(2, outlier_indices, k) = data_uv(2, outlier_indices, k) + multiplier * randn(size(outlier_indices));
    data_uv(1, outlier_indices, k) = data_uv(1, outlier_indices, k) + multiplier * randn(size(outlier_indices));
    data_uv(2, outlier_indices, k) = data_uv(2, outlier_indices, k) + multiplier * randn(size(outlier_indices));
    data_xy(1:3, :, k) = points3D(1:3, :);
    
    % ----------- Plot the new frame in the same figure -----------
    plot(pixels_aux(1,:), pixels_aux(2,:),'.', 'Color',colors(k,:), 'MarkerSize',12);
    plot(data_uv_noise(1,:, k), data_uv_noise(2,:, k),'*', 'Color',colors(k,:), 'MarkerSize',15);
    % Optional text
    text(50, 50, sprintf('Iteration %d', k), 'Color',colors(k,:),'FontSize',12);
    
    drawnow;         % update the figure
    pause(0.1);      % short pause to slow down animation
end


%% Section to compue the homography
X = data_xy;
U = data_uv;

figure('Name','Board Frames in Inertial Frame', 'Position', [100, 100, 1200, 1200]);
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Chessboard Frames in the Camera Coordinate System');

% Optionally plot the inertial frame axes at the origin:
plot3([0 0.05],[0 0],[0 0],'r','LineWidth',2); % X-axis (red)
plot3([0 0],[0 0.05],[0 0],'g','LineWidth',2); % Y-axis (green)
plot3([0 0],[0 0],[0 0.05],'b','LineWidth',2); % Z-axis (blue)
text(0, 0, 0, 'Inertial', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment','left');
% Set axis limits: x from -0.3 to 0.3, y from -0.3 to 0.3, and z from 0 to 0.8
xlim([-1, 1]);
ylim([-1, 1]);
zlim([-1, 3]);
% Adjust the view
view(13,20);
% Decide how "long" you want each axis to appear:
scale = 0.1;  % Tweak to suit

% Number of frames:
numFrames = size(R_plane,3);

% Create a colormap with a unique color for each frame.
colors = lines(numFrames);  % Alternatively, try jet(numFrames) or hsv(numFrames)

for k = 1:numFrames
    % Translation of frame k (relative to inertial)
    tx = t_plane(1,k);
    ty = t_plane(2,k);
    tz = t_plane(3,k);

    % Rotation matrix for frame k
    Rk = R_plane(:,:,k);

    % Extract (scaled) basis vectors for plotting the local axes
    xAxis = Rk(:,1) * scale;
    yAxis = Rk(:,2) * scale;
    zAxis = Rk(:,3) * scale;

    % Translation of frame k (relative to inertial)
    tx_camera = t_camera(1,k);
    ty_camera = t_camera(2,k);
    tz_camera = t_camera(3,k);

    % Rotation matrix for frame k
    Rk_camera = R_camera(:,:,k);

    % Extract (scaled) basis vectors for plotting the local axes
    xAxis_camera = Rk_camera(:,1) * scale;
    yAxis_camera = Rk_camera(:,2) * scale;
    zAxis_camera = Rk_camera(:,3) * scale;

    % Build transformation matrix from camera frame to inertial frame
    T_matrix = [Rk, [tx; ty; tz]; 0 0 0 1];

    % Transform points from the camera frame to the inertial frame.
    % Here, X(1:2, :, k) are the 2D coordinates in the camera frame (assumed to lie on Z=0).
    points2D = X(1:3, :, k);
    numPoints = size(points2D, 2);
    % Augment with zeros for Z (since points are on the image plane) and ones for homogeneous coordinates
    points_homogeneous = [points2D; ones(1, numPoints)];
    points3D_I = T_matrix * points_homogeneous;  
    
    % Choose a unique color for this frame
    currentColor = colors(k, :);
    % Plot the transformed points in the inertial frame
   
    plot3(points3D_I(1,:), points3D_I(2,:), points3D_I(3,:), 'o', ...
          'MarkerSize', 2, 'MarkerFaceColor', currentColor);

    % Plot each axis for the camera frame
    % X-axis in red
    plot3([tx, tx + xAxis(1)], [ty, ty + xAxis(2)], [tz, tz + xAxis(3)], 'r', 'LineWidth', 2);
    % Y-axis in green
    plot3([tx, tx + yAxis(1)], [ty, ty + yAxis(2)], [tz, tz + yAxis(3)], 'g', 'LineWidth', 2);
    % Z-axis in blue
    plot3([tx, tx + zAxis(1)], [ty, ty + zAxis(2)], [tz, tz + zAxis(3)], 'b', 'LineWidth', 2);

    plot3([tx_camera, tx_camera + xAxis_camera(1)], [ty_camera, ty_camera + xAxis_camera(2)], [tz_camera, tz_camera + xAxis_camera(3)], 'r', 'LineWidth', 2);
    % Y-axis in green
    plot3([tx_camera, tx_camera + yAxis_camera(1)], [ty_camera, ty_camera + yAxis_camera(2)], [tz_camera, tz_camera + yAxis_camera(3)], 'g', 'LineWidth', 2);
    % Z-axis in blue
    plot3([tx_camera, tx_camera + zAxis_camera(1)], [ty_camera, ty_camera + zAxis_camera(2)], [tz_camera, tz_camera + zAxis_camera(3)], 'b', 'LineWidth', 2);

    % Optionally label the frame
    text(tx_camera, ty_camera, tz_camera, sprintf('Camera %d', k), 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment','left');
    drawnow;         % update the figure
    pause(0.1);      % short pause to slow down animation

end

%% Computing transformation and plot approximations
% Create a figure for live plotting
figure('Name','Image Plane Animation','NumberTitle','off');
axis([0 cam.nu 0 cam.nv]);   % x from [0..1280], y from [0..1024]
axis ij;                     % flip y-axis to match image coords
hold on; grid on;
xlabel('u (pixels)'); ylabel('v (pixels)');
title('Projected Points in the Image Plane');
colors = lines(samples);  % Alternatively, try jet(numFrames) or hsv(numFrames)
for k=1:size(data_uv,3)-1
    uv_1 = data_uv(:, :, k);
    uv_2 = data_uv(:, :, k+1);

    %% Computing  the initial value for the optimizer her we are computing the elemental matrix
    F_normalization = fundamental_analytical(uv_1, uv_2, cam.K);

    %% Since we have tha calibration matrices we can direclty compute the elemental matrix
    %uv_1_aux = pinv(cam.K)*uv_1;
    %uv_2_aux = pinv(cam.K)*uv_2;
    % 
    %% Computing solution based on casadi. This solution is not robust to outliers
    [F_optimization] = FundamentalCasadi(uv_1, uv_2, F_normalization);
    
    % % 
    % % %% Computing solutiion using ransac for outliers
    [best_model, num_iterations] = fitRansac(uv_1, uv_2, 8, F_normalization, 0.01);
    [U, S, V ] = svd(best_model);
    vector = [S(1,1), S(2, 2), 0];
    F_rank = U * diag(vector) *V';
    
    % 
    [R, t] = recoverPoseFromFundamental(F_optimization, cam.K, uv_1(1:2, :)', uv_2(1:2, :)');
    R_opti_only(:, :, k) = R;
    t_opti_only(:, k) = R'*t;

    [R_ransac, t_ransac] = recoverPoseFromFundamental(best_model, cam.K, uv_1(1:2, :)', uv_2(1:2, :)');
    R_opti_ransac(:, :, k) = R_ransac;
    t_opti_ransac(:, k) = R_ransac'*(-t_ransac);
    
    %% Computing points world
    P1 = cam.K * [eye(3, 3), zeros(3, 1)];
    P2 =  cam.K * [R_ransac, t_ransac];
    pts3D_4xN = triangulatePoints(uv_1, uv_2, P1, P2);
    pts3D_4xN = pts3D_4xN./pts3D_4xN(4, :);
    

    H1 = [eye(3, 3), zeros(3,1); zeros(1, 3), 1];
    H2 = [R_ransac, t_ransac; zeros(1, 3), 1];
    

    %% Compute a better points based on casadi non-linear optimization
    x_init  = init_optimization_variables(t_ransac, R_ransac, pts3D_4xN);
    [x_vector_opt, x_trans_opt, R_quaternion_opt, distortion] = cameraCalibrationCasADi(pts3D_4xN, uv_2, cam.K, x_init);
    H2_casadi = [R_quaternion_opt, x_trans_opt; zeros(1, 3), 1];
    pts3D_4xN_casadi = [x_vector_opt; ones(1, size(pts3D_4xN,2))];

    [pixels_aux_1] = projection_values(H1, pts3D_4xN, k1, k2, K);
    [pixels_aux_2] = projection_values(H2, pts3D_4xN, k1, k2, K);
    [pixels_aux_2_casadi] = projection_values(H2_casadi, pts3D_4xN_casadi, distortion(1), distortion(2), K);


    %% Computing the error of the approximation
    error = norm(uv_2(1:2, :) - pixels_aux_2, 2)
    error_casadi =  norm(uv_2(1:2, :) - pixels_aux_2_casadi, 2)

    %% Plot results image
    plot(uv_2(1, :), uv_2(2, :),'.', 'Color',colors(k+1,:), 'MarkerSize',15);
    plot(pixels_aux_2(1,:), pixels_aux_2(2,:),'x', 'Color','g', 'MarkerSize',15);
    plot(pixels_aux_2_casadi(1,:), pixels_aux_2_casadi(2,:),'*', 'Color','b', 'MarkerSize',15);

    % Optional text
    text(50, 50, sprintf('Iteration %d', k), 'Color',colors(k,:),'FontSize',12);
    
    drawnow;         
    pause(0.1);      
end

%% getting elements of the images
scale = x_from_camara(3, 2)/t(3);
scale_2 = x_from_camara(3, 2)/t_ransac(3);
scale_casadi = x_from_camara(3, 2)/ H2_casadi(3, 4);

t = t *scale
t_ransac = t_ransac *scale_2
t_casadi = x_trans_opt * scale_casadi
x_from_camara

R'
R_ransac'
R_quaternion_opt'
R_camera

%% Check results
figure('Name','Validation Points', 'Position', [100, 100, 1200, 1200]);
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Validation Points');

% Optionally plot the inertial frame axes at the origin:
plot3([0 0.05],[0 0],[0 0],'r','LineWidth',2); % X-axis (red)
plot3([0 0],[0 0.05],[0 0],'g','LineWidth',2); % Y-axis (green)
plot3([0 0],[0 0],[0 0.05],'b','LineWidth',2); % Z-axis (blue)
text(0, 0, 0, 'Inertial', 'FontSize', 14, 'Color', 'k', 'HorizontalAlignment','left');
% Set axis limits: x from -0.3 to 0.3, y from -0.3 to 0.3, and z from 0 to 0.8
% Adjust the view
view(13,20);
plot3(pts3D_4xN(1,:), pts3D_4xN(2,:), pts3D_4xN(3,:), 'o', ...
          'MarkerSize', 2, 'MarkerFaceColor', currentColor);
plot3(pts3D_4xN_casadi(1,:), pts3D_4xN_casadi(2,:), pts3D_4xN_casadi(3,:),'*', 'Color','b', 'MarkerSize',15);
