function R = quatTorot(quat)
    import casadi.*
    % Function to transform a quaternion to a rotation matrix.
    % The quaternion is normalized by dividing by (q' * q) -- be aware
    % that this is equivalent to 1 / (||q||^2), which is not the usual 
    % 1 / ||q||.  If you actually want a unit quaternion, you may need 
    % to replace (q'*q) with sqrt(q'*q).

    % Normalize
    q = quat / (quat' * quat);

    % Extract components
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    % Build the rotation matrix (3x3) using CasADi operations
    R = vertcat( ...
         horzcat( q0^2 + q1^2 - q2^2 - q3^2,  2*(q1*q2 - q0*q3),       2*(q1*q3 + q0*q2) ), ...
         horzcat( 2*(q1*q2 + q0*q3),        q0^2 + q2^2 - q1^2 - q3^2, 2*(q2*q3 - q0*q1) ), ...
         horzcat( 2*(q1*q3 - q0*q2),        2*(q2*q3 + q0*q1),        q0^2 + q3^2 - q1^2 - q2^2 ) ...
    );
end