Functions and brief description
1) qDemo % a simple demo of quaternions in action -- rotating a point around an arbitrary axis -- animation
2) Q = qConj( Q1 )  % quaternion conjugation
3) R = qGetR( Qrotation ) % get rotation matrix equivalent to rotation quaternion
4) Q = qGetRotQuaternion( teta, rotationVector ) % creates quaternion that rotates by teta angle around vector rotationVector
5) Q = qInv( Q1 )   % inverse of Q1 quaternion, i.e. Q = Q1^(-1)
6) d = qLength( Q ) % norm (length) of Q
7) Q = qMul( Q1, Q2, Q3, ..., Q10 )    % quaternion multiplication Q = Q1*Q2*Q3*....
8) Q = qNormalize( Q1 ) % normalize Q1, i.e. Q = Q1/||Q1||, ||Q||=1
9) Protated = qRotatePoint( P, Qrotation ) % rotates point P according to quaternion Qrotation
10) Q = qGetQ( R ) - get a quaternion equivalent to a 3x3 rotation matrix