function [ stressTensor ] = getExactStressTensor( traction, coords )
% Get the stress tensor in cartesian coordinates for the exact solution of
% a hole in an infinite plate subject to a traction at infinity

r=norm(coords);
theta=atan2(coords(2),coords(1));

%% Hole problem stresses
% We assume that the origin of the problem is that the centre of the hole
% 
% a=1;
% sigma_rr = traction/2 * ( 1 - (a^2)/(r^2) ) + traction / 2 * ( 1 - 4 * (a^2)/(r^2) ...
%           + 3 * (a^4)/(r^4) ) * cos(2*theta);
% sigma_rt = -traction/2 * ( 1 + 2 * (a^2)/(r^2) - 3 * (a^4)/(r^4) ) * sin(2*theta);
% sigma_tt = traction/2 * ( 1 + (a^2)/(r^2) ) - traction/2 * ( 1 + 3 * (a^4)/(r^4) ) ...
%             * cos(2*theta);
% 
% polarMat=[sigma_rr sigma_rt; sigma_rt sigma_tt];
% transMat=[cos(theta) sin(theta); -sin(theta) cos(theta)];
% 
% stressTensor=transMat'*polarMat*transMat;

%% The L-plate parameters

alpha=3*pi/4;
gamma = 2.565819161212361/(2 * alpha);
gamma2=4.281342942588585/(2*alpha);
an=1; bn=0;     % an=1,bn=0 gives pure mode 1 loading and vice versa for mode 2

[sigmaxx sigmayy sigmaxy]=getWedgeStresses( r, theta, gamma, gamma2, alpha, an, bn );

stressTensor = [sigmaxx sigmaxy; sigmaxy sigmayy];
end

