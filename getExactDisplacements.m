function [ disp ] = getExactDisplacements( r, theta, traction )

disp=zeros(length(theta),2);

%% L-plate displacements

alpha=3*pi/4;
gamma = 2.565819161212361/(2 * alpha);
gamma2=4.281342942588585/(2*alpha);
an=1; bn=0;     % an=1,bn=0 gives pure mode 1 loading and vice versa for mode 2

[disp(:,1) disp(:,2)] = getWedgeDisplacements( r, theta, gamma, gamma2, alpha, an, bn );

%% Hole plate displacements

% global shearMod kappa

% a=1;
% disp=zeros(length(theta),2);
% 
% ur = traction * r .* cos(2*theta) / ( 4 * shearMod ) .* ( 1 + ( kappa + 1 ) .* (a^2)./(r.^2) - (a^4)./(r.^4)  ) ...
%     + traction * r / ( 8 * shearMod ) .* ( ( kappa - 1) + 2 .* (a^2)./(r.^2) );
% 
% utheta = -traction * r .* sin(2*theta) / ( 4 * shearMod) .* ( 1 + ( kappa - 1) * (a^2)./(r.^2) + (a^4)./(r.^4) );
% 
% disp(:,1)=cos(theta).*ur -sin(theta).*utheta;
% disp(:,2)= sin(theta).*ur + cos(theta).*utheta;

end

