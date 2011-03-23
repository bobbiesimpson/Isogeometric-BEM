function [ R dR ] = NURBSbasis( i, p, xi, knot, weight )
% We pass in 
%   i = basis function number
%   p = order or b-spline functions
%   xi = the coordinate we are evaluating at (in parameter space)
%   knot = the vector of knot values
%   weights = a vector of all the 'weights' required for NURBS

[N_ip dN_ip]=BsplineBasis(i,p,xi,knot);

numBasisFns=length(knot)-1-p;
W_xi=0; dW_xi=0;

for ihat=1:numBasisFns
    [N dN]=BsplineBasis(ihat,p,xi,knot);
    W_xi=W_xi + N*weight(ihat);
    dW_xi=dW_xi + dN*weight(ihat);
end


R=N_ip*weight(i)/W_xi;
dR = weight(i) * ( W_xi*dN_ip - dW_xi*N_ip )/ (W_xi^2);

end

