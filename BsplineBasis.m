function [ N dN ] = BsplineBasis( i,p,xi, knot )
% calculate the non-rational Bspline basis function

% i=basis fn we want
% p=order (0=const, 1=linear, 2=quadratic)
% xi=parametric coordinate along b-spline
% knot= the knot vector (eg. [0 0 0 1 2 3 4 4 5 5 5]) in parameter space

basis=zeros(p+1,p+1);   % this will contain the basis functions of everything of lower order

% [ N_1,0 N_2,0 N_3,0 ..
%   N_2,0 N_2,1 ...
%   ETC.

% first calculate the basis fns for order 0

% The first if statement is the conventional interval for the zeroth order
% B-spline functions. But if we want the curve to be interpolatory at the
% end of the parameter space, we need to put another check in. Basically,
% we want the zeroth order basis function for the max non-zero interval to
% be equal to one at the end of the interval.
for c=0:p
    if(xi>=knot(i+c) &&  xi<knot(i+c+1))   % do we lie within the interval xi_i<=xi<xi_(i+!)?
        basis(1,c+1)=1;
    elseif xi==max(knot) && knot(i+c+1)==max(knot) && knot(i+c)~=knot(i+c+1)
        basis(1,c+1)=1;
    end
    if (p==0)           % if the basis order is zero we have the fn directly
        N=basis(1,c+1);
        dN=0;
        return
    end
end

% and now for the higher order fns
for b=2:p+1
    start=p+1-b+1;
    finish=1;
    for c=start:-1:finish
        basis(b,c)=BsplineBasisHighOrder(basis(b-1,c:(c+1)),xi, knot, i+c-1, b-1);
    end
end

N=basis(p+1,1);

temp=p / ( knot(i+p) - knot(i) ) * basis(p,1);
if isnan(temp) 
    temp=0;
end
temp2=p / ( knot(i+p+1) - knot(i+1) ) * basis(p,2);
if isnan(temp2)
    temp2=0;
end
dN = temp-temp2;

end

