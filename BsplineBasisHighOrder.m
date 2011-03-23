function [ basis ] = BsplineBasisHighOrder( lowerBasis, xi, knot, i, p )
% the recursive function that gives us the basis functions for any p>0

if(lowerBasis(1) < eps) basis=0;
else
    basis=(xi-knot(i))/(knot(i+p)-knot(i)) * lowerBasis(1);
end

if(lowerBasis(2) > eps) basis=basis +...
        (knot(i+p+1) - xi)/(knot(i+1+p)-knot(i+1)) * lowerBasis(2);
end

end

