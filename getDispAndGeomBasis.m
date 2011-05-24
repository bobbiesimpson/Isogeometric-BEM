function [ shape geomShape geomShapeDeriv] = getDispAndGeomBasis( xi_param, NURBScurve, element )

% ---- Given the coordinate in the parameter space, find the basis
% ---- functions that we use the interpolate displacement -----
% ------------------ and the geometry -------------------------

% ------------------------
% ------- Inputs ---------
% ------------------------

% xi_param := the coordinate at which we wish to evaluate the basis functions (in parameter space)
% NURBScurce := the struct which defines a curve, be it a NURBS or crack 
% element := the local element defined in our curve range = [1,num_els_on_curve]

% ------------------------
% ------- Outputs --------
% ------------------------

% shape = displacement basis function
% geomShape = NURBS basis functions (geometry)
% geomShapeDeriv = NURBS basis function derivatives 

global p

% first check coordinate is in element, otherwise error
elRange = NURBScurve.elRange(element,:);
if (xi_param < elRange(1)) || (xi_param > elRange(2))
    shape=zeros(1,3);
    return
end

% We always calculate the NURBS basis function since all the boundary
% geometry is described by NURBS

bsFnConn = NURBScurve.bsFnConn(element,:);
numBasisFns=length(bsFnConn);
geomShape=zeros(1,numBasisFns);
geomShapeDeriv = zeros(1,numBasisFns);

for lclBasis=1:numBasisFns
    i=bsFnConn(lclBasis);
    [geomShape(lclBasis) geomShapeDeriv(lclBasis)] = NURBSbasis(i, p, xi_param, NURBScurve.knotVec, NURBScurve.controlPts(:,3)' );
end
        
switch NURBScurve.curveType
    case 'crack'
        shape = zeros(1,3);
        
        xi = 2 * (xi_param - elRange(1)) / ( elRange(2) - elRange(1) ) - 1;
        
        g = 2/3;    % this creates discontinuous fns as used by Aliabadi et al.
        switch element
            % --- we have a semi-discontinuous element (at end) ----
            case 1
                shape(1) = xi * 1/(1 + g) .* (-g+xi);
                shape(2) = 1 + xi / (g + g^2) .* ( -(1-g^2) + xi * (-1 - g));
                shape(3) = xi * 1 / (g + g^2) .* (1 + xi);
                
            % --- we have a semi-discontinuous element (at start) --
            case NURBScurve.ne
                shape(1) = xi * 1 / (g + g^2) .* (-1 + xi);
                shape(2) = 1 + xi / (g + g^2) .* ( (1-g^2) + xi * (-1 - g));
                shape(3) = xi * 1/(1 + g) .* (g + xi);
            % --- we have a discontinuous element --
            otherwise
                shape(1) = 1.125 * xi * (xi- 2/3);
                shape(2) = (1 - 1.5 * xi) * (1 + 1.5 * xi);
                shape(3) = 1.125 * xi * (xi + 2/3);
                
                %shape(1) = xi .* ( -1/(2*g) + xi * 1 / (2 * g^2) );
                %shape(2) = 1 - xi.^2 * 1 / g^2;
                %shape(3) = xi .* (1 / (2 * g) + xi * 1 / (2 * g^2));
        end
        
    % --- we must be on a NURBS curve ----
    otherwise
        shape = geomShape;
end

end


