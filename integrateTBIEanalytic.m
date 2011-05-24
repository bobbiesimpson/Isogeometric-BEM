function [ Hsubmatrix ] = integrateTBIEanalytic( NURBScurve, element, srcXi_param, collocNormal, collocCoords )

% -----------------------------------------------------------------------
% --------- Analytical integration of the hypersingular TBIE integral ---
% -----------------------------------------------------------------------

% --------------------
% ------ Inputs ------
% --------------------

% NURBScurve    :       The NURBS curve we are integrating on
% element       :       The element (local to the curve) we are integrating over
% scrXi_param   :       The collocation point in paramters space
% collocNormal  :       The normal at the collocation point

global E mu

intConst = E / (4 * pi * (1 - mu^2));

range = NURBScurve.elRange(element,:);                  % range of element in parameter space
bsFnConn = NURBScurve.bsFnConn(element,:);
elcoords = NURBScurve.controlPts(bsFnConn,1:2);

srcXi=convertToParentCoordSpace(srcXi_param, range);    % coordinate of src pt in parent space
    
% get the normal of the element (assume it is flat and normal is a constant
% throughout

[dispShape geomShape geomShapeDeriv] = getDispAndGeomBasis(srcXi_param, NURBScurve, element);

[jacob_xi, normals] = getKernelParameters( elcoords, collocCoords, geomShape, geomShapeDeriv );
el_length = jacob_xi * (range(2) - range(1));

I1 = 0.75 * ( ( 3.0 * srcXi - 1.0 ) * log( abs( ( 1.0 - srcXi ) / ( 1 + srcXi ) ) ) +...
( 6.0 * srcXi * srcXi - 2.0 * srcXi - 3.0 ) / ( srcXi * srcXi - 1.0 ) );
I2 = 0.5 * ( 9.0 * srcXi * log ( abs( ( 1.0 + srcXi ) / ( 1.0 - srcXi ) ) ) -...
( 18.0 * srcXi * srcXi - 13.0 ) / ( srcXi * srcXi - 1.0 ) );
I3 = 0.75 * ( ( 3.0 * srcXi + 1.0 ) * log( abs( ( 1.0 - srcXi ) / ( 1.0 + srcXi ) ) ) ...
+ ( 6.0 * srcXi * srcXi + 2.0 * srcXi - 3.0 ) / ( srcXi * srcXi - 1.0 ) );

Hsubmatrix = intConst * 2 / el_length * dot(collocNormal,normals) .* [I1 0 I2 0 I3 0;
                                                                      0 I1 0 I2 0 I3];


end

