function [ Hsubmatrix ] = integrateDBIEanalytic( NURBScurve, element, srcXi_param, jumpTerm )

% ------ We wish to integrate an element analytically using the DBIE -----
% ------ This should only get called if we are on a crack surface --------
% ------ since the analytical integrals ARE NOT VALID for NURBS ----------

global const3 const4

% --- if we are not on a crack --------
if strcmp(NURBScurve.curveType,'crack') ~= 1
    return
    
% ------- we are on a crack ----------
else
    
    integConst = const3 * const4;           % integration constant
    
    range = NURBScurve.elRange(element,:);  % range of element in parameter space
    
    srcXi=convertToParentCoordSpace(srcXi_param, range);    % coordinate of src pt in parent space
    srcNdisp = getDispAndGeomBasis(srcXi_param, NURBScurve, element);
    
    h = [0 -1; 1 0] .* integConst;
    I1 = 0.75 * ( 0.5 * srcXi * ( 3.0 * srcXi - 2.0 ) * log( abs( ( 1.0 - srcXi ) / ( 1.0 + srcXi ) ) )...
        + 3.0 * srcXi - 2.0 );
    I2 = 0.5 * ( 0.5 * ( 3.0 * srcXi - 2.0 ) * ( 3.0 * srcXi + 2.0 ) ...
        * log( abs( ( 1.0 + srcXi ) / ( 1.0 - srcXi ) ) ) - 9.0 * srcXi );
    I3 = 0.75 * ( 0.5 * srcXi * ( 3.0 * srcXi + 2.0 ) * log( abs( ( 1.0 - srcXi ) / ( 1.0 + srcXi ) ) )...
        + 3.0 * srcXi + 2.0 );
    Hsubmatrix = [h.*I1 h.*I2 h.*I3] + ...
                 [srcNdisp(1).*jumpTerm srcNdisp(2).*jumpTerm srcNdisp(3).*jumpTerm];
    
end

end

