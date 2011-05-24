function [ Hsubmatrix ] = integrateHsubmatrixSST( ngp, NURBScurve, element, collocCoords, srcXi_param, jumpTerm)

% A function which calculate the H matrix for the singular case using the
% subtraction of singularity technique by Guiggiani and Casalini

% ---------------------------------
% ----------- INPUTS --------------
% ---------------------------------

% ngp:              # of gauss points
% NURBScurve:       the NURBS curve we are on
% element:          the element on the current NURBS curve we are integrating over
% collocCoords:     the coordinates of the collocation points (do we need this?)
% srcXi_param:      the coordinate of our source point in parameter space
% jumpTerm          the previously calculated jump term.

% ---------------------------------

global const3 const4

[gpt gwt]=lgwt(ngp,-1,1);               % get the gauss points

Hsubmatrix=zeros(2,6);

range = NURBScurve.elRange(element,:);
bsFnConn = NURBScurve.bsFnConn(element,:);
elcoords = NURBScurve.controlPts(bsFnConn,1:2);
    
if srcXi_param==range(1)            % Annoying, but I need to evaluate the parameters 
    nudgedXi=srcXi_param+eps;       % at the source point a small distance away from the actual point
elseif srcXi_param==range(2)
    nudgedXi=srcXi_param-eps;
else 
    nudgedXi=srcXi_param;
end

srcXi=convertToParentCoordSpace(srcXi_param, range);
[srcNdisp srcNgeom srcdNgeom] = getDispAndGeomBasis(nudgedXi, NURBScurve, element);

hterm=[0 -const4*const3; const4*const3 0];
htermMatrix=[hterm.*srcNdisp(1) hterm.*srcNdisp(2) hterm.*srcNdisp(3)];

jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space

for pt=1:ngp    % integrate using Gaussian quadrature
    xi=gpt(pt);
    xi_param=convertToParamSpace( xi, range);    % get the gauss point in parameter space
            
    [Ndisp Ngeom dNgeom] = getDispAndGeomBasis(xi_param, NURBScurve, element);
    
    [jacob_xi,normals, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, Ngeom, dNgeom );
    jacob=jacob_xi*jacob_param;     % the final jacobian we use
    
    Ttemp=zeros(2,2);
    
    for i=1:2
        for j=1:2
            Ttemp(i,j)=DBIEkernels(i,j,r,dr,drdn,normals);
        end
    end
    
    tempMatrix=[Ndisp(1)*Ttemp Ndisp(2)*Ttemp Ndisp(3)*Ttemp]...
        *(xi-srcXi)*jacob;
    tempMatrix=tempMatrix-htermMatrix;
    tempMatrix=tempMatrix./(xi-srcXi);
    Hsubmatrix=Hsubmatrix + tempMatrix*gwt(pt);
        
end

% and add the analytical integral on at the end
if abs(srcXi)<(1-100*eps)
    Hsubmatrix=Hsubmatrix+htermMatrix*log(abs((1-srcXi)/(1+srcXi)));
else
    jacob_xi=getKernelParameters(elcoords, collocCoords, srcNgeom, srcdNgeom);
    jacob_s=jacob_xi*jacob_param;
    beta_m=1/jacob_s;
    Hsubmatrix=Hsubmatrix+htermMatrix*log(abs(2/beta_m))*sign(xi-srcXi);
end

% and add the jump terms at the end
if abs(srcXi_param-range(2)) > 100*eps
    jumpMatrix=[jumpTerm.*srcNdisp(1) jumpTerm.*srcNdisp(2) jumpTerm.*srcNdisp(3)];
    Hsubmatrix=Hsubmatrix+jumpMatrix;
end

end

