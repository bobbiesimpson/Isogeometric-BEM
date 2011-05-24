function [Gsubmatrix] = integrateGsubmatrix_Telles(ngp, NURBScurve, element, collocCoords, srcXi_param)
    % integrate the G submatrix using Telles' transformation
    
    [gpt gwt]=lgwt(ngp,-1,1);

    Gsubmatrix=zeros(2,6);
    
    % ------- Curve information ---------
    range = NURBScurve.elRange(element,:);
    bsFnConn = NURBScurve.bsFnConn(element,:);
    elcoords = NURBScurve.controlPts(bsFnConn,1:2);
    
    srcXi=convertToParentCoordSpace(srcXi_param, range);
    xiStar=srcXi^2 - 1;
    gamBar=nthroot( srcXi*xiStar + abs(xiStar) ,3) + nthroot(srcXi*xiStar - abs(xiStar),3) + srcXi;
    jacob_param=(range(2)-range(1)) / 2;    % jacobian from parent to parameter space
    
    for pt=1:ngp    % integrate using Gaussian quadrature
        gXi=gpt(pt);
        xi=((gXi-gamBar).^3 + gamBar*(gamBar.^2+3)) / (1+3*gamBar.^2);
        xi_param=convertToParamSpace( xi, range);    % get the gauss point in parameter space
        
        [Ndisp Ngeom dNgeom] = getDispAndGeomBasis(xi_param, NURBScurve, element);
        
        [jacob_xi, normals, r, dr, drdn] = getKernelParameters( elcoords, collocCoords, Ngeom, dNgeom );
        jacob=jacob_xi*jacob_param;     % the final jacobian we use
        
        Ttemp=zeros(2,2); Utemp=zeros(2,2);
        for i=1:2
            for j=1:2
                [Ttemp(i,j) Utemp(i,j)]=DBIEkernels(i,j,r,dr,drdn,normals);
            end
        end
        jacobTelles=3* ((gXi-gamBar).^2) / (1+3*gamBar.^2);
        Gsubmatrix=Gsubmatrix + [Ndisp(1)*Utemp Ndisp(2)*Utemp Ndisp(3)*Utemp]*jacob*gwt(pt)*jacobTelles;

    end
end

