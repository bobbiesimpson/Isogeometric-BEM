function [ jacob, normals, r, dr, drdn ] = getKernelParameters( elcoords, collocCoords, N, dN )
 % calculate the parameters we need to evaluate the kernels

dxydxi=dN*elcoords;             % the geometry derivatives
jacob=norm(dxydxi);             % jacobian
normals=1/jacob * [ dxydxi(2) -dxydxi(1) ];

fieldPt=N*elcoords;
relDist=fieldPt-collocCoords;
r=norm(relDist);
dr=1/r * relDist;
drdn=dr*normals';

end

