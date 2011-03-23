% This is a BEM isogeometric analysis code

% RNS 2011

clc
clear all
close all

addpath C_files/

%profile on

% ------------------------------------
% ------- Global constants -----------
% ------------------------------------

global E mu const3 const4 const2 const1 shearMod kappa
global p knotVec controlPts elRange

% assume plane strain for this example

E=1e5;
mu=0.3;
shearMod=E/(2*(1+mu));
const4=(1-2*mu);
const3=1/(4*pi*(1-mu));
const2=(3-4*mu);
const1=1/(8*pi*shearMod * (1-mu));
kappa=3-4*mu;

infinitePlate=0;     % if this flag is set, we apply the infinite plate tractions.
                        % otherwise we apply uniform traction (tractionX)
                        
% these parameters are for the plate with a hole problem                       
tractionAtInfinity=100;
exactTracInterval=[2/3 1];

% and this is for the L-plate problem
%exactTracInterval = [0 1/2; 5/6 1];     % the knot interval over which exact tractions are applied

if infinitePlate
    tractionX=0;
    tractionY=0;
else
    tractionX=0;   % presribed traction on upper and left surfaces
    tractionY=-10;
end

% ------------------------------------
% ------- Mesh generation  -----------
% ------------------------------------

p=2;        % degree of basis functions

numMeshes=2; inc=5;

L2relNorm=zeros(ceil(numMeshes/inc),2);

meshCounter=0;

for mesh=2:inc:numMeshes
    meshCounter=meshCounter+1;
    refinement=mesh;
    [ controlPts, knotVec, collocPts, collocCoords, bsFnConn, dispConn, tracConn, elRange, tracDispConn ]=generateBEMmesh( p, refinement );
    
    nDof=length(collocPts)*2;    % number of Dof (2 x number of control points)
    ne=size(dispConn,1);         % number of 'elements'
    nPts=length(collocPts);      % number of control points (the first and last are shared)
    
    % ---------------------------------------
    % ------- Boundary conditions -----------
    % ---------------------------------------
    
    dispDofs=(1:max(dispConn(end,:))*2)';     % all the DOFs (each node has x and y components)
    tracDofs=(1:max(tracConn(end,:))*2)';
    
    [ presDispDOFs, presTracDOFs, dirichletVals, nonZeroXTracDOFs, nonZeroYTracDOFs ]=assignDirichletAndNeumannNodes(refinement, dispConn, tracConn);
    
    unknownDispDofs=setxor(dispDofs,presDispDOFs);  % and all those that aren't prescribed are unknown
    unknownTracDofs=setxor(tracDofs,presTracDOFs);
    
    
    % -----------------------------------
    % ------- IsoBEM analysis -----------
    % -----------------------------------
    
    ngp_s=12;       % # gauss points for singular integrals
    ngp_r=6;        % # gauss points for regular integrals
    
    H=zeros(nDof,nDof);             % initialise our global matrices
    A=zeros(nDof,nDof);
    G=zeros(nDof,max(tracConn(:,p+1))*2);
    z=zeros(nDof,1);
    
    collocNormals=findNormalsAtCollocationPoints(nPts, ne, elRange, bsFnConn, collocPts, dispConn, controlPts);
    
    for c=1:nPts
        srcXi_param=collocPts(c);               % the local coordinate of collocation point (parameter space)
        
        for element=1:ne
            range=elRange(element,:);
            glbBsFnConn=bsFnConn(element,:);    % connectiviy of NURBS basis fns
            dElConn=dispConn(element,:);        % element connectivity for displacement nodes
            tElConn=tracConn(element,:);        % element connectivity for traction nodes
            
            elcoords=controlPts(dElConn,1:2);        % coordinates of element nodes
            collocGlbPt=collocCoords(c,:);                   % coordinates of collocation point
            
            if (srcXi_param <= range(2)) && (srcXi_param >= range(1) || (element==ne && srcXi_param==0))
                
                if element==ne && srcXi_param==0
                    srcXi_param=1;
                end
                
                % we have a singular integral
                jumpTerm=calculateJumpTerm(collocNormals(c,:,1), collocNormals(c,:,2));
                
                Hsubmatrix=integrateHsubmatrixSST( ngp_s, elcoords, glbBsFnConn, collocGlbPt, srcXi_param, range, jumpTerm);
                Gsubmatrix=integrateGsubmatrix_Telles(ngp_s, elcoords, glbBsFnConn, collocGlbPt, srcXi_param, range);
            else

                % we have a regular integral
                [Hsubmatrix, Gsubmatrix] = integrateHGsubmatrices_GLQ(ngp_r, elcoords, glbBsFnConn, collocGlbPt, range);
            end
            
            % apply the submatrices to the global matrices
            
            sctrVec(1:2:5)=dElConn*2-1;
            sctrVec(2:2:6)=dElConn*2;
            rowSctrVec=[c*2-1 c*2];         % scatter vector for rows
            H(rowSctrVec,sctrVec)=H(rowSctrVec,sctrVec)+Hsubmatrix;
            
            sctrVec(1:2:5)=tElConn*2-1;
            sctrVec(2:2:6)=tElConn*2;
            G(rowSctrVec,sctrVec)=G(rowSctrVec,sctrVec)+Gsubmatrix;
            
            % and in the case of applying the infinite plate model, let's
            % calculate the exact tractiont terms
            
            if infinitePlate
                
                if (range(1) >= exactTracInterval(1,1) && range(2) <= exactTracInterval(1,2)) || ...
                   (range(1) >= exactTracInterval(2,1) && range(2) <= exactTracInterval(2,2))
                    
                    exactzterm=integrateExactTractionTerm(12, elcoords, glbBsFnConn, collocGlbPt, srcXi_param, range, tractionAtInfinity);
                    z(rowSctrVec)=z(rowSctrVec)+exactzterm;
                end
            end
            
        end
    end

    SF=abs(trace(H)/trace(G(1:nDof,1:nDof)));   % scale factor
    
    A(:,unknownDispDofs)=H(:,unknownDispDofs);
    globalTracUnknownDOF=setxor(dispDofs,unknownDispDofs);
    
    % we need to map the traction DOF to the displacement DOF
    tracDispConnDOF=zeros(length(tracDispConn)*2,1);
    tracDispConnDOF(1:2:end)=tracDispConn*2-1;
    tracDispConnDOF(2:2:end)=tracDispConn*2;
    mappedTractionDofs=tracDispConnDOF(unknownTracDofs);
        
    for i=1:length(mappedTractionDofs)
        A(:,mappedTractionDofs(i))=A(:,mappedTractionDofs(i))-G(:,unknownTracDofs(i))*SF;
    end
    
    % the prescribed Tractions
    presXTracs = ones(length(nonZeroXTracDOFs),1) * tractionX;
    presYTracs = ones(length(nonZeroYTracDOFs),1) * tractionY;
    
   % plotPrescribedTractions( nonZeroXTracDOFs, nonZeroYTracDOFs, tracDispConnDOF,presXTracs, presYTracs)
    
    % put all the knowns on the right hand side
    z=z-H(:,presDispDOFs)*dirichletVals';
    z=z+G(:,nonZeroXTracDOFs) * presXTracs;
    z=z+G(:,nonZeroYTracDOFs) * presYTracs;
    
    soln=A\z;   % and solve
    
    displacements=zeros(nDof,1);
    displacement(presDispDOFs)=dirichletVals;
    displacement(unknownDispDofs)=soln(unknownDispDofs);

    soln(globalTracUnknownDOF)=soln(globalTracUnknownDOF)*SF;  % multiply the tractions by the scale factor
    
    % plot the deformed profile
    plotDeformedProfile( displacement, nPts, controlPts, tractionAtInfinity )

    L2relNorm(meshCounter,1)=nDof;
    L2relNorm(meshCounter,2)=calculateL2BoundaryNorm( displacement, dispConn, bsFnConn, tractionAtInfinity);
end
% profile viewer
% profsave(profile('info'),'profile_results')
% 
%figure(2); hold on
%loglog(L2relNorm(1:meshCounter,1), L2relNorm(1:meshCounter,2), 'k+-')

