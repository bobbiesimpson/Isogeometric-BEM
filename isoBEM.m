% This is a BEM isogeometric analysis code

% RNS 2011

clc
clear all
%close all

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

% ---------------------------------------
% ------- Boundary Conditions -----------
% ---------------------------------------

infinitePlate=0;     % if this flag is set, we apply the infinite plate tractions.                                          
tractionAtInfinity=100;

if infinitePlate        % We have an infintie plate problem
    tractionX=0;
    tractionY=0;
else                    % We have a finite problem
    tractionX=0;   
    tractionY=200;
end

% create a struct which contains all boundary condition information
BoundCondition = struct('presXTracRange', [0 3; 4 6]./7, ...
                        'presYTracRange', [1 7]./7,...
                        'nonZeroXTracRange',[],...  
                        'nonZeroYTracRange',[2 3]./7,... 
                        'zeroXDispRange', [3 4; 6 7]./7,...  
                        'zeroYDispRange', [0 1/7]);
                    
% ------------------------------------
% ------- crack parameters -----------
% ------------------------------------

crackRange = [4 5; 5 6]./7;      % define ranges of crack [ upper ; lower]

% ------------------------------------
% ------- Mesh generation  -----------
% ------------------------------------

p=2;        % degree of basis functions

numMeshes=10; inc=5;

L2relNorm=zeros(ceil(numMeshes/inc),2);

meshCounter=0;
    
for mesh=10:inc:numMeshes
    meshCounter=meshCounter+1;
    refinement=mesh;
    [ NURBScurve ]=generateBEMmesh( p, refinement );
    
    % create a global matrix of all the collocation points 
    collocPts = [NURBScurve(1).collocPts NURBScurve(2).collocPts NURBScurve(3).collocPts];
    dispConn = [NURBScurve(1).dispConn; NURBScurve(2).dispConn; NURBScurve(3).dispConn];
    tracConn = [NURBScurve(1).tracConn; NURBScurve(2).tracConn; NURBScurve(3).tracConn];
    tracDispConn = [NURBScurve(1).tracDispConn; NURBScurve(2).tracDispConn(2:end); NURBScurve(3).tracDispConn(2:end)];
    elRange = [NURBScurve(1).elRange; NURBScurve(2).elRange; NURBScurve(3).elRange];
    collocCoords = [NURBScurve(1).collocCoords; NURBScurve(2).collocCoords; NURBScurve(3).collocCoords];
    
    nDof=length(collocPts)*2;    % number of Dof (2 x number of control points)
    ne=size(dispConn,1);         % number of 'elements'
    nPts=length(collocPts);      % number of control points (the first and last are shared)
    numCurves = size(NURBScurve,2);
    
    % ---------------------------------------
    % ------- Boundary conditions -----------
    % ---------------------------------------
    
    dispDofs=(1:max(dispConn(end,:))*2)';     % all the DOFs (each node has x and y components)
    tracDofs=(1:max(tracConn(end,:))*2)';

    [ presDispDOFs, presTracDOFs, dirichletVals, nonZeroXTracDOFs, nonZeroYTracDOFs ]=...
        assignDirichletAndNeumannNodes(refinement, dispConn, tracConn, elRange, BoundCondition);
    
    unknownDispDofs=setxor(dispDofs,presDispDOFs);  % and all those that aren't prescribed are unknown
    unknownTracDofs=setxor(tracDofs,presTracDOFs);
    
    % -----------------------------------
    % ------- IsoBEM analysis -----------
    % -----------------------------------
    
    ngp_s=12;       % # gauss points for singular integrals
    ngp_r=12;       % # gauss points for regular integrals
    
    % ------- Initialise ---------
    H=zeros(nDof,nDof);                      % H matrix
    A=zeros(nDof,nDof);                      % A matrix
    G=zeros(nDof,max(tracConn(:,p+1))*2);    % G matrix
    z=zeros(nDof,1);                         % z matrix
    
    % ------- Find normals @ colloc pts ---------
    collocNormals=findNormalsAtCollocationPoints(nPts, collocPts, NURBScurve);
    
    % ------- loop over collocation points ----
    for c=1:nPts
        srcXi_param=collocPts(c);               % colloc pt (parameter space)
        jumpTerm=calculateJumpTerm(collocNormals(c,:,1), collocNormals(c,:,2));
        
        % ------- loop over NURBS curves ----
        for curve=1:numCurves
            
            ne = NURBScurve(curve).ne;                  % number of elements (on curve)
            knotVec = NURBScurve(curve).knotVec;        % knot vector for curve
            controlPts = NURBScurve(curve).controlPts;  % control points for curve
            elRange = NURBScurve(curve).elRange;        % element ranges
            bsFnConn = NURBScurve(curve).bsFnConn;      % basis function connectivity for curve
        
            % ------- loop over elements on curve ----
            for element=1:ne
                range=elRange(element,:);
                glbBsFnConn=bsFnConn(element,:);                 % connectiviy of NURBS basis fns
                dElConn=NURBScurve(curve).dispConn(element,:);   % element connectivity for displacement nodes
                tElConn=NURBScurve(curve).tracConn(element,:);   % element connectivity for traction nodes
                
                elcoords=controlPts(glbBsFnConn,1:2);            % coordinates of element nodes
                collocGlbPt=collocCoords(c,:);                   % coordinates of collocation point
                
                mirrorRange = getMirrorElement(range, crackRange);  % get "mirror" elmt range, if applicable
                
                % -------- check for collocation on upper crack surface --------
                if srcXi_param > crackRange(1,1) && srcXi_param < crackRange(1,2)
                    
                    Gsubmatrix = zeros(2,6);
                    
                    % --- we are in the same element ------
                    if (srcXi_param <= range(2)) && (srcXi_param >= range(1))
                       % fprintf('Upper crack, DBIE singular, analytical, same element\n') 
                        
                        Hsubmatrix = integrateDBIEanalytic(NURBScurve(curve), element, srcXi_param, jumpTerm);
                        
                    % --- we are in the "mirror" element ---
                    elseif (srcXi_param <= mirrorRange(2)) && (srcXi_param >= mirrorRange(1))
                        %fprintf('Upper crack, DBIE singular, analytical, opposite element\n')
                        
                        mirror_srcXi_param = (range(1) - srcXi_param) + mirrorRange(2);
                        Hsubmatrix = integrateDBIEanalytic(NURBScurve(curve), element, mirror_srcXi_param, jumpTerm);
                        
                    % --- we have a non-singular integral
                    else
                       % fprintf('Upper crack, DBIE nonsingular, G_L\n')
                        [Hsubmatrix Gsubmatrix] = integrateHGsubmatrices_GLQ(ngp_r, NURBScurve(curve), element, collocGlbPt);
                    end
                    
                % -------- check for collocation on lower crack surface --------
                elseif srcXi_param > crackRange(2,1) && srcXi_param < crackRange(2,2)
                    
                    Gsubmatrix = zeros(2,6);
                    
                    % --- we are in the same element ------
                    if (srcXi_param <= range(2)) && (srcXi_param >= range(1))
                        %fprintf('Lower crack, TBIE analytical singular, same element\n')

                        Hsubmatrix = integrateTBIEanalytic( NURBScurve(curve), element, srcXi_param, collocNormals(c,:,1), collocGlbPt);
                        
                        % --- we are in the "mirror" element ---
                    elseif (srcXi_param <= mirrorRange(2)) && (srcXi_param >= mirrorRange(1))
                        %fprintf('Lower crack, TBIE analytical singular, opposite element\n')

                        mirror_srcXi_param = (range(2) - srcXi_param) + mirrorRange(1);
                        Hsubmatrix = integrateTBIEanalytic( NURBScurve(curve), element, mirror_srcXi_param, collocNormals(c,:,1), collocGlbPt);
                        
                        % --- we have a non-singular integral
                    else
                        %fprintf('Lower crack, TBIE nonsingular G_L. c=%d element=%d\n', c, element)
                        
                        [Hsubmatrix, Gsubmatrix] = integrateHGsubmatricesTBIE_GLQ( ngp_r, NURBScurve(curve), element, collocGlbPt, collocNormals(c,:,1));
                        
                        
                    end
                    
                % -------- we're collocating on a non-crack surface
                % --------   
                else
                    if (srcXi_param <= range(2)) && (srcXi_param >= range(1) || (element==ne && srcXi_param==0 && curve==numCurves))
                        %fprintf('Non-crack, DBIE singular, SST, Telles\n')
                        if element==ne && srcXi_param==0 && curve == numCurves
                           srcXi_param=1;
                        end
                        
                        % we have a singular integral
                        Hsubmatrix=integrateHsubmatrixSST(ngp_s, NURBScurve(curve), element, collocGlbPt, srcXi_param, jumpTerm);
                        Gsubmatrix=integrateGsubmatrix_Telles(ngp_s, NURBScurve(curve), element, collocGlbPt, srcXi_param);
                    else
                        %fprintf('Non-crack, DBIE non-singular, G_L\n')
                        
                        % we have a regular integral
                        [Hsubmatrix, Gsubmatrix] = integrateHGsubmatrices_GLQ(ngp_r, NURBScurve(curve), element, collocGlbPt);
                    end
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
    
   %plotPrescribedTractions( nonZeroXTracDOFs, nonZeroYTracDOFs, tracDispConnDOF,presXTracs, presYTracs)
    
    % put all the knowns on the right hand side
    z=z-H(:,presDispDOFs)*dirichletVals';
    z=z+G(:,nonZeroXTracDOFs) * presXTracs;
    z=z+G(:,nonZeroYTracDOFs) * presYTracs;
    
    soln=A\z;   % and solve
    
    displacement=zeros(nDof,1);
    displacement(presDispDOFs)=dirichletVals;
    displacement(unknownDispDofs)=soln(unknownDispDofs);

    soln(globalTracUnknownDOF)=soln(globalTracUnknownDOF)*SF;  % multiply the tractions by the scale factor
    
    % plot the deformed profile
    plotDeformedProfile( displacement, NURBScurve)

    %L2relNorm(meshCounter,1)=nDof;
    %L2relNorm(meshCounter,2)=calculateL2BoundaryNorm( displacement, dispConn, bsFnConn, tractionAtInfinity);
        
    fprintf('Condition number of A matrix %2.2f\n', cond(A))
end

% profile viewer
% profsave(profile('info'),'profile_results')
% 
%figure(2); hold on
%loglog(L2relNorm(1:meshCounter,1), L2relNorm(1:meshCounter,2), 'k+-')

