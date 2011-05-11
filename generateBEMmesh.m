function [ NURBScurve ] = generateBEMmesh( p, refinement )

% This function generate the isogeom BEM mesh and returns
% controlPts    :   matrix of the control points
% knotVec       :   knot vector
% collocPts     :   collocation points in parametric space
% elConn        :   connectivity matrix which gives non-zero basis fns in element
% tracConn      :   traction connectivity matrix - specifies 'traction DOF' for each element
% elRange       :   the lower and upper range of the element in the parametric space ie [xi_lower, xi_upper]

% Assume p=2    :   the degree of the NURBS basis
% refinement    :   1,2,3 etc. As this number increases, the h-refinement increases



% Plate with a hole geometry

% knotVec=[0 0 0 1/9 1/9 1/3 1/3 2/3 2/3 7/9 7/9 1 1 1];
% controlPts=[-4 0 1; -2.5 0 1; -1 0 1; -1 1 1/sqrt(2); 0 1 1; 0 2.5 1; 0 4 1; -2 4 1; -4 4 1; -4 2 1; -4 0 1];



% L-plate problem
% knotVec = [0 0 0 1/6 1/6 2/6 2/6 3/6 3/6 4/6 4/6 5/6 5/6 1 1 1];
% controlPts = [0 -sqrt(2) 1; 1/sqrt(2) -1/sqrt(2) 1; sqrt(2) 0 1; 1/sqrt(2) 1/sqrt(2) 1; ...
%               0 sqrt(2) 1; -1/(2*sqrt(2)) 3/(2*sqrt(2)) 1; -1/sqrt(2) 1/sqrt(2) 1; -1/(2*sqrt(2)) 1/(2*sqrt(2)) 1;...
%               0 0 1; -1/(2*sqrt(2)) -1/(2*sqrt(2)) 1; -1/sqrt(2) -1/sqrt(2) 1; -1/(2*sqrt(2)) -3/(2*sqrt(2)) 1; 0 -sqrt(2) 1];



% Spanner problem
% knotVec = [0 0 0 1 1 2 2 3 3 4 5 6 6 7 7 8 9 10 10 11 11 12 12 13 13 13];
% controlPts = [ 0 0; -0.5 -0.5; -1 -1; -1.5, -1; -2 -1;...
%                -2 -1.25; -2 -1.5; -1 -3; 1 -1; 5 -1; 10 -1; 10 0;...
%                10 1; 5 1; 1 1; -1 3; -2 1.5; -2 1.25; -2 1; ...
%                -1.5 1; -1 1; -0.5 0.5; 0 0];
%

% Uniaxial tension

% knotVec = [ 0 0 0 1 1 2 2 3 3 4 4 4];
% controlPts = [ 0 0; 5 0; 10 0; 10 5; 10 10; 5 10; 0 10; 0 5;
%                 0 0];

% Crack problem
delta = 2;
crackUpper = 5 + delta;
crackLower = 5 - delta;
upperCPt = (10 + crackUpper)/2; lowerCPt = (0 + crackLower)/2;
middleUCP = (5 + crackUpper)/2; middleLCP = (5 + crackLower)/2;

knotVec = [ 0 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 7];
controlPts = [ 0 0; 5 0; 10 0; 10 5; 10 10; 5 10; 0 10; 0 upperCPt; 0 crackUpper; 2.5 middleUCP; 5 5; 
               2.5 middleLCP; 0 crackLower; 0 lowerCPt; 0 0];

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));

controlPts = [controlPts weights'];

uniqueKnots=unique(knotVec);
n=length(knotVec)-1-p;
weightedPts=[controlPts(:,1).*controlPts(:,3) controlPts(:,2).*controlPts(:,3) controlPts(:,3)];

% -------------------------------------
% ------------ h-refinement -----------
% -------------------------------------

% the idea here is to take the original knot vector and add knots at point
% in between existing points

if refinement>1
    oldKnotVec=knotVec;
    newKnots=zeros((length(uniqueKnots)-1),refinement-1);
    
    for i=1:size(newKnots,1)
        distributed=linspace(uniqueKnots(i),uniqueKnots(i+1),refinement+1);
        newKnots(i,:)=distributed(2:end-1);
    end
    [m n]=size(newKnots);
    newKnots=sort(reshape(newKnots,1,m*n));

    
    % and now redefine the control points
    for knot=1:length(newKnots)                                % loop over all the new knots
        newControlPts=zeros(length(weightedPts)+1,3);        % increase the array of control points
        knot_bar=newKnots(knot);
        kval=find(knotVec>knot_bar, 1 )-1;
        for i=1:length(newControlPts)
            if i<= (kval - p)
                alpha=1;
            elseif i>=(kval-p+1) && i<= kval
                alpha= ( newKnots(knot) - oldKnotVec(i) ) / ( oldKnotVec(i+p) - oldKnotVec(i) );
            else
                alpha=0;
            end 
            newPoints=zeros(1,3);
            if i~=1
                newPoints=(1-alpha).*weightedPts(i-1,:);
            end
            if i~= length(newControlPts)
               newPoints=newPoints + alpha * weightedPts(i,:);
            end
            newControlPts(i,:)=newPoints;
        end
        weightedPts=newControlPts;
        knotVec=sort([knotVec knot_bar]);
        oldKnotVec=knotVec;
    end

controlPts=[weightedPts(:,1)./weightedPts(:,3) weightedPts(:,2)./weightedPts(:,3) weightedPts(:,3)];
end

% --------------------------------------------------------
% ------------- Split geometry for crack analysis --------
% --------------------------------------------------------

if refinement==1
    % we split the element just before and after the crack entrance by
    % inserting a knot three times
    
    crackEntIndex = find(knotVec==4/7,1);       % point where crack begins
    crackExIndex = find(knotVec==6/7,1,'last'); % point where crack ends
    
    newKnots = [mean(knotVec(crackEntIndex-1:crackEntIndex))*ones(1,3) ...
        mean(knotVec(crackExIndex:crackExIndex+1))*ones(1,3)];
else
    % we decrease the continuity at the element boundary just before and
    % after the crack entrance by inserting a knot twice
    
    xiIndex1 = find(knotVec < 4/7,1,'last');       % point where crack begins
    xiIndex2 = find(knotVec > 6/7,1); % point where crack ends
    
    newKnots = [knotVec(xiIndex1)*ones(1,2) ...
        knotVec(xiIndex2)*ones(1,2)];

end

[ knotVec controlPts ] = knotInsertion( knotVec, controlPts, newKnots, p );

uniqueKnots=unique(knotVec);

% --------------------------------------------------------
% ------------- Create structs of NURBScurves ------------
% --------------------------------------------------------

NURBScurve = struct('knotVec', [], 'controlPts', [], 'ne', {}, 'elRange', [],...
                    'elConn', [], 'tracConn', [], 'tracDispConn', []);
NURBScurve(3).knotVec = [];

index1 = 1:find(knotVec==newKnots(1),1,'last');   % 
ctPtIndex1 = 1:length(index1)-p-1;

index2 = find(knotVec==newKnots(1),1):find(knotVec==newKnots(end),1,'last');   % 
ctPtIndex2 = ctPtIndex1(end)+1:ctPtIndex1(end)+length(index2)-p-1;

index3 = find(knotVec==newKnots(end),1):length(knotVec);   % 
ctPtIndex3 = ctPtIndex2(end)+1:size(controlPts,1);

NURBScurve(1).knotVec = knotVec(index1);
NURBScurve(1).controlPts = controlPts(ctPtIndex1,:);

NURBScurve(2).knotVec = knotVec(index2);
NURBScurve(2).controlPts = controlPts(ctPtIndex2,:);

NURBScurve(3).knotVec = knotVec(index3);
NURBScurve(3).controlPts = controlPts(ctPtIndex3,:);

% --------------------------------------------------------
% ------------- Define element connectivities ------------
% --------------------------------------------------------

% loop over each NURBS curve


previousKnotVal = 0;
dispDofCounter = 0;
tracDofCounter = 0;

for curve=1:size(NURBScurve,2)
    
    uniqueKnots = sort(unique(NURBScurve(curve).knotVec));
    
    ne = length(uniqueKnots)-1;
    knotVec = NURBScurve(curve).knotVec;
    
    NURBScurve(curve).ne=ne;   % number of elements
    NURBScurve(curve).elRange=zeros(ne,2);        % initialise matrices
    NURBScurve(curve).elConn=zeros(ne,p+1);
    elKnotIndices=zeros(ne,2);
    NURBScurve(curve).tracConn=zeros(ne,p+1);
    
    element = 1;
    
    % first determine our element ranges
    for i=1:length(knotVec)
        currentKnotVal=knotVec(i);
        if knotVec(i)~=previousKnotVal
            NURBScurve(curve).elRange(element,:)=[previousKnotVal currentKnotVal];
            elKnotIndices(element,:)=[i-1 i];
            element=element+1;
        end
        previousKnotVal=currentKnotVal;
    end
    
    
    % determine our element ranges and the corresponding knot indices
    numRepeatedKnots=0;

    % generate our connectivity matrices for disps and tractions
    for e=1:ne
        indices=(elKnotIndices(e,1)-p+1):elKnotIndices(e,1);
        previousKnotVals=knotVec(indices);
        currentKnotVals=ones(1,p)*knotVec(elKnotIndices(e,1));
        if isequal(previousKnotVals,currentKnotVals) && e~=1;
            numRepeatedKnots=numRepeatedKnots+1;
        end
        NURBScurve(curve).elConn(e,:)=((elKnotIndices(e,1)-p):elKnotIndices(e,1)) + dispDofCounter;
        NURBScurve(curve).tracConn(e,:)=((elKnotIndices(e,1)-p):elKnotIndices(e,1)) + ...
                                        numRepeatedKnots + tracDofCounter;
    end
    NURBScurve(curve).bsFnConn=NURBScurve(curve).elConn;
    
    dispDofCounter = NURBScurve(curve).elConn(end,end)-1;
    tracDofCounter = NURBScurve(curve).tracConn(end,end);
    
    % create connectivity matrix for traction -> disp DOF
    totalTractionDOF=length(knotVec)-p-1+numRepeatedKnots;
    NURBScurve(curve).tracDispConn=zeros(totalTractionDOF,1);
    for e=1:ne
        NURBScurve(curve).tracDispConn(NURBScurve(curve).tracConn(e,:)) = ...
                                                        NURBScurve(curve).elConn(e,:);
    end
    
end

NURBScurve(curve).elConn(end,end)=1;  % the last point is equal to the first point

% ------------------------------------------
% ------- Generate the collocation points--- 
% ------------------------------------------

% we loop over curves and generate collocation points where we assume that
% the next curve will generate the end collocation point of the present
% curve

for curve=1:size(NURBScurve,2)      % loop over curves
        
    knotVec = NURBScurve(curve).knotVec;
    controlPts = NURBScurve(curve).controlPts;
    NURBScurve(curve).collocPts = zeros(1,length(knotVec) - p -2);
    NURBScurve(curve).collocCoords = zeros(length(knotVec) - p -2,2);
    
    for i=1:length(knotVec)-p-2
        
        xi = sum(knotVec((i+1):(i+p)))/p;
        NURBScurve(curve).collocPts(i) = xi;
        
        NURBScurve(curve).collocCoords(i,1) = NURBSinterpolation(xi, p, knotVec, controlPts(:,1)', controlPts(:,3)');
        NURBScurve(curve).collocCoords(i,2) = NURBSinterpolation(xi, p, knotVec, controlPts(:,2)', controlPts(:,3)');

    end
    
end

% -----------------------------------
% ------- Plot the mesh - -----------
% -----------------------------------

figure(1); grid on
figure(2); grid on

for curve=1:size(NURBScurve,2)
    
    knotVec = NURBScurve(curve).knotVec;
    controlPts = NURBScurve(curve).controlPts;
    collocPts = NURBScurve(curve).collocPts;
    collocCoords = NURBScurve(curve).collocCoords;
    xi=linspace(min(knotVec),max(knotVec),300);
    
    n=size(xi,2);
    numBasisFns=length(knotVec)-1-p;
    NURBSvalues=zeros(n,numBasisFns);
    NURBSderivVals=zeros(n,numBasisFns);

    for i=1:numBasisFns
        
        for point=1:n
            [N dN] = NURBSbasis(i,p,xi(point),knotVec,controlPts(:,3)');
            NURBSvalues(point,i)=N;
            NURBSderivVals(point,i)=dN;
        end
        figure(1)
        hold on
        plot(xi,NURBSvalues(:,i), 'k-')
        hold off
        
    end
    
    
    NURBSCoords=NURBSvalues*controlPts;
    uniqueKnots = unique(NURBScurve(curve).knotVec);
    knotCoords=zeros(length(uniqueKnots),2);
    
    for point=1:(length(knotCoords)-1)
        knotCoords(point,1)=NURBSinterpolation(uniqueKnots(point), p, knotVec, controlPts(:,1)', controlPts(:,3)');
        knotCoords(point,2)=NURBSinterpolation(uniqueKnots(point), p, knotVec, controlPts(:,2)', controlPts(:,3)');
    end
    figure(2); hold on

    plot(NURBSCoords(:,1), NURBSCoords(:,2), 'k-', controlPts(:,1), controlPts(:,2), 'ro')
    plot(collocCoords(:,1), collocCoords(:,2), 'kx', 'MarkerSize', 10)
    plot(knotCoords(:,1), knotCoords(:,2), 'ks', 'MarkerSize', 6)
    hold off
    
end

axis equal

%legend('NURBS basis functions')


% 
% save 'dat_files/NURBS_geometry.dat' NURBSCoords -ASCII
% save 'dat_files/controlPts.dat' controlPts -ASCII
% save 'dat_files/collocPts.dat' collocCoords -ASCII
% save 'dat_files/elementEdges.dat' knotCoords -ASCII


% plot the exact circle geometry
% x=-1:0.01:0;
% y=sqrt(1-x.^2);
% figure(2); hold on
% plot(x,y,'g--')
% hold off; axis equal

%legend('Original geometry','Control points', 'Collocation points', 'Element edges')


end

