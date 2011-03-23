function [ controlPts, knotVec, collocPts, collocCoords, bsFnConn, elConn, tracConn, elRange, tracDispConn ] = generateBEMmesh( p, refinement )

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
% knotVec = knotVec/max(knotVec);
% weights = ones(1,size(controlPts,1));
% weights(9)=2;
% weights(15)=2;
% controlPts = [controlPts weights'];
% 
% uniqueKnots=unique(knotVec);
% n=length(knotVec)-1-p;
% weightedPts=[controlPts(:,1).*controlPts(:,3) controlPts(:,2).*controlPts(:,3) controlPts(:,3)];

%% Outer boundary of dam
knotVec=[0 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 16];
knotVec = knotVec/max(knotVec);
controlPts = [-3.000000 -20.000000;   97.000249/2 -20.0000015;     100.000249 -20.000003;     (100.000249+99.999847)/2 -10.0000015; 
          99.999847 0.000000;   (99.999847+90.359371)/2 0.000000;    90.359371 0.000000;     88.359375 0.000000; 
          88.359375 2.000002;   88.359375 (2.000002+5.078125)/2;     88.359375 5.078125;     (88.359375+71.173198)/2 5.078125; 
          71.173198 5.078125;   71.173198-7.550835 5.078125;  58.331384, 10.465140; (58.331384+43.671875)/2, (10.465140+25.390625)/2
          43.671875, 25.390625;   (43.671875+31.484375)/2, (25.390625+33.515625)/2;  31.484375, 33.515625; (31.484375+23.359375)/2, (33.515625+32.5000000)/2; 
          23.359375, 32.5000000;   (23.359375+24.176103)/2, (32.5000000+21.683710)/2;  24.176103, 21.683710;
          24.412440, 18.553823;   23.678620, 15.502011;(23.678620+20.227456)/2,(15.502011+1.149322)/2;
          20.227456, 1.149322; 19.951101028311882,-8.537042651468596e-06; 18.769013,0;
          15.769013/2,0;-3,0;-3,-10;-3,-20];
      
weights = ones(1,size(controlPts,1));
weights(8) = 1/sqrt(2);
weights(14) =0.9221;
weights(24) =0.987908;
weights(27) =0.785423;

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
uniqueKnots=unique(knotVec);
end

% --------------------------------------------------------
% ------------- Define element connectivities ------------
% --------------------------------------------------------


ne=length(uniqueKnots)-1;   % number of elements
elRange=zeros(ne,2);        % initialise matrices
elConn=zeros(ne,p+1);
elKnotIndices=zeros(ne,2);
tracConn=zeros(ne,p+1);

% determine our element ranges and the corresponding knot indices
element=1;
previousKnotVal=0;

for i=1:length(knotVec) 
    currentKnotVal=knotVec(i);
    if knotVec(i)~=previousKnotVal
        elRange(element,:)=[previousKnotVal currentKnotVal];
        elKnotIndices(element,:)=[i-1 i];
        element=element+1;
    end
    previousKnotVal=currentKnotVal;
end

numRepeatedKnots=0;
for e=1:ne
    indices=(elKnotIndices(e,1)-p+1):elKnotIndices(e,1);
    previousKnotVals=knotVec(indices);
    currentKnotVals=ones(1,p)*knotVec(elKnotIndices(e,1));
    if isequal(previousKnotVals,currentKnotVals) && length(nonzeros(previousKnotVals))>1;
        numRepeatedKnots=numRepeatedKnots+1;
    end
    elConn(e,:)=(elKnotIndices(e,1)-p):elKnotIndices(e,1);
    tracConn(e,:)=elConn(e,:)+numRepeatedKnots;
end
bsFnConn=elConn;
elConn(end,end)=1;  % the last point is equal to the first point

% create connectivity matrix for traction -> disp DOF
totalTractionDOF=length(knotVec)-p-1+numRepeatedKnots;
tracDispConn=zeros(totalTractionDOF,1);
for e=1:ne
    tracDispConn(tracConn(e,:))=elConn(e,:);
end

% -----------------------------------
% ------- Plot the mesh - -----------
% -----------------------------------


xi=linspace(0,max(knotVec),300);

n=size(xi,2);
numBasisFns=length(knotVec)-1-p;
NURBSvalues=zeros(n,numBasisFns);
NURBSderivVals=zeros(n,numBasisFns);
collocPts=zeros(1,numBasisFns);
r=zeros(n,numBasisFns);

figure(1); grid on

for i=1:numBasisFns

    collocPts(i)=sum(knotVec((i+1):(i+p)))/p;
    
    for point=1:n
        [N dN] = NURBSbasis(i,p,xi(point),knotVec,controlPts(:,3)');
        NURBSvalues(point,i)=N;
        NURBSderivVals(point,i)=dN;
        r(point,i)=abs(collocPts(i)-xi(point));
    end
    
    hold on
    plot(xi,NURBSvalues(:,i), 'k-', collocPts(i), 0, 'kx', uniqueKnots,0,'ks')
    hold off

end


legend('NURBS basis functions', 'collocation points', 'element edges')


% ---------------------------------------------------------
% ------- Get the coords of collocation points ------------
% ---------------------------------------------------------
% 

collocCoords=zeros(numBasisFns,2);
for point=1:(numBasisFns)  
    collocCoords(point,1)=NURBSinterpolation(collocPts(point), p, knotVec, controlPts(:,1)', controlPts(:,3)');
    collocCoords(point,2)=NURBSinterpolation(collocPts(point), p, knotVec, controlPts(:,2)', controlPts(:,3)');
end

NURBSCoords=NURBSvalues*controlPts;

% trim off last entries
collocCoords=collocCoords(1:(end-1),:);
collocPts=collocPts(1:(end-1));

% -------------------------------------------
% ------- Get the element points ------------
% -------------------------------------------

knotCoords=zeros(length(uniqueKnots),2);
for point=1:(length(knotCoords)-1)   
    knotCoords(point,1)=NURBSinterpolation(uniqueKnots(point), p, knotVec, controlPts(:,1)', controlPts(:,3)');
    knotCoords(point,2)=NURBSinterpolation(uniqueKnots(point), p, knotVec, controlPts(:,2)', controlPts(:,3)');

end

% ---------------------------------------------------------
% ------- Plot the mesh and collocation points ------------
% ---------------------------------------------------------
% 


figure(2); grid on
hold on
plot(NURBSCoords(:,1), NURBSCoords(:,2), 'k-', controlPts(:,1), controlPts(:,2), 'ro') 
plot(collocCoords(:,1), collocCoords(:,2), 'kx', 'MarkerSize', 10)
plot(knotCoords(:,1), knotCoords(:,2), 'ks', 'MarkerSize', 6)
hold off
axis equal
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

legend('Original geometry','Control points', 'Collocation points', 'Element edges')


end

