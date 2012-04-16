function [ presDispDOFs, presTracDOFs, dirichletVals, nonZeroXTracDOFs, nonZeroYTracDOFs ] = assignDirichletAndNeumannNodes(n, elConn, tracConn)

% we work out the nodes where disp is specified and trac specified

% n: number of elements per line
% elConn: the displacement element connectivity matrix
% tracConn: the traction element Connectivity matrix
% traction: the traction that is applied on the upper and left edges
% X- Displacements are specified on line 5 while y-displacements are
% specified on line 1

%% Hole in infinite plate problem

% % y-disps prescribed on first line
% elms=1:n;    
% sctrY=elConn(elms,:);
% sctrY=unique(reshape(sctrY,1,numel(sctrY)));
% 
% % x-disps
% elms=(2*n+1):(3*n);
% sctrX=elConn(elms,:);
% sctrX=unique(reshape(sctrX,1,numel(sctrX)));
% 
% presDispDOFs=sort([sctrX*2-1 sctrY*2]);
% dirichletVals=zeros(1,length(presDispDOFs));
% 
% 
% % and now for the tractions
% 
% % x-tracs
% elms=[1:2*n (3*n+1):5*n];
% sctrX=tracConn(elms,:);
% sctrX=unique(reshape(sctrX,1,numel(sctrX)));
% 
% % y-tracs
% elms=(n+1):(5*n);
% sctrY=tracConn(elms,:);
% sctrY=unique(reshape(sctrY,1,numel(sctrY)));
% 
% presTracDOFs=sort([sctrX*2-1 sctrY*2]);
% 
% % work out the traction DOFs where we specify a nonzero traction
% 
% % nonzero y-tractions
% elms=(3*n+1):4*n;
% sctrY=tracConn(elms,:);
% sctrY=reshape(sctrY,1,numel(sctrY));
% nonZeroYTracDOFs=unique(sctrY*2);
% 
% % nonzero x-tractions
% elms=(4*n+1):5*n;
% sctrX=tracConn(elms,:);
% sctrX=unique(reshape(sctrX,1,numel(sctrX)));
% nonZeroXTracDOFs=sctrX*2-1;

%% L-plate problem

% The only displacements that we prescribe are the nodes in the notch and
% at the right hand corner

% elms=[4]*n + 1;
% sctrX = elConn(elms,1);
% sctrX=reshape(sctrX,1,numel(sctrX));
% 
% elms=[1 4]*n + 1;
% sctrY = elConn(elms,1);
% sctrY=reshape(sctrY,1,numel(sctrY));
% 
% presDispDOFs=sort([sctrX*2-1 sctrY*2]);
% dirichletVals=zeros(1,length(presDispDOFs));
% 
% 
% % Now let's prescribe our tractions - this is everywhere except the two
% % dirichlet nodes
% tracEndElms = [4]*n;
% tracStartElms = tracEndElms + 1;
% sctrX = [tracConn(tracEndElms,end) tracConn(tracStartElms,1)];
% sctrX = reshape(sctrX,1, numel(sctrX));
% 
% tracEndElms = [1 4]*n;
% tracStartElms = tracEndElms + 1;
% sctrY = [tracConn(tracEndElms,end) tracConn(tracStartElms,1)];
% sctrY = reshape(sctrY,1, numel(sctrY));
% 
% allTracDOFs = 1:2*max(tracConn(:,end));
% presTracDOFs = setxor(allTracDOFs, [sctrX*2-1 sctrY*2]);
% 
% % and our non-zero trac DOF 
% 
% % and now get the nonzero x and y tractions
% nonZeroYTracDOFs = presTracDOFs(find(mod(presTracDOFs,2)==0));
% nonZeroXTracDOFs = presTracDOFs(find(mod(presTracDOFs,2)==1));

% nonzero tractions
% elms=[1:3*n 5*n+1:6*n];
% sctr=tracConn(elms,:);
% sctr=reshape(sctr,1,numel(sctr));
% nonZeroYTracDOFs=unique(sctr*2);
% nonZeroXTracDOFs=unique(sctr*2-1);

%% spanner example

% prescribe our Dirichlet nodes

% elms=[1 12] *n + 1;
% sctrX = elConn(elms,1);
% sctrX=reshape(sctrX,1,numel(sctrX));
% 
% elms=[1 2 11 12]*n + 1;
% sctrY = elConn(elms,1);
% sctrY=reshape(sctrY,1,numel(sctrY));
% 
% presDispDOFs=sort([sctrX*2-1 sctrY*2]);
% dirichletVals=zeros(1,length(presDispDOFs));
% 
% % Now let's prescribe our tractions - this is everywhere except the two
% % dirichlet nodes
% tracEndElms = [1 12]*n;
% tracStartElms = tracEndElms + 1;
% sctrX = [tracConn(tracEndElms,end) tracConn(tracStartElms,1)];
% sctrX = reshape(sctrX,1, numel(sctrX));
% 
% tracEndElms = [1 2 11 12]*n;
% tracStartElms = tracEndElms + 1;
% sctrY = [tracConn(tracEndElms,end) tracConn(tracStartElms,1)];
% sctrY = reshape(sctrY,1, numel(sctrY));
% 
% allTracDOFs = 1:2*max(tracConn(:,end));
% presTracDOFs = setxor(allTracDOFs, [sctrX*2-1 sctrY*2]);
% 
% % and our non-zero trac DOF 
% prevElmt=6*n;
% tracElms = 6*n+1:7*n;
% procElmt=7*n+1;
% rightSideDOF = unique(reshape(tracConn(tracElms,:),1,numel(tracConn(tracElms,:))));
% sctrY = [tracConn(prevElmt,end) rightSideDOF tracConn(procElmt,1)];
% 
% % and now get the nonzero x and y tractions
% nonZeroYTracDOFs = sctrY*2;
% nonZeroXTracDOFs = [];

%% uniaxial tension

% first prescribe out dirichlet nodes
elms=1:n;    
sctrY=elConn(elms,:);
sctrY=unique(reshape(sctrY,1,numel(sctrY)));

% x-disps
elms=3*n+1:4*n;
sctrX=elConn(elms,:);
sctrX=unique(reshape(sctrX,1,numel(sctrX)));

presDispDOFs=sort([sctrX*2-1 sctrY*2]);
dirichletVals=zeros(1,length(presDispDOFs));


% and now for the tractions

% x-tracs
elms=1:3*n;
sctrX=tracConn(elms,:);
sctrX=unique(reshape(sctrX,1,numel(sctrX)));

%firstNode=1;
%sctrX = setxor(sctrX,firstNode);

% y-tracs
elms=(n+1):(4*n);
sctrY=tracConn(elms,:);
sctrY=unique(reshape(sctrY,1,numel(sctrY)));

%lastNode=tracConn(end,end);
%sctrY = setxor(sctrY,lastNode);

presTracDOFs=sort([sctrX*2-1 sctrY*2]);

% work out the traction DOFs where we specify a nonzero traction

% nonzero y-tractions
elms=(2*n+1):3*n;
sctrY=tracConn(elms,:);
sctrY=reshape(sctrY,1,numel(sctrY));
nonZeroYTracDOFs=unique(sctrY*2);

% nonzero x-tractions
nonZeroXTracDOFs=[];


end

