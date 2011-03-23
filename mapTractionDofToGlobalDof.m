function [ globalTractionDof ] = mapTractionDofToGlobalDof( tractionDof, tracConn, dispConn )
% Since tractions are defined over elements, it is necessary at some point
% to map them to the global DOF which correspond to the nodes. 

% What this function does is take a vector of Traction DOF and outputs a
% vector of the same size but with the global DOF in place

globalTractionDof=zeros(1,length(tractionDof));

elements=ceil(tractionDof/6);   % the elements of each of Dof

for index=1:length(tractionDof)
    value=tractionDof(index);
    localNode=ceil((value - (elements(index)-1)*6)/2);  % the local node number
    dir=value - ((elements(index)-1)*6 + (localNode-1)*2);  % the direction (x=1,y=2)
    
    globalNode=dispConn(elements(index),localNode); % the globalNode number
    globalDof=(globalNode-1)*2 + dir;   % and the DOF we want
    
    globalTractionDof(index)=globalDof;
    
end

end

