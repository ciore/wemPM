% This file is part of wemPM, a code to compute the vibroacoustic
% response in a poroelasic material using a wave expansion discretisation
% method
% 
% Copyright (C) 2017 Ciaran O'Reilly <ciaran@kth.se>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%
% initialise boundary conditions structure

function bcs=initbcs(nodesT,elemsT,facesT,ndof)

for g=1:size(nodesT,2)
  nodes=nodesT{g};
  ndim=size(nodes,2);
  elems=elemsT{g};
  faces=facesT{g};
  nodes=[nodes zeros(size(nodes,1),3-size(nodes,2))];
  bcnds=[];
  for i=unique(faces(:,4))'
    ii=find(faces(:,4)==i);
    bcsfaces=elemface(elems(faces(ii,1),:),faces(ii,3));
    bci=unique(bcsfaces);
    if log(size(elems,2))/log(2)==1
      conn=nodeconn(elems);
      bcnds=[bcnds; bci angle((nodes(bci,1)-nodes(nonzeros(unique(conn(bci,:))),1))+1i*(nodes(bci,2)-nodes(nonzeros(unique(conn(bci,:))),2))) angle(sqrt((nodes(bci,1)-nodes(nonzeros(unique(conn(bci,:))),1)).^2+(nodes(bci,2)-nodes(nonzeros(unique(conn(bci,:))),2)).^2)+1i*(nodes(bci,3)-nodes(nonzeros(unique(conn(bci,:))),3))) i*ones(size(bci))];
    else
      bcsnorms=facenorm(nodes,bcsfaces);
      [bcsnodes,bcsfaces]=compressnodes(nodes,bcsfaces);
      bcnds=[bcnds; bci angle(elem2node(bcsnodes,bcsfaces,bcsnorms(:,1))+1i*elem2node(bcsnodes,bcsfaces,bcsnorms(:,2))) angle(sqrt(elem2node(bcsnodes,bcsfaces,bcsnorms(:,1)).^2+elem2node(bcsnodes,bcsfaces,bcsnorms(:,2)).^2)+1i*elem2node(bcsnodes,bcsfaces,bcsnorms(:,3))) i*ones(length(unique(bcsfaces)),1)];
    end
  end  
  bcs.ndsi(g)={bcnds(:,1)};
  bcs.face(g)={bcnds(:,4)};
  bcs.type(g)={zeros(size(bcnds(:,1)))};
  bcs.coef(g)={zeros(size(bcnds(:,1),1),ndof)};
  bcs.rhds(g)={zeros(size(bcnds(:,1)))};
  n=[cos(bcnds(:,2)).*cos(bcnds(:,3)) sin(bcnds(:,2)).*cos(bcnds(:,3)) sin(bcnds(:,3))];
  bcs.norm(g)={n(:,1:ndim)};
  bcs.nbnd(g)=size(bcnds,1);
end