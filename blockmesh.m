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
% generate a regular (mutli-)block mesh
%
% usage:
% mesh=blockmesh(nd,dim,os,facen)
% nd is the dimension of the grid
% dim is the dimension of the volume
% os is the offset from +x +y +z
% facen is the boundary number given as [x- x+ y- y+ z- z+]

function mesh=blockmesh(ndT,dimT,osT,facenT)

for g=1:size(dimT,1)
  
  dim=dimT(g,:);
  os=osT(g,:);
  nd=ndT(g,:);
  facen=facenT(g,:);
  
  % prelims
  ndim=size(dim,2);
  dim=[dim zeros(1,3-ndim)];
  os=[os zeros(1,3-ndim)];
  nd=[nd ones(1,3-ndim)];
  facen=[facen zeros(1,2*(3-ndim))];
  
  % generate nodes
  for i=1:nd(1)
    for j=1:nd(2)
      for k=1:nd(3)
        if (nd(1)>1)
          nodes(i,j,k,1)=dim(1)/(nd(1)-1)*(i-1)+os(1);
        else
          nodes(i,j,k,1)=0;
        end
        if (nd(2)>1)
          nodes(i,j,k,2)=dim(2)/(nd(2)-1)*(j-1)+os(2);
        else
          nodes(i,j,k,2)=0;
        end
        if (nd(3)>1)
          nodes(i,j,k,3)=dim(3)/(nd(3)-1)*(k-1)+os(3);
        else
          nodes(i,j,k,3)=0;
        end
      end
    end
  end
  nodes=[reshape(nodes(:,:,:,1),prod(nd),1) reshape(nodes(:,:,:,2),prod(nd),1) reshape(nodes(:,:,:,3),prod(nd),1)];
  clear i j k
  
  % generate elements
  [ix,iy,iz]=meshgrid(1:nd(1),1:nd(2),1:nd(3));
  ix=permute(ix,[2 1 3]);
  iy=permute(iy,[2 1 3]);
  iz=permute(iz,[2 1 3]);
  inodes=[reshape(ix,prod(nd),1) reshape(iy,prod(nd),1) reshape(iz,prod(nd),1)];
  if ndim==1
    i=(1:prod(nd));
    i=setdiff(i,find(inodes(:,1)==nd(1)));
    elems=[1 2];
    i=(i-1)'*ones(1,length(elems));
    elems=ones((nd(1)-1),1)*elems;
    elems=elems+i;
  elseif ndim==2
    i=(1:prod(nd));
    i=setdiff(i,find(inodes(:,1)==nd(1)));
    i=setdiff(i,find(inodes(:,2)==nd(2)));
    elems=[1 2 nd(1)+2 nd(1)+1];
    i=(i-1)'*ones(1,length(elems));
    elems=ones((nd(1)-1)*(nd(2)-1),1)*elems;
    elems=elems+i;
  elseif ndim==3
    i=(1:prod(nd));
    i=setdiff(i,find(inodes(:,1)==nd(1)));
    i=setdiff(i,find(inodes(:,2)==nd(2)));
    i=setdiff(i,find(inodes(:,3)==nd(3)));
    elems=[1 2 nd(1)+2 nd(1)+1 nd(1)*nd(2)+1 nd(1)*nd(2)+2 nd(1)*nd(2)+nd(1)+2 nd(1)*nd(2)+nd(1)+1 ];
    i=(i-1)'*ones(1,length(elems));
    elems=ones((nd(1)-1)*(nd(2)-1)*(nd(3)-1),1)*elems;
    elems=elems+i;
  end
  clear i ix iy iz inodes
  
  % generate bcels
  ndim=length(find(nd>1));
  if ndim==1
    face(1,1)={[1 1 1]};
    face(2,1)={[(nd(1)-1) 1 2]};
  elseif ndim==2
    face(1,1)={[(1:(nd(1)-1):(nd(1)-1)*(nd(2)-1))' 2*ones(nd(2)-1,1) 4*ones(nd(2)-1,1)]};
    face(2,1)={[(1:(nd(1)-1):(nd(1)-1)*(nd(2)-1))'+(nd(1)-2) 2*ones(nd(2)-1,1) 2*ones(nd(2)-1,1)]};
    face(1,2)={[(1:nd(1)-1)' 2*ones(nd(1)-1,1) 1*ones(nd(1)-1,1)]};
    face(2,2)={[(1:nd(1)-1)'+(nd(1)-1)*(nd(2)-2) 2*ones(nd(1)-1,1) 3*ones(nd(1)-1,1)]};
  elseif ndim==3
    face(1,1)={[reshape([(1:(nd(1)-1):(nd(1)-1)*(nd(2)-1))'*ones(1,nd(3)-1)+ones(nd(2)-1,1)*(1:(nd(2)-1)*(nd(1)-1):(nd(2)-1)*(nd(1)-1)*(nd(3)-1))-1],(nd(2)-1)*(nd(3)-1),1) 4*ones((nd(2)-1)*(nd(3)-1),1) 4*ones((nd(2)-1)*(nd(3)-1),1)]};
    face(2,1)={[reshape([(1:(nd(1)-1):(nd(1)-1)*(nd(2)-1))'*ones(1,nd(3)-1)+ones(nd(2)-1,1)*(1:(nd(2)-1)*(nd(1)-1):(nd(2)-1)*(nd(1)-1)*(nd(3)-1))-1],(nd(2)-1)*(nd(3)-1),1)+(nd(1)-2) 4*ones((nd(2)-1)*(nd(3)-1),1) 2*ones((nd(2)-1)*(nd(3)-1),1)]};
    face(1,2)={[reshape([(1:nd(1)-1)'*ones(1,nd(3)-1)+ones(nd(1)-1,1)*(1:(nd(2)-1)*(nd(1)-1):(nd(2)-1)*(nd(1)-1)*(nd(3)-1))-1],(nd(1)-1)*(nd(3)-1),1) 4*ones((nd(1)-1)*(nd(3)-1),1) 1*ones((nd(1)-1)*(nd(3)-1),1)]};
    face(2,2)={[reshape([(1:nd(1)-1)'*ones(1,nd(3)-1)+ones(nd(1)-1,1)*(1:(nd(2)-1)*(nd(1)-1):(nd(2)-1)*(nd(1)-1)*(nd(3)-1))-1],(nd(1)-1)*(nd(3)-1),1)+(nd(2)-2)*(nd(1)-1) 4*ones((nd(1)-1)*(nd(3)-1),1) 3*ones((nd(1)-1)*(nd(3)-1),1)]};
    face(1,3)={[(1:(nd(2)-1)*(nd(1)-1))' 4*ones((nd(2)-1)*(nd(1)-1),1) 5*ones((nd(2)-1)*(nd(1)-1),1)]};
    face(2,3)={[(1:(nd(2)-1)*(nd(1)-1))'+(nd(2)-1)*(nd(1)-1)*(nd(3)-2) 4*ones((nd(2)-1)*(nd(1)-1),1) 6*ones((nd(2)-1)*(nd(1)-1),1)]};
  end
  faces=[];
  facen=facen(:,1:ndim*2);
  for ibc=unique(facen)
    for i=find(ismember(facen,ibc))'
      for j=1:numel(i)
        faces=[faces;face{i(j)} ibc*ones(size(face{i(j)},1),1)];
      end
    end
  end
  
  % add to structure
  mesh.nodes(g)={nodes(:,1:ndim)};
  mesh.elems(g)={elems};
  mesh.faces(g)={faces};
  mesh.nnds(g)=size(nodes,1);
  mesh.nels(g)=size(elems,1);
  mesh.nfcs(g)=size(faces,1);
  mesh.ndim=ndim;
  mesh.ngrd=size(dimT,1);
  clear nodes elems faces
  
end