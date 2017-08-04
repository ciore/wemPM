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
% interpolate data from element centroids to nodes 
% 
% usage:
% data=elem2node(nodes,elems,data)

function data=elem2node(nodes,elems,data)

ielems=elemconn(elems);
elemscen=node2elem(elems,nodes);
relems=zeros(size(ielems));
delems=zeros(size(ielems));
for i=1:size(ielems,2)
  ii=find(ielems(:,i));
  v=nodes(ii,:)-elemscen(ielems(ii,i),:);
  r=sqrt(dot(v,v,2));
  relems(ii,i)=r;
  delems(ii,i)=data(ielems(ii,i),:);
end
warning off MATLAB:divideByZero
relemsinv=1./relems;
relemsinv(find(isinf(relemsinv)))=0;
data=sum(delems.*relemsinv,2)./sum(relemsinv,2);

