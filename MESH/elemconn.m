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
% determine node to element connectivities
%
%usage:
%conn=elemconn(elems[,iblanks])

function conn=elemconn(conn,iblanks)

nnel=size(conn,2);
[n,i]=sort(reshape(conn',numel(conn),1));
node=unique(n);
conn=ceil(i/nnel);
[~,i]=unique([n,conn],'rows');
n=n(i);
conn=conn(i);
a=find([n(2:end); n(end)+1]-n);
a=a-[min(n)-1; a(1:end-1)];
for i=1:max(a)-1
  ni=find(a==max(a)-i);
  n=[n; reshape(node(ni)*ones(1,i),numel(ni)*i,1)];
  conn=[conn; zeros(numel(ni)*i,1)];
end
[n,i]=sort(n);
conn=conn(i);
conn=reshape(conn,max(a),length(conn)/max(a))';
conn(node,:)=conn;
if exist('iblanks')==1
  conn(find(iblanks==0),:)=0;
end