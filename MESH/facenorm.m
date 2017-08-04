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
% evaluate face norms
%
% usage:
% facesnorm=facenorm(nodes,faces)

function facesnorm=facenorm(nodes,faces)

ndim=log(size(faces,2))/log(2);
v1=nodes(faces(:,2),:)-nodes(faces(:,1),:);
if ndim==1
  v2=ones(size(v1,1),1)*[0 0 1];
elseif ndim==2
  v2=nodes(faces(:,3),:)-nodes(faces(:,1),:);
end
facesnorm=cross(v1,v2,2);
facesnorm=facesnorm./(sqrt(dot(facesnorm,facesnorm,2))*ones(1,size(v2,2)));
