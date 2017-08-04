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
% interpolate data from nodes to element centroids
% ONLY WORKS FOR HOMOGENEOUS ELEMENT TYPES
% 
% usage:
% data=node2elem(elems,data)

function data=node2elem(elems,data)

elemstype=elemtype(elems);
for tels=unique(elemstype)'
  if tels==1 %line
    data=(data(elems(:,1),:)+data(elems(:,2),:))/2;
  elseif tels==2 %quadrilateral
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:)+data(elems(:,4),:))/4;
  elseif tels==3 %triangle
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:))/3;
  elseif tels==4 %block
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:)+data(elems(:,4),:)+data(elems(:,5),:)+data(elems(:,6),:)+data(elems(:,7),:)+data(elems(:,8),:))/8;
  elseif tels==5 %wedge
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:)+data(elems(:,4),:)+data(elems(:,5),:)+data(elems(:,7),:))/6;
  elseif tels==6 %tetrahedral
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:)+data(elems(:,5),:))/4;
  elseif tels==7 %pyramid
    data=(data(elems(:,1),:)+data(elems(:,2),:)+data(elems(:,3),:)+data(elems(:,4),:)+data(elems(:,5),:))/5;
  end
end
