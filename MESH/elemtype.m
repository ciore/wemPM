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
% definition of element types
% uses fluent/gambit numbering, see elemfacedef 
% 
% usage:
% elemstype=elemtype(elems)
% where elems is a n x 4 matrix for n 2d elements and a n x 8 matrix for n 3d elements 
% element types are:
% 1 -> line
% 2 -> quadrilateral
% 3 -> triangle
% 4 -> block
% 5 -> wedge
% 6 -> tetrahedral
% 7 -> pyramid

function elemstype=elemtype(elems)

if size(elems,2)==2
  elemstype=1*ones(size(elems,1),1);
elseif size(elems,2)==4
  elemstype=2*ones(size(elems,1),1);
  elemstype(find(elems(:,4)==elems(:,3)))=3;
else
  elemstype=4*ones(size(elems,1),1);
  elemstype(find((elems(:,4)~=elems(:,3))&(elems(:,6)==elems(:,5))&(elems(:,7)~=elems(:,6))))=5;
  elemstype(find((elems(:,4)==elems(:,3))&(elems(:,6)==elems(:,5))))=6;
  elemstype(find((elems(:,4)~=elems(:,3))&(elems(:,6)==elems(:,5))&(elems(:,7)==elems(:,6))))=7;
end

