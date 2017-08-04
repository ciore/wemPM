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
% definition of element face numbering
% uses fluent/gambit numbering
% 
% usage:
% elemface=defelemsfaces(i,r)
% where
% i = 1 -> line
% i = 2 -> quadrilateral
% i = 3 -> triangle
% i = 4 -> block
% i = 5 -> wedge
% i = 6 -> tetrahedral
% i = 7 -> pyramid
% and
% r returns only the r'th face entry

function elemface=elemfacedef(i,r)

elemface=[{[1; 2]};...
  {[1 2; 2 3; 3 4; 4 1]};...
  {[1 2; 2 3; 3 1]};...
  {[1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 2 1 4 3; 5 6 7 8]};...
  {[1 5 7 4; 5 2 3 7; 2 1 4 3; 1 2 5 5; 4 7 3 3]};...
  {[2 1 3 3; 1 2 5 5; 2 3 5 5; 3 1 5 5]};...
  {[1 4 3 2; 1 2 5 5; 2 3 5 5; 3 4 5 5; 4 1 5 5]}];
elemface=elemface{i};
if exist('r')==1
  elemface=elemface(r,:);
end
