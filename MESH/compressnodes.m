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
% remove nodes which have no connectivity and renumber elems
% 
% usage:
% [nodes,elems]=compressnodes(nodes,elems)

function [nodes,elems]=compressnodes(nodes,elems)

n(:,1)=unique(elems);
d=zeros(size(nodes,1),1);
d(n)=(1:length(n))'-n;
elems=elems(:,:)+reshape(d(elems(:,:)),size(elems,1),size(elems,2));
nodes=nodes(n,:);
