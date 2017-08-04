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
% add a boundary condition to a face

function bcs=addbcface(bcs,g,facei,numi)

bcs.face(g)={[bcs.face{g}; repmat(numi,numel(find(bcs.face{g}==facei)),1)]};
fnames=fieldnames(bcs);
for f=1:numel(fnames)
  if iscell(eval(['bcs.',fnames{f}]))&not(strcmp(['bcs.',fnames{f}],'bcs.face'))
    eval(['bcs.',fnames{f},'(g)={[bcs.',fnames{f},'{g}; bcs.',fnames{f},'{g}(bcs.face{g}==facei,:)]};']);
  end
end
bcs.nbnd(g)=numel(bcs.ndsi{g});

