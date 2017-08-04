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
% analytical to equivalent fluid model

function [pair,pefm]=analytical(xair,xefm,freq,spd0,rho0,sigma,d1,anec,scale)

xair=xair-xefm(end);
xefm=xefm-xefm(end);
[Zefm,kefm,~,~,R]=efm(freq,spd0,rho0,sigma,d1,0,anec); %equivalent fluid wave number
if anec==1
  Rend=0;
else
  Rend=1;
end
pefm=exp(-1i*kefm*xefm)+Rend*exp(1i*kefm*xefm);
kair=2*pi*freq/spd0;
pair=exp(-1i*kair*(xair-xair(end)))+R*exp(1i*kair*(xair-xair(end)));
A=pefm(1)/pair(end);
pair=A*pair;
pefm=pefm*scale/pair(1);
pair=pair*scale/pair(1);