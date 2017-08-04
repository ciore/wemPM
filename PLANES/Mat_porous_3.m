% Mat_porous_3.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
%
% This file is part of PLANES.
%
% PLANES (Porous LAum NumErical Simulator) is a software to compute the
% vibroacoustic response of sound packages containing coupled
% acoustic/elastic/porous substructures. It is mainly based on the
% Finite-Element Method and some numerical methods developped at
% LAUM (http://laum.univ-lemans.fr).
%
% You can download the latest version of PLANES at
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
%
% For any questions or if you want to
% contribute to this project, contact Olivier.
%
% PLANES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%

function PEM=Mat_porous_3(PEM)

%% Porous material model parameters
PEM.eqf='JCA';
PEM.frame='structural';
PEM.aniso='no';
PEM.phi=0.95;
PEM.sig=42000;
PEM.alpha=1.100;
PEM.LCV=1.50E-05;
PEM.LCT=4.500E-05;
PEM.rho_1=126.000;
PEM.young=694400E+00;
PEM.nu=0.24000E+00;
PEM.eta=0.05;
PEM.N=(PEM.young)/(2*(1+PEM.nu));
PEM.A_hat=(PEM.young*PEM.nu)/((1+PEM.nu)*(1-2*PEM.nu));

%%
PEM=orderfields(PEM);

% %%
% fnames=fieldnames(PEM);
% for f=1:numel(fnames)
%   eval([fnames{f},'=PEM.',fnames{f},';']);
% end

