% properties_PEM.m
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

function PEM=properties_PEM(PEM,air,freq)

omega=2*pi*freq;

%% Biot densities with tortuosity effects
PEM.rho_12=-PEM.phi*air.rho*(PEM.alpha-1);
PEM.rho_11=PEM.rho_1-PEM.rho_12;
PEM.rho_2=PEM.phi*air.rho;
PEM.rho_22=PEM.rho_2-PEM.rho_12;
PEM.rho_22_til=PEM.phi^2*PEM.rho_eq_til;
PEM.rho_12_til=PEM.rho_2-PEM.rho_22_til;
PEM.rho_11_til=PEM.rho_1-PEM.rho_12_til;
PEM.rho_til=PEM.rho_11_til-((PEM.rho_12_til.^2)./PEM.rho_22_til);

PEM.gamma_til=PEM.phi*(PEM.rho_12_til./PEM.rho_22_til-(1-PEM.phi)/PEM.phi);
PEM.rho_s_til=PEM.rho_til+PEM.gamma_til.^2.*PEM.rho_eq_til;

%% Biot in-vacuo elastic coefficients
PEM.N=(PEM.young)/(2*(1+PEM.nu));
PEM.A_hat=(PEM.young*PEM.nu)/((1+PEM.nu)*(1-2*PEM.nu));
switch PEM.frame
  case{'anelastic'}
    PEM.structural_loss=1+(b_hat*(1i*omega/beta_hat)^alpha_hat)/(1+(1i*omega/beta_hat)^alpha_hat);
  case{'structural'}
    PEM.structural_loss=1+1i*PEM.eta;
end
PEM.N=PEM.N*PEM.structural_loss;
PEM.A_hat=PEM.A_hat*PEM.structural_loss;
PEM.P_hat=PEM.A_hat+2*PEM.N;

%% Biot 1956 elastic coefficients
PEM.R_til=PEM.K_eq_til*PEM.phi^2;
PEM.Q_til=((1-PEM.phi)/PEM.phi)*PEM.R_til;
PEM.P_til=PEM.P_hat+PEM.Q_til.^2./PEM.R_til;

%% Parameters for energies
PEM.xi=(1-PEM.phi)/PEM.phi;
PEM.b_til=(1i*omega)*(PEM.rho_12-PEM.rho_12_til);
PEM.b_r=real(PEM.b_til);
PEM.rho_f_til=PEM.rho_til-PEM.rho_1;

%%
if strcmp(PEM.aniso,'yes')
  PEM.C_hat=C_hat_conservative*PEM.structural_loss;
end

%%
PEM=orderfields(PEM);
