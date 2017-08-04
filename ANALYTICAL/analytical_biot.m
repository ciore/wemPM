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
% analytical solution to 1D Biot

function soln=analytical_biot(freq,x)

nx=1;
ny=0;
L=max(x)-min(x);

%%
addpath('PLANES')
air=air_properties_generic;
medium=4003;
PEM.name_mat=['Mat_porous_' num2str(medium-1000*floor(medium/1000))];
PEM.typ_mat=floor(medium/1000);
eval(['PEM=Mat_porous_' num2str(medium-1000*floor(medium/1000)),'(PEM);'])
PEM=properties_JCA(PEM,air,freq);
PEM=properties_PEM(PEM,air,freq);
[delta,mu]=compute_Biot_waves(PEM,freq);

%%
% 3 forward waves
SV=Phi_Biot_vector(nx,ny,delta(1),delta(2),delta(3),mu(1),mu(2),mu(3),PEM.N,PEM.A_hat,PEM.K_eq_til,freq);
% 3 backward waves
SV=[SV Phi_Biot_vector(-nx,-ny,delta(1),delta(2),delta(3),mu(1),mu(2),mu(3),PEM.N,PEM.A_hat,PEM.K_eq_til,freq)];

% Cancelation of unnecessary fields for this 1D problem  (v_y^s,v_y^t,sigma_xy and sigma_yy)
SV([2 4 6 7],:)=[];
% Cancelation of Biot wave #3 which is not excited
SV(:,[3 6])=[];

% extraction of equations in x=0;
temp=SV*diag([exp(-1i*delta(1)*0) exp(-1i*delta(2)*0) exp(1i*delta(1)*0) exp(1i*delta(2)*0)]);
Mat(1:2,1:4)=temp([3 4],:);
RHS(1:2,1)=[0;1];

% extraction of equations in x=L;
temp=SV*diag([exp(-1i*delta(1)*L) exp(-1i*delta(2)*L) exp(1i*delta(1)*L) exp(1i*delta(2)*L)]);
Mat(3:4,1:4)=temp([1 2],:);
RHS(3:4,1)=[0;0];

% Resolution of the linear system to obtain the wave amplitudes
X=Mat\RHS;

%% Solution
soln.us=SV(1,1)*X(1)*exp(-1i*delta(1)*x)+SV(1,2)*X(2)*exp(-1i*delta(2)*x)+SV(1,3)*X(3)*exp(1i*delta(1)*x)+SV(1,4)*X(4)*exp(1i*delta(2)*x);
soln.ut=SV(2,1)*X(1)*exp(-1i*delta(1)*x)+SV(2,2)*X(2)*exp(-1i*delta(2)*x)+SV(2,3)*X(3)*exp(1i*delta(1)*x)+SV(2,4)*X(4)*exp(1i*delta(2)*x);
soln.s=SV(3,1)*X(1)*exp(-1i*delta(1)*x)+SV(3,2)*X(2)*exp(-1i*delta(2)*x)+SV(3,3)*X(3)*exp(1i*delta(1)*x)+SV(3,4)*X(4)*exp(1i*delta(2)*x);
soln.p=SV(4,1)*X(1)*exp(-1i*delta(1)*x)+SV(4,2)*X(2)*exp(-1i*delta(2)*x)+SV(4,3)*X(3)*exp(1i*delta(1)*x)+SV(4,4)*X(4)*exp(1i*delta(2)*x);

soln.dpdx=SV(4,1)*X(1)*-1i*delta(1)*exp(-1i*delta(1)*x)+SV(4,2)*X(2)*-1i*delta(2)*exp(-1i*delta(2)*x)+SV(4,3)*X(3)*1i*delta(1)*exp(1i*delta(1)*x)+SV(4,4)*X(4)*1i*delta(2)*exp(1i*delta(2)*x);
soln.Res=abs(mean(1i*2*pi*freq*PEM.rho_eq_til*(PEM.gamma_til*soln.us+soln.ut)+soln.dpdx));
