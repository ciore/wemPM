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
clear

%% inputs
freq=10^3;%*(2.^([-20:20]/3)); %frequency

%% setup domain and mesh
domain.dim=[0.05 0.025]; %dimensions of domain
domain.off=[0 0]; %offset of domain
domain.fce=[1 2 3 4]; %face ids (for boundary conditions) on [-x +x -y +y]
nnd=domain.dim./0.0025+1; %number of nodes in each direction
mesh=blockmesh(nnd,domain.dim,domain.off,domain.fce);

%% loop for frequency
for f=1:length(freq)
  fprintf(['computing frequency ',num2str(f),' of ',num2str(length(freq)),' ... \n']);
  
  %% setup boundary conditions
  addpath('MESH')
  bcs=initbcs(mesh.nodes,mesh.elems,mesh.faces,8);
  bcs=addbcface(bcs,1,1,5);
  bcs=addbcface(bcs,1,2,6);
  bcs=addbcface(bcs,1,3,7);
  bcs=addbcface(bcs,1,4,8);
  bcs=addbcface(bcs,1,1,9);
%   bcs=addbcface(bcs,1,2,10);
  bcs=addbcface(bcs,1,3,11);
  bcs=addbcface(bcs,1,4,12);
  bcs.type{1}(:)=2;
  bcs.coef{1}(bcs.face{1}==1,5)=1; %face 1
  bcs.coef{1}(bcs.face{1}==2,1)=1; %face 2
  bcs.coef{1}(bcs.face{1}==3,2)=1; %face 3
  bcs.coef{1}(bcs.face{1}==4,2)=1; %face 4
  bcs.coef{1}(bcs.face{1}==5,6)=1; %face 5
  bcs.coef{1}(bcs.face{1}==6,3)=1; %face 6
  bcs.coef{1}(bcs.face{1}==7,4)=1; %face 7
  bcs.coef{1}(bcs.face{1}==8,4)=1; %face 8
  bcs.coef{1}(bcs.face{1}==9,8)=1; bcs.rhds{1}(bcs.face{1}==9)=1; %face 9
%   bcs.coef{1}(bcs.face{1}==10,1)=1; %face 10
  bcs.coef{1}(bcs.face{1}==11,6)=1; %face 11
  bcs.coef{1}(bcs.face{1}==12,6)=1; %face 12
  
  %% setup physics
  addpath('PLANES')
  air=air_properties_generic;
  medium=4003;
  PEM.name_mat=['Mat_porous_' num2str(medium-1000*floor(medium/1000))];
  PEM.typ_mat=floor(medium/1000);
  eval(['PEM=Mat_porous_' num2str(medium-1000*floor(medium/1000)),'(PEM);'])
  PEM=properties_JCA(PEM,air,freq(f));
  PEM=properties_PEM(PEM,air,freq(f));
  physics=PEM;

  %% assemble
  [stiff,force,err,tol,condH]=assemble(mesh,bcs,physics,freq(f));
   
  %%  solve
  fprintf('solving ... \n');
  operationtime=cputime;
  fprintf(' solving directly ... ');
  q=full(stiff\force);
  fprintf('done\n');
  
  %% postprocess
  nnds=mesh.nnds;
  soln.usx(:,f)=q((1:nnds)+nnds*0);
  soln.usy(:,f)=q((1:nnds)+nnds*1);
  soln.utx(:,f)=q((1:nnds)+nnds*2);
  soln.uty(:,f)=q((1:nnds)+nnds*3);
  soln.sxx(:,f)=q((1:nnds)+nnds*4);
  soln.sxy(:,f)=q((1:nnds)+nnds*5);
  soln.syy(:,f)=q((1:nnds)+nnds*6);
  soln.p(:,f)=q((1:nnds)+nnds*7);
%   soln.cond(:,f)=condest(stiff);
  
  us=mean(reshape(soln.usx(:,f),nnd(1),nnd(2)),2);
  ut=mean(reshape(soln.utx(:,f),nnd(1),nnd(2)),2);
  s=mean(reshape(soln.sxx(:,f),nnd(1),nnd(2)),2);
  p=mean(reshape(soln.p(:,f),nnd(1),nnd(2)),2);
  
  Z=p(1)/ut(1);
  R=(Z-air.rho*air.c)/(Z+air.rho*air.c);
  alpha(:,f)=1-abs(R)^2;
  
  x=mesh.nodes{1}(1:nnd(1),1);
  addpath('ANALYTICAL')
  ana=analytical_biot(freq(f),x);
  x0=linspace(min(x),max(x),201);
  ana0=analytical_biot(freq(f),x0);
  Z=ana.p(1)/ana.ut(1);
  R=(Z-air.rho*air.c)/(Z+air.rho*air.c);
  alphaana(:,f)=1-abs(R)^2;  
  
  Err(1,f)=sqrt(sum(abs(us-ana.us).^2))./sqrt(sum(abs(ana.us.^2)));
  Err(2,f)=sqrt(sum(abs(ut-ana.ut).^2))./sqrt(sum(abs(ana.ut.^2)));
  Err(3,f)=sqrt(sum(abs(s-ana.s).^2))./sqrt(sum(abs(ana.s.^2)));
  Err(4,f)=sqrt(sum(abs(p-ana.p).^2))./sqrt(sum(abs(ana.p.^2)));
  
  PinvRes(:,f)=sqrt(sum(err.^2));
  Tol(:,f)=mean(tol);
  CondH(:,f)=mean(condH);

  
  %% visualise

  figure(1)
  
  subplot(2,2,1)
  plot(x0,real(ana0.us),'b',x0,imag(ana0.us),'r',x,real(us),'ob',x,imag(us),'or')
  xlabel('x [m]'),ylabel('u_s [m/s]')
  title(['Error = ',num2str(Err(1,f))])
  legend({'Re';'Im'},'location','best')

  subplot(2,2,2)
  plot(x0,real(ana0.ut),'b',x0,imag(ana0.ut),'r',x,real(ut),'ob',x,imag(ut),'or')
  xlabel('x [m]'),ylabel('u_t [m/s]')
  title(['Error = ',num2str(Err(1,f))])
  legend({'Re';'Im'},'location','best')
  
  subplot(2,2,3)
  plot(x0,real(ana0.s),'b',x0,imag(ana0.s),'r',x,real(s),'ob',x,imag(s),'or')
  xlabel('x [m]'),ylabel('s [Pa]')
  title(['Error = ',num2str(Err(1,f))])
  legend({'Re';'Im'},'location','best')
  
  subplot(2,2,4)
  plot(x0,real(ana0.p),'b',x0,imag(ana0.p),'r',x,real(p),'ob',x,imag(p),'or')
  xlabel('x [m]'),ylabel('p [Pa]')
  title(['Error = ',num2str(Err(1,f))])
  legend({'Re';'Im'},'location','best')
  
  pause(.1)

end

% figure(2)
% clf
% semilogx(freq,alphaana,'b',freq,alpha,'ob')
% xlabel('frequency [Hz]')
% ylabel('Absorption Coefficient')
% 
% figure(3)
% clf
% loglog(freq,Err)
% xlabel('frequency [Hz]')
% ylabel('Relative Error')
% legend({'us';'ut';'s';'p'},'location','southeast')
% 
% figure(4)
% clf
% loglog(freq,PinvRes,'b')
% xlabel('frequency [Hz]')
% ylabel('Residual of pinv (AA^+A - A)')
% 
% figure(5)
% f=1;
% clf
% subplot(1,2,1)
% dpatch(mesh.nodes{:},mesh.elems{:},real(soln.p(:,f)));
% title('Re(p)')
% shading interp
% colorbar;
% colormap(jet);
% subplot(1,2,2)
% dpatch(mesh.nodes{:},mesh.elems{:},abs(soln.p));
% title('Abs(p)')
% shading interp
% colorbar;
% colormap(jet);
% 
