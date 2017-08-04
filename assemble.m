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
% assemble the stiffness matrix and forcing vector

function [stiff,force,pinverr,tol,condH]=assemble(mesh,bcs,physics,freq)

fprintf('assembling...\n');
operationtime=cputime;

%% define wave directions
ndim=mesh.ndim;
if ndim==1
  nwav=2;
  t=0:pi:pi;
  t=[cos(t)' sin(t)' zeros(size(t))'];
elseif ndim==2
  nwav=11;
  t=2*pi/nwav:2*pi/nwav:2*pi;
  t=[cos(t)' sin(t)' zeros(size(t))'];
elseif ndim==3
  [~,t]=bucky;
  nwav=size(t,1);
end

%% assemble stiffness matrix
fprintf(' looping through assembly ... ');
nnds=mesh.nnds;
ngrd=mesh.ngrd;
stiff=[];
force=[];
ndof=8;

%% loop through grids
for g=1:ngrd
  nodes=mesh.nodes{g};
  elems=mesh.elems{g};
  bndsi=bcs.ndsi{g};
  bcoef=bcs.coef{g};
  brhds=bcs.rhds{g};
  bnorm=bcs.norm{g};
  robi=find(bcs.type{g}==2);
%   diri=find(bcs.type{g}==1);
%   inti=find(bcs.type{g}==9);
  
  %% loop through point
  for i=1:size(nodes,1)
    
    %% compute Biot waves
    addpath('PLANES')
    [delta,mu]=compute_Biot_waves(physics,freq);
    k=[delta(1)*ones(1,nwav) delta(2)*ones(1,nwav) delta(3)*ones(1,nwav)];
    tta=[t; t; t];
    for ti=1:size(t,1)
      vec(:,ti:size(t,1):3*size(t,1))=Phi_Biot_vector(t(ti,1),t(ti,2),delta(1),delta(2),delta(3),mu(1),mu(2),mu(3),physics.N,physics.A_hat,physics.K_eq_til,freq);
    end
    
    %% find neighbours
    [ii,~]=find(elems==i);
    iconn=setdiff(unique(nonzeros(elems(ii,:))),i);
    
    %    %internal bc
    %     if sum(ismember(bcint(:,1),i))>0 %discontinuity/internal bc
    %       ii=find(bcint(:,1)==i);
    % %       ii=ii(1);
    %       iconn=[iconn; bcint(ii,5)];
    %     end

    x0=nodes(i,:)-nodes(i,:); %local coordinates
    xm=nodes(iconn,:)-(ones(length(iconn),1)*nodes(i,:)); %local coordinates


    h0=exp(-1i*(x0(:,1)*(k.*tta(:,1).')+x0(:,2)*(k.*tta(:,2).')));
    H0=exp(-1i*(xm(:,1)*(k.*tta(:,1).')+xm(:,2)*(k.*tta(:,2).')));
    h=[];
    H=[];
    for ii=1:ndof
      h=[h; vec(ii,:).*h0];
      H=[H; repmat(vec(ii,:),numel(iconn),1).*H0];
    end
    if sum(ismember(bndsi(robi),i))>0
      ii=robi(bndsi(robi)==i);
      for iii=ii'
        H=[H; vec(find(bcoef(iii,:)),:).*h0];
      end
    end
    Hplus=pinv(H);%,1e-3);
    pinverr(i,:)=sum(sum(abs(H*Hplus*H-H).^2,1),2);
    tol(i,:)=max(size(H))*eps(norm(H));
    condH(i,:)=norm(H)*norm(Hplus);
    hHplusL=h*Hplus(:,1:length(iconn)*ndof);
    hHplusR=h*Hplus(:,length(iconn)*ndof+1:end);
    stiff=[stiff; [(i:nnds:nnds*ndof)'; reshape(((i:nnds:nnds*ndof)'*ones(1,length(iconn)*ndof))',length(iconn)*ndof^2,1)]...
      [(i:nnds:nnds*ndof)'; reshape(reshape(iconn*ones(1,ndof)+ones(length(iconn),1)*(nnds*(0:ndof-1)),length(iconn)*(ndof),1)*ones(1,ndof),length(iconn)*(ndof)^2,1)]...
      [ones(ndof,1); reshape(-hHplusL.',numel(hHplusL),1)]...
      ];
    if sum(ismember(bndsi(robi),i))>0
      ii=robi(bndsi(robi)==i);
      force=[force; [(i:nnds:nnds*ndof)']...
        [ones(size(hHplusR,1),1)]...
        [reshape(hHplusR*brhds(ii),numel(hHplusR(:,1)),1)]...
        ];
    end
    
%     %%
%     if sum(ismember(bcint(:,1),i))>0
%       ii=find(bcint(:,1)==i);
%       ii=ii(1);
%       stiff=[stiff; (i:nnds:nnds*ndof)' ((0:ndof-1)'*nnds)+bcint(ii,5) -1*ones(ndof,1)];
%     end
  end
  
end

%%
stiff=sparse(stiff(:,1),stiff(:,2),stiff(:,3));
force=sparse(force(:,1),force(:,2),force(:,3),size(stiff,1),1);

%%
fprintf('done\n');

operationtime=cputime-operationtime;
fprintf(['assembling took ',num2str(operationtime),' seconds\n']);
