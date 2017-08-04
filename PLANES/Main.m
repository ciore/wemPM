clear

freq=1500;
x=linspace(0,0.05,201);

nx=1;
ny=0;
L=max(x)-min(x);
% omega=2*pi*freq;

%%
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

SV(4,:)=-SV(4,:); %%Added by Ciaran

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
us=SV(1,1)*X(1)*exp(-1i*delta(1)*x)+SV(1,2)*X(2)*exp(-1i*delta(2)*x)+SV(1,3)*X(3)*exp(1i*delta(1)*x)+SV(1,4)*X(4)*exp(1i*delta(2)*x);
ut=SV(2,1)*X(1)*exp(-1i*delta(1)*x)+SV(2,2)*X(2)*exp(-1i*delta(2)*x)+SV(2,3)*X(3)*exp(1i*delta(1)*x)+SV(2,4)*X(4)*exp(1i*delta(2)*x);
s=SV(3,1)*X(1)*exp(-1i*delta(1)*x)+SV(3,2)*X(2)*exp(-1i*delta(2)*x)+SV(3,3)*X(3)*exp(1i*delta(1)*x)+SV(3,4)*X(4)*exp(1i*delta(2)*x);
p=SV(4,1)*X(1)*exp(-1i*delta(1)*x)+SV(4,2)*X(2)*exp(-1i*delta(2)*x)+SV(4,3)*X(3)*exp(1i*delta(1)*x)+SV(4,4)*X(4)*exp(1i*delta(2)*x);

dpdx=SV(4,1)*X(1)*-1i*delta(1)*exp(-1i*delta(1)*x)+SV(4,2)*X(2)*-1i*delta(2)*exp(-1i*delta(2)*x)+SV(4,3)*X(3)*1i*delta(1)*exp(1i*delta(1)*x)+SV(4,4)*X(4)*1i*delta(2)*exp(1i*delta(2)*x);
Res=abs(mean(1i*2*pi*freq*PEM.rho_eq_til*(PEM.gamma_til*us+ut)+dpdx))

%% plots

figure(5)
clf

subplot(2,2,1)
plot(x,real(us),x,imag(us))
grid
ylabel('us')

subplot(2,2,2)
plot(x,real(ut),x,imag(ut))
grid
ylabel('ut')

subplot(2,2,3)
plot(x,real(s),x,imag(s))
grid
ylabel('sig')
xlabel('x')

subplot(2,2,4)
plot(x,real(p),x,imag(p))
grid
ylabel('p')
xlabel('x')