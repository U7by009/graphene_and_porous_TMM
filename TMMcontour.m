function TMMcontour(N,P,wavelength,indices,layers,theta)
%% INSTRUCTIONS.
% This function allows to visualize reflectance contourplots for a
% dielectric system of n-layers. The layer n-1 is considered a porous
% media, that means, the effective index of that leyers is evaluated in
% terms of its porosity and the last layer (n-layer). The 2nd layer is 
% assumed to be defined by N-dielectric (n+ik) graphene monolayers [1]. 
%
% Figure 1 shows the reflectance for TE polarization, figure 2 depicts the
% reflectance for TM polarization and figure 3 represents the airthmetic 
% mean. The method used in this Matlab function is the Matrix Transfer 
% Method [2].
%
% The input variables are:
% 
% * Number of graphene monolayers *N* .
% * Porosity *P* of the porous (end-1) layer.
% * *wavelength* is the operational wavelength in nm.
% * *indices* vector [indices(1) indices(2) indices(3), ... ,
% indices(end)] of the refractive index of each layer, must be noticed that
% the refractive index of the porous media must be written as the
% non-porous material index (P=0).
% * *layers* vector [thickness(1) thickness(2) thickness(3), ... ,
% thickness(end)] of the thickness of each layer (the first and the 
% last layer must be introduced as NaN).
% * *theta* is a vector [theta(1) theta(2) theta(3), ... , theta(end)] of
% different angles of incidence of the TE/TM ligh (in radiands)
%
% example:TMM_contour(4,linspace(0,1,100),632.8,[1.61,2.6+1.6i,2.56,1.333],[NaN,0.34,1000,NaN],linspace(0*pi/180,85*pi/180,1000))

%% Initialize the time counter.
% This may be useful for a comparative between numerical methods.
tic
%% Transfer Matrix Method.
% Applying Transfer Matrix Method for a stratified material.

% The 2nd layer is assumed N-dielectric (n+ik) graphene monolayers.
layers(2)=N*0.34;
% Evaluating porous media effective index with Bruggerman formula. 
for q2=1:length(P)
Indices(q2,:)=indices;
Bruggerman=@(x) (1-P(q2))*((Indices(1,end-1)^2-x.^2)/(Indices(1,end-1)^2+2*x.^2))+P(q2)*((Indices(1,end)^2-x.^2)/(Indices(1,end)^2+2*x.^2));
Indices(q2,end-1)=fzero(Bruggerman,[Indices(1,end-1) Indices(1,end)]);

% Defining wavenumbers, dielectric permitivity and dielectric permeability.
for m=1:length(theta)
for n=1:size(Indices,2)
kx=2*pi*Indices(q2,1)*sin(theta(m))./wavelength;
kz(n)=(2*pi/wavelength)*sqrt(Indices(q2,n).^2-Indices(q2,1).^2*sin(theta(m)).^2);
mu(n)=1;
ep(q2,n)=Indices(q2,n).^2;
end

% Defining coefficients.
for n=1:size((Indices),2)-1
a(n)=(kz(n+1)./kz(n))*mu(n)./mu(n+1);
b(n)=(kz(n+1)./kz(n))*ep(q2,n)./ep(q2,n+1);
end
% Transfer matrix from 1-to-2 and N-1-to-N.
M1TE=0.5*[1+a(1),(1-a(1)).*exp(1i*kz(2).*layers(2));1-a(1),(1+a(1)).*exp(1i*kz(2).*layers(2))];
M1TM=0.5*[1+b(1),(1-b(1)).*exp(1i*kz(2).*layers(2));1-b(1),(1+b(1)).*exp(1i*kz(2).*layers(2))];
MNTE=0.5*[(1+a(end)).*exp(-1i*kz(end-1).*layers(end-1));(1-a(end))];
MNTM=0.5*[(1+b(end)).*exp(-1i*kz(end-1).*layers(end-1));(1-b(end))];
% Transfer matrix is defined as layer 1-to-layer 2 transfer matrix. 
MTE=M1TE;
MTM=M1TM;
% Multiplying matrices in a recursive approach.
for n=2:size((Indices),2)-2
MLITE=0.5*[(1+a(n))*exp(-1i*kz(n)*layers(n)),(1-a(n))*exp(-1i*kz(n)*layers(n))*exp(1i*kz(n+1)*layers(n+1));1-a(n),(1+a(n))*exp(1i*kz(n+1).*layers(n+1))];
MTE=mtimes(MTE,MLITE);
MLITM=0.5*[(1+b(n))*exp(-1i*kz(n)*layers(n)),(1-b(n))*exp(-1i*kz(n)*layers(n))*exp(1i*kz(n+1)*layers(n+1));1-b(n),(1+b(n))*exp(1i*kz(n+1).*layers(n+1))];
MTM=mtimes(MTM,MLITM);

end
% Transfer matrix.
MTE=mtimes(MTE,MNTE);
MTM=mtimes(MTM,MNTM);
% Incident amplitude is set to 1
A1i=1;
% Reflected TE amplitude is evaluated in terms of the transfer matrix
A1rTE(q2,m)=MTE(2,1)/MTE(1,1);
A1rTM(q2,m)=MTM(2,1)/MTM(1,1);
end

end
%% Ploting graphs
% Assembling variables for a contour plot
X1=repmat(theta,length(P),1)*180/pi;
Y1=repmat(P,length(theta),1)'*100;
% Plotting TE reflectance map
figure(1)
hold on
surf(X1,Y1,abs(A1rTE).^2)
shading interp
set(gca,'FontSize',20)
xlim([0 85])
ylim([0 100])
title('Reflectance TE','FontSize',30,'interpreter','tex')
xlabel('Angle of incidence (deg)','FontSize',30,'interpreter','tex')
ylabel('Porosity (%)','FontSize',30,'interpreter','tex')
% Plotting TM reflectance map
figure(2)
hold on
surf(X1,Y1,abs(A1rTM).^2)
shading interp
set(gca,'FontSize',20)
xlim([0 85])
ylim([0 100])
title('Reflectance TM','FontSize',30,'interpreter','tex')
xlabel('Angle of incidence (deg)','FontSize',30,'interpreter','tex')
ylabel('Porosity (%)','FontSize',30,'interpreter','tex')
% Plotting the produc of TE and TM reflectance. This is made in order to
% avoid dividing by 0 out of the total internal reflection region.
figure(3)
hold on
surf(X1,Y1,0.5*(abs(A1rTE).^2+abs(A1rTM).^2))
shading interp
set(gca,'FontSize',20)
xlim([0 85])
ylim([0 100])
caxis([0 1])
title('Reflectance (TE vs. TM)','FontSize',30,'interpreter','tex')
xlabel('Angle of incidence (deg)','FontSize',30,'interpreter','tex')
ylabel('Porosity (%)','FontSize',30,'interpreter','tex')
%% Stop the time counter
toc
%% References
% [1] Ye, Q., Wang, J., Liu, Z., Deng, Z. C., Kong, X. T., Xing, F., ... 
% & Tian, J. G. (2013). Polarization-dependent optical absorption of 
% graphene under total internal reflection. Applied Physics Letters,
% 102(2), 021912.
%
% [2] Born, M., & Wolf, E. (2013). Principles of optics: electromagnetic
% theory of propagation, interference and diffraction of light. Elsevier.
end