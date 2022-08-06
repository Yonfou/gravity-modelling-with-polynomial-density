clc; clear all; close all

% equally-spaced sampling points in the longitudinal direction for the source region
lon = [-179.5:179.5].'/180*pi;  nsl = length(lon);
dlmwrite('source_longitude.dat',lon,'delimiter','\t','newline','pc','precision','%30.18f')

% latitutes of the center of each element (SF) and the interval (dSF) in the latitudinal direction for the source region. 
% Each rwo has 2 elements and there are NSF rows.
lat = [-89.5:89.5].'/180*pi;  nsf = length(lat);
dlat = ones(nsf,1)/180*pi;
dlmwrite('source_latitude.dat',[lat,dlat],'delimiter','\t','newline','pc','precision','%30.18f')

% the total number of sampling points in the radial direction at different latitudes for the source region
nr = 10*ones(nsf,1);  nsr = max(nr);
dlmwrite('source_number_in_depth.dat',nr,'delimiter','\t','newline','pc')

% radial coordinates of the center of each element (SR) and the interval (dSR) in the radial direction 
% at different latitudes for the source region. The outer loop is r. Note that the sampling points are the same along the longitudinal direction. 
% Each row has 2 elements and there are NSR(=MAX(NR))¡ÁNSF rows
sr0 = 6271e3;  sr1 = 6371e3;
dsr_= (sr1-sr0)/nr(1);
sr_= (sr0+0.5*dsr_):dsr_:(sr1-0.5*dsr_);
sr = zeros(nsf,nsr);
dsr = dsr_*ones(nsf,nsr);
for ii=1:nsf
    sr(ii,:)=sr_(:);
end
dlmwrite('source_depth.dat',[sr(:),dsr(:)],'delimiter','\t','newline','pc','precision','%30.18f')

% the degree of polynomial function used in each element along the radial direction. 
% Each row has NSF elements and there are NSR rows.
np=zeros(nsr,nsf);
dlmwrite('polynomial_degree.dat',np,'delimiter','\t','newline','pc')


% the location (r,¦Õ,¦Ë) of the observation points, the inner loop is ¦Ë and the outermost loop is r. 
% Note that the interval in the longitudinal direction is the same as those in the source region. 
% Each rwo has 3 elements and there are NPL*NPF*NPR rows.
npl=nsl;  npf=nsf;  npr=1;
dat=zeros(npl*npf*npr,3);

i0=0;
for ii=1:npr
for jj=1:npf
for kk=1:npl
    i0=i0+1;
    dat(i0,1)=sr1+10.0e3;
    dat(i0,2)=lat(jj);
    dat(i0,3)=lon(kk);
end
end
end

dlmwrite('observation.dat',dat,'delimiter','\t','newline','pc','precision','%30.18f')

% File that contains the dentisity of the model. 
% Each row has 1 elements and there are NSR*NSL*NSF (only for homogeneous density) or (NSR*MAXVAL(NP)+1)*NSL*NSF rows. 
% The inner loop is ¦Ë and the outermost loop is r.
rho=1e3*ones(nsl,nsf,nsr,max(max(np))+1);
rho=rho(:);
dlmwrite('density_model.dat',rho,'delimiter','\t','newline','pc','precision','%30.18f')