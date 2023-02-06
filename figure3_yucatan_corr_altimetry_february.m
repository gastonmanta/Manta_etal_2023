%%%%%%load data%%%%%%%%%%%%%%

clear all; close all

ncdisp('sea_level_gulf_of_mexico_c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg.nc')

nc_filename = 'sea_level_gulf_of_mexico_c3s_obs-sl_glo_phy-ssh_my_twosat-l4-duacs-0.25deg.nc'; 

ncid=netcdf.open(nc_filename,'nowrite'); 

% Get information about the contents of the file.
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

disp(' '),disp(' '),disp(' ')
disp('________________________________________________________')
disp('^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~')
disp(['VARIABLES CONTAINED IN THE netCDF FILE: ' nc_filename ])
disp(' ')
for i = 0:numvars-1
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,i);
    disp(['--------------------< ' varname ' >---------------------'])
    flag = 0;
    for j = 0:numatts - 1
        attname1 = netcdf.inqAttName(ncid,i,j);
        attname2 = netcdf.getAtt(ncid,i,attname1);
        disp([attname1 ':  ' num2str(attname2)])
        if strmatch('add_offset',attname1)
            offset = attname2;
        end
        if strmatch('scale_factor',attname1)
            scale = attname2;
            flag = 1;
        end        
    end
    disp(' ')
    
    if flag
        eval([varname '= double(double(netcdf.getVar(ncid,i))*scale + offset);'])
    else
        eval([varname '= double(netcdf.getVar(ncid,i));'])   
    end
end
disp('^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~')
disp('________________________________________________________')
disp(' '),disp(' ')


time=double(time)+datenum('1950-01-01 00:00:00');


% 22 N between 87W and 84W


%[x y1]=find_close_value(latitude,22)

%[x x1]=find_close_value(longitude,-87)

%[x x2]=find_close_value(longitude,-84)


%lon_yuc=longitude(x1:x2);


%adt_yuc=squeeze(adt(x1:x2,y1,:));

%v_yuc=squeeze(vgos(x1:x2,y1,:));
%u_yuc=squeeze(ugos(x1:x2,y1,:));


%% howmollers

% 
% figure('Renderer', 'painters', 'Position', [50 50 800 600])
% pcolor(lon_yuc,time,adt_yuc');shading interp
% datetick('y','yyyy');axis tight
% title('ADT (m)  altimetry through Yucatan')
% colorbar;
% 
% 
% 
% figure('Renderer', 'painters', 'Position', [50 50 800 600])
% pcolor(lon_yuc,time,v_yuc');shading interp
% datetick('y','yyyy');axis tight
% title('Vgos altimetry through Yucatan')
% colorbar;
% 
% cmocean('balance');caxis([-1.2 1.2])
% set(gcf,'color','w');
% 
%  print(gcf,'-dpng','-r300','Howmoller_vgos_yuc')
% 
% 
% 
% 
% figure
% 
% plot(time,movmean(squeeze(nmean(v_yuc)),7))
% 
% title('Vgos mean altimetry through Yucatan')
% datetick('x','yyyy');axis tight
% set(gcf,'color','w');
% 
% 
% xlim([datenum('2013-07-01') datenum('2014-07-01')])
% 
% 
% print(gcf,'-dpng','-r300','Howmoller_vgos_yuc')
% 
% 
% twhole_direct =ncread(nc_file3,'OCEAN_VOLUME_TRANSPORT_WHOLE_SECTION_FROM_DIRECT_METHOD');
% 
% tlayers_geos=ncread(nc_file3,'OCEAN_VOLUME_TRANSPORT_IN_LAYERS_FROM_GEOSTROPHIC_METHOD');
% 
% tlayers_direct=ncread(nc_file3,'OCEAN_VOLUME_TRANSPORT_IN_LAYERS_FROM_DIRECT_METHOD');
% 
% LAYER_CENTER_LOCATION_VERTICAL=ncread(nc_file3,'LAYER_CENTER_LOCATION_VERTICAL');


%% I will grid the data into the YUC1 to YUC9

%load the grid and bathymetry polygon of the mexican section
load  pgon_yuc

lon_yucs=x_grid;
lat_yucs=y_grid;

for i =1:length(adt)

adt_yucs(:,i)=griddata(longitude,latitude,adt(:,:,i)',lon_yucs,lat_yucs);

v_yucs(:,i)=griddata(longitude,latitude,vgos(:,:,i)',lon_yucs,lat_yucs);

u_yucs(:,i)=griddata(longitude,latitude,ugos(:,:,i)',lon_yucs,lat_yucs);

end


%make v normal to the section
i=-10;

un= (u_yucs*cosd(i))+(v_yucs*sind(i));
vn= -(u_yucs*sind(i))+(v_yucs*cosd(i));



load('vmerged')
vdaily=vmerged;
tdaily=time_merged;lon=X;z=Z;

% correlation in situ with itself

for i=1:length(vdaily(:,1,1))

for j=1:length(vdaily(1,:,1))

i1=squeeze(vdaily(i,j,:))';

i2=squeeze(vdaily(1,j,:));

[rhos(i,j),pvals(i,j)]=corr(i1',i2,'rows','pairwise');

%rmsees(i,j)=rmse(i1',i2');

end

end

rhoos=rhos;
ind=find(pvals>0.01); rhoos(ind)=NaN;


% correlation of altimerty with in situ data at different smoothings

beg=find(tdaily(1)==time);

endd=find(tdaily(end)==time);

for i=1:length(vdaily(:,1,1))

for j=1:length(vdaily(1,:,1))

i1=squeeze(vdaily(i,j,:))';

i2=squeeze(vn(j,beg:endd));

[rho(i,j),pval(i,j)]=corr(i1',i2','rows','pairwise');

%rmsee(i,j)=rmse(i1',i2');

end

end

rhoo=rho;
ind=find(pval>0.01); rhoo(ind)=NaN;

% 
% for i=1:length(vdaily(:,1,1))
% 
% for j=1:length(vdaily(1,:,1))
% 
% i1=squeeze(vdaily(i,j,:))';
% 
% i2=squeeze(vn(j,beg:endd));
% 
% [rho30(i,j),pval30(i,j)]=corr(i1',i2','rows','pairwise');
% 
% %rmsee(i,j)=rmse(i1',i2');
% 
% end
% 
% end
% 
% rhoo30=rho30;
% ind30=find(pval30>0.001); rhoo(ind)=NaN;
% 




% movmean as centered moving average

k=50;
vdailyk=zeros(length(vdaily(:,1,1)),length(vdaily(1,:,1)),length(vdaily),k)*NaN;

for k=1:50
vdailyk(:,:,:,k)=movmean(vdaily,k,3);
end


for k=1:50

for i=1:length(vdaily(:,1,1))

for j=1:length(vdaily(1,:,1))


i1=squeeze(vdailyk(i,j,:,k))';

i2=squeeze(vn(j,beg:endd));

[rho30(i,j,k),pval30(i,j,k)]=corr(i1',i2','rows','pairwise');

% vva=squeeze(vdailyanom(i,j,1:length(tdaily)-k+1));
% lag_corra(i,j,k)=corr(vva,ln5cut(k:length(tdaily))');

end
end
end



rhoo30=rho30;
ind=find(pval30>0.01); rhoo30(ind)=NaN;


amax=squeeze(nmax(nmax(rho30)))

amean=squeeze(nmean(nmean(rho30)))

akmax=squeeze(nmax(nmax(rho30)))

akmean=squeeze(nmean(nmean(rho30)))


% [qq,qw]=max(max(rho30(:,:,21)));
% 
% [qqa,qwa]=max(max(rho30(:,:,23)));
% 
% [ii,jj]=find(rho30(:,:,21)==max(akmax(21)));
% 
% [iia,jja]=find(rho30(:,:,23)==max(akmax(22)));



% just lon lat positions of the moorings 
load Anclajes




figure('Renderer', 'painters', 'Position', [0 0 800 600])
subplot(2,2,1)
hold on

contourf(lon,z,rhoos,[0:.05:1],'linestyle','none');%shading interp
contour(lon,z,rhoos,[0:.1:1],'k','showtext','on','Labelspacing',440);%shading interp

for i=1:length(Lon)
xline(Lon(i),'--','color',[0.81818 0.77647 0.70909])
%plot(Lon(i),-20,'dk')
end

%plot(pgon_yuc,'FaceColor',[.7 .7 .7])

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])

xlim([lon(1) lon(end)]);ylim([-2200 0]);
box on
%colormap(rednblue)
%cmocean('balance');
 %caxis([-120 120])
%title('RHO geostrophic (only p<0.001)')

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'});
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
%ylabel('Depth (m)');

caxis([0 1])

box on
ylabel('Depth (m)');

subplot(2,2,2)
hold on
contourf(lon,z,rho30(:,:,1),[0:.05:1],'linestyle','none');%shading interp
contour(lon,z,rho30(:,:,1),[0:.1:1],'k','showtext','on','Labelspacing',440);%shading interp

for i=1:length(Lon)
xline(Lon(i),'--','color',[0.81818 0.77647 0.70909])
%plot(Lon(i),-20,'dk')
end

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909]);
xlim([lon(1) lon(end)]);ylim([-2200 0]);
box on
%colormap(rednblue)
%cmocean('balance');
 %caxis([-120 120])


xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'});
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

set(gcf,'color','w');;set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on


subplot(2,2,3)
hold on
contourf(lon,z,rho30(:,:,21),[0:.05:1],'linestyle','none');%shading interp
contour(lon,z,rho30(:,:,21),[0:.1:1],'k','showtext','on','Labelspacing',440);%shading interp

for i=1:length(Lon)
xline(Lon(i),'--','color',[0.81818 0.77647 0.70909])
%plot(Lon(i),-20,'dk')
end

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909]);
xlim([lon(1) lon(end)]);ylim([-2200 0]);
box on
%colormap(rednblue)
%cmocean('balance');
 %caxis([-120 120])


xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'});
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on




[ss,aa]=max(akmax)

subplot(2,2,4)
plot(akmax,'k')

hold on

plot(aa,akmax(aa),'+r')

ylabel(' Max. Correlation')

xlabel('Moving average (days)')

grid on


set(findall(gcf,'-property','FontSize'),'FontSize',12);


% savefig('figure3_corr_altimetry_february')



% savefig('figure3_corr_altimetry_february')
% print(gcf,'-painters','-depsc2','-r600','figure3_corr_altimetry_february')
%  print(gcf,'-dpng','-r600','figure3_corr_altimetry_february')


rho21=rho30(:,:,21)



load rhoe30 




ww=abs(rho30)-abs(rho21);

figure


hold on
contourf(lon,z,ww,[-1:.1:1],'linestyle','none');%shading interp
contour(lon,z,ww,[-1:.1:1],'k','showtext','on','Labelspacing',440);%shading interp

for i=1:length(Lon)
xline(Lon(i),'--','color',[0.81818 0.77647 0.70909])
%plot(Lon(i),-20,'dk')
end

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909]);
xlim([lon(1) lon(end)]);ylim([-2200 0]);
box on
%colormap(rednblue)
cmocean('balance');caxis([-.8 .8])
 %caxis([-120 120])


xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'});
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

set(gcf,'color','w');;set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on


% aa=squeeze(vmerged(73,57,:))
% 

% figure;
% 
% plot(nmean(v_yucs,2))
% 
% 
% [C,ia,ib] = intersect(time,time_merged);
% 
% vnn=vn(:,ia);
% 
% vyucanom=vnn-nmean(vnn,2);
% 
% 
% %Si tenes una matriz ff de dim (lt,ly,lx) haces:
% ff=(permute(vyucanom,[2 1]));
% 
% ff=ff(:,:);   %pasa a tener dim  (lt,ly*lx)
% mapf=ff(1,:);    % Save geographical dimension to use in mask.
% ff(:,any(isnan(ff)))=[];%Remove NaNs
% [lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)
% 
% 
% 
% %V33 = clusterdata(ff,'linkage','ward','savememory','on','maxclust',3);
% 
% Vf = clusterdata(ff,'linkage','ward','savememory','on','maxclust',4);
% 
% 
% 
% 
% 
% figure;hold on
% ind=find(Vf==1);
% 
% 
% plot(nmean(vyucanom(:,ind),2))
% 
% ind2=find(Vf==2);
% plot(nmean(vyucanom(:,ind2),2))
% 
% ind3=find(Vf==3);
% plot(nmean(vyucanom(:,ind3),2))
% 
% ind4=find(Vf==4);
% plot(nmean(vyucanom(:,ind4),2))
% 
% 
% 
% figure;hold on
% ind=find(Vf==1);
% 
% 
% stem(nmean(vyucanom(:,ind),2))
% 
% ind2=find(Vf==2);
% stem(nmean(vyucanom(:,ind2),2))
% 
% ind3=find(Vf==3);
% stem(nmean(vyucanom(:,ind3),2))
% 
% ind4=find(Vf==4);
% stem(nmean(vyucanom(:,ind4),2))
% 
% 
% 
% load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat')
% 
% 
% load lc_rem.mat
% 
% lcnorm5=movmedian(lc_norm,5);
% lcnorm5=lcnorm5(ia);
% 
% 
% 
% figure
% % plot(tdaily,ln5,'k'); hold on
% %ylim([500 3000]);yyaxis right
% %plot(tdaily,lcnorm5,'k'); hold on
% scatter(tdaily,lcnorm5,20,Vf,'filled');cmocean('thermal',4);colorbar
% 
% ylim([500 3000])
% hold on
% 
% for i=1:length(Shedding)
% xline(Shedding(i),'r','linewidth',.65)
% end
% 
% 
% for i=1:length(Detached)
% xline(Detached(i),'color', [0 .4 0],'linewidth',.65)
% end
% 
% 
% for i=1:length(Reattached)
% xline(Reattached(i),'color', [1.0000    0.8398 0],'linewidth',.65)
% end
% 
% 
% datetick('x')
% %xlim([datenum('1993-1-1') datenum('2022-1-1')])
% xlim([time_merged(1) time_merged(end)])
% 
% xticklabels([])
% ylabel('LC extension (km)')
% box on; grid on
% 
% 
% 
% legend('0.17m','< V >','Eddy Shedding','Detached','Reattached','orientation','horizontal','box','off')
% 
%  set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
% 
% 
% ylabel('Depth (m)');
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12);
% 
% 
% 
% % save Vf Vf
% % 
% 
% 
%  [C,ia,ib] = intersect(Shedding,time_merged);
% 
% 
% 
% figure;hist(Vf(ib),4)
% 
% 
% 
% 
% 
% 



