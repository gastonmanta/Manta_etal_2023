
%load('Yuc_Vel_2012_2018_gridded.mat')

load('vmerged.mat')


load('Yuc_eofs_f7_n800b.mat','bp')


load pgon_yuc
% v=zeros(length(X(:,1)),length(Z(1,:)),length(Velh(1,:))).*NaN;
% 
% for i=1:length(Velh(1,:))
% 
% vv=zeros(length(X(:,1)),length(Z(1,:))).*NaN;vv(celda)=Velh(:,i);
% 
% v(:,:,i)=vv;
% 
% end
% 
% lon=X(1,:);z=Z(:,1);
% 
% j=1;
% for i=1:24:length(v)
%     vdaily(:,:,j)=nmean(v(:,:,i:i+11),3);
%     j=j+1;
% end
% 




mean_julio=nmean(vmerged,3);

std_julio=nstd(vmerged,1,3);

lon=X(1,:);

z=Z(:,1);

load('Anclajes.mat')








%%%%%%load data%%%%%%%%%%%%%%

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

 load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat', 'Shedding')

load mask_HYCOM_25.mat
load LC_VmaxContour_V2.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat



% 22 N between 87W and 84W

[x y1]=find_close_value(latitude,22);

[x x1]=find_close_value(longitude,-87);

[x x2]=find_close_value(longitude,-84);

lon_yuc=longitude(x1:x2);

adt_yuc=squeeze(adt(x1:x2,y1,:));
v_yuc=squeeze(vgos(x1:x2,y1,:));
u_yuc=squeeze(ugos(x1:x2,y1,:));



% mask for GOM

adtcut=adt(:,:,1);

[lx,ly]=meshgrid(longitude,latitude);ly=ly';lx=lx';

ind=find(ly<21.2 & lx>-89);
adtcut(ind)=NaN;

ind=find(ly<22.6 & lx>-84.6);
adtcut(ind)=NaN;

ind=find(ly<23 & lx>-84);
adtcut(ind)=NaN;

maskadt=adtcut./adtcut;


% 
% figure;pcolor(lx,ly,maskadt);shading interp

adt=adt.*maskadt;

%level for LC extension
lev = 0.17;


for i=1:length(time)

adtt=nmean(nmean(adt(:,:,i)));

C = contourc(longitude,latitude,(adt(:,:,i)-adtt)',[lev lev]);

[ws nr ln ar] = ssh17cmtck_rev(C);

ws=ws(1);nr=nr(1);ln=ln(1);ar=ar(1);

ws1(i)=ws;
nr1(i)=nr;
ln1(i)=ln;
ar1(i)=ar;


end

% ln1 is the LC extension. ln5 is just to remove spikes
ln5=movmedian(ln1,5);
nr5=movmedian(nr1,5);



load lc_rem
lcnorm5=movmedian(lc_norm,5);



load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat', 'Shedding')





um=nmean(ugos,3);

vm=nmean(vgos,3);


adtm=nmean(adt,3);
adtt=nmean(nmean(adtm));
adtm=adtm-adtt;

spacing=3;


load('Anclajes.mat')



figure('Renderer', 'painters', 'Position', [0 0 1100 350])
subplot(1,2,1)
m_proj('mercator','lon',[longitude(1) longitude(end)], 'lat',[latitude(1) latitude(end)]);hold on

[CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -200 0],'linestyle','none');

[CS,CH]=m_etopo2('contour',[-7000:1000:-1000 -200 0],'color',[.7 .7 .7]);

%[CS,CH]=m_etopo2('contourf',[-6000:1000:-1000 -200],'showtext','on','labelspacing',150);


colormap([m_colmap('blues',128)]);



um(33,8)=1;
vm(33,8)=0;

hold on
m_contour(longitude,latitude,adtm',[.17 .17],'r','linewidth',2)

%cmocean('balance',16)
%caxis([-.8 .8])
hold on
%m_plot(longitude(ix),latitude(iy),'m+')
%colorbar
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');

m_line([Lon(1) Lon(end)],[Lat(1) Lat(end)],'color','g','linewidth',2)
m_quiver(lx(1:spacing:end),ly(1:spacing:end),um(1:spacing:end),vm(1:spacing:end),'k','autoscalefactor',1)
m_text(longitude(33),latitude(10),'1 m.s^-^1')

% title('Pearson correlation between Vgos and LC extension 1 month after (magenta >prctile99.5)')
set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12);
m_grid('linestyle','none');

%ax=m_contfbar([.05 .13],0.92,CS,CH,'endpiece','no','axfrac',.02);

%ax.XLabel.String='Depth (m)';

set(gcf,'color','w');





figure('Renderer', 'painters', 'Position', [0 0 1100 350])
subplot(1,2,1)
hold on
contourf(lon,z,mean_julio,[-.20:.05:1.20],'linestyle','none');%shading interp

[C,hContour] = contour(lon,z,mean_julio,[-.20:.10:.20 .30 .70 1.10],'k','showtext','on','Labelspacing',440);%shading interp
clabel(C,hContour,'FontSize',12,'Color','k');

%contour(lon,z,mean_julio,[-25:10:125],'k');%shading interp

% for i=1:length(Lon)
% %xline(Lon(i),'--','color',[.8 .8 .8])
% plot(Lon(i),-20,'dk')
% end

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909]);

xlim([lon(1) lon(end)]);ylim([-2200 -20]);
box on
colormap(rednblue)
%cmocean('balance');
 caxis([-1 1])

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');
box on
%title('Mean velocity (m.s^-^1)')

subplot(1,2,2)
hold on
contourf(lon,z,std_julio,[0:.05:.50],'linestyle','none');%shading interp

[C,hContour] =contour(lon,z,std_julio,[0:.10:.50],'k','showtext','on','Labelspacing',1440);%shading interp
clabel(C,hContour,'FontSize',12,'Color','k');

%  for i=1:length(Lon)
% %xline(Lon(i),'--','color',[.8 .8 .8])
% plot(Lon(i),-20,'dk')
% end

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909]);
xlim([lon(1) lon(end)]);ylim([-2200 -20]);
box on

colormap(rednblue)

caxis([-.40 .40])


xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');


%title('Standard deviation (m.s^-^1)')

set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12);
box on; set(gcf,'color','w');

% % 
% % print(gcf,'-painters','-depsc2','-r600','yuc_mean_and_std')
% % print(gcf, '-dpng','-r600','yuc_mean_and_std')



% print(gcf,'-painters','-depsc2','-r600','fig1a_yuc')





