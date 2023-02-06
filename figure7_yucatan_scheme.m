
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


% aw=dir('/Users/gaston/Desktop/Scripps/GoM_DT2018/GoM_DT2018/Eddies/adt*');
% 
% vort_merged=zeros(length(X),length(Y),length(dir));
% 
% 
% for i=1:length(aw)
% 
% date_of_interest=datenum('1993-1-1')+i-1;
% 
% datestr(date_of_interest,'yyyymmdd');
% load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
% 
% vort_merged(:,:,i)=Vorticity;
% 
% end
% 
% 
% 
% timee=datenum('1993-1-1'):datenum('2021-12-31')

load vort_merged

time=double(time)+datenum('1950-01-01 00:00:00');

 load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat', 'Shedding')

load mask_HYCOM_25.mat
load LC_VmaxContour_V2.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat
load vmerged
Xv=X;Zv=Z;
load pgon_yuc
load mask_HYCOM_25.mat
load LC_VmaxContour_V2.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies.mat
load LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat
 
load lc_rem



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

lag=30; %set the lag (in days)

smooth=24;


%in this way the moving average is actually the previous days

vdaily=vmerged-nmean(vmerged,3);

vdaily30=movmean(vdaily,[0 smooth],3);

vdaily30m=movmean(vmerged,[0 smooth],3);


cd /Users/gaston/Desktop/Scripps

input_dir='/Users/gaston/Desktop/Scripps/GoM_DT2018/GoM_DT2018/Eddies'


output_dir='/Users/gaston/Desktop/Scripps/video_loop_current_and_yucatan_flow_january_smooth_3/'

tdaily=time_merged;

time_step=tdaily(1):1:tdaily(end);

load Anclajes

% time_step=tdaily(601):1:tdaily(end);


%time_step=Shedding

%vdaily30=vdaily;

load lc_rem
lcnorm5=movmedian(lc_norm,1);

[C,ia,ib] = intersect(time_remi,time_merged);

lcnorm5=lcnorm5(ia);


ln5=movmedian(ln1,5);
% 
% for i=1:length(time_step)
% 
% date_of_interest=time_step(1)+i-1;

date_of_interest=datenum('2018-10-09');% retracted example

date_of_interest=datenum('2019-01-23');% retracted example
date_of_interest=datenum('2019-04-18');% retracted example


date_of_interest=datenum('2019-07-03');% cluster shedding

date_of_interest=datenum('2019-07-03');% cluster shedding

date_of_interest=datenum('2019-03-08');%the other cluster 
date_of_interest=datenum('2013-09-21');%the other cluster 



figure

set(gcf,'Position', [82         406        1012         207]);

subplot(1,3,1)

ind2=find(tdaily==date_of_interest);

pcolor(Xv,Zv,vdaily30(:,:,ind2)); hold on;%colorbar

cmocean('balance',16*2);caxis([-.8 0.8])
shading interp

contour(Xv,Zv,vdaily30(:,:,ind2),[-1.4:.1:1.4],'k','showtext','on','Labelspacing',360)

hold on
%axis ij
%colorbar; axis ij
%contour(X,Z,pposition,[0:5:50],'k','showtext','on','Labelspacing',360)

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -20]);
box on

 %title('Across channel velocity anomalies m.s^-^1')

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

hold on; yline(-20,'g','linewidth',2)




subplot(1,3,3)

ind2=find(tdaily==date_of_interest);

pcolor(Xv,Zv,vdaily30(:,:,ind2)); hold on;%colorbar

cmocean('balance',16*2);caxis([-.8 0.8])
shading interp

contour(Xv,Zv,vdaily30(:,:,ind2),[-1.4:.1:1.4],'k','showtext','on','Labelspacing',360)

hold on
%axis ij
%colorbar; axis ij
%contour(X,Z,pposition,[0:5:50],'k','showtext','on','Labelspacing',360)

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -20]);
box on

 %title('Across channel velocity anomalies (m.s^-^1)')

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

hold on; yline(-20,'g','linewidth',2)

subplot(1,3,1)

ind2=find(tdaily==date_of_interest);

pcolor(Xv,Zv,vdaily30m(:,:,ind2)); hold on;%colorbar

cmocean('balance',14*2);caxis([-1.4 1.4])
shading interp

contour(Xv,Zv,vdaily30m(:,:,ind2),[-1.4:.2:1.4],'k','showtext','on','Labelspacing',360)

hold on
%axis ij
%colorbar; axis ij
%contour(X,Z,pposition,[0:5:50],'k','showtext','on','Labelspacing',360)

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])


xlim([x_grid(1) x_grid(end)]);ylim([-2200 -20]);
box on

% title('Across channel velocity anomalies m.s^-^1')

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
ylabel('Depth (m)');

hold on; yline(-20,'g','linewidth',2)






subplot(1,3,2)



hold on

WL_plot=longitude(1)+360;EL_plot=longitude(end)+360;SL_plot=latitude(1);NL_plot=latitude(end);

load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
%X=X-360;
[lonx,latx]=meshgrid(X,Y);


m_proj('mercator','lon',[longitude(1)+360 longitude(end)+360], 'lat',[latitude(1) latitude(end)]);hold on

% [CS,CH]=m_etopo2('contourf',[-7000:1000:-1000 -200 0],'linestyle','none');

%[CS,CH]=m_etopo2('contour',[-7000:1000:-1000 -200 0],'color',[.7 .7 .7]);

%m_pcolor(X,Y,Vorticity');shading interp;cmocean('balance',21);caxis([-3e-5 3e-5])


speed=sqrt(U.*U+V.*V);

m_pcolor(X,Y,speed');shading interp;




m_etopo2('contour',[-200 -200],'color',[.7 .7 .7]);


[CS,CH]=m_etopo2('contour',[-200 -200],'color',[.7 .7 .7]);

%um(33,8)=1;
%vm(33,8)=0;

%hold on
%m_contour(longitude,latitude,adtm',[.17 .17],'r','linewidth',2)

%cmocean('balance',16)
%caxis([-.8 .8])
hold on
%m_plot(longitude(ix),latitude(iy),'m+')
%colorbar

% spacing=2;
% 
% %m_line([Lon(1) Lon(end)],[Lat(1) Lat(end)],'color','g','linewidth',2)
% m_quiver(lx(1:spacing:end),ly(1:spacing:end),um(1:spacing:end),vm(1:spacing:end),'k','autoscalefactor',1)
% 
% m_text(longitude(33),latitude(10),'1 m.s^-^1')

% title('Pearson correlation between Vgos and LC extension 1 month after (magenta >prctile99.5)')


%ax=m_contfbar([.05 .13],0.92,CS,CH,'endpiece','no','axfrac',.02);

%ax.XLabel.String='Depth (m)';



%[~,hc]=m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.6 .6 .6]);
m_contour(X,Y,ADT'*100,[-40:5:120],'color',[.2 .2 .2]);


% t=1;
% m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.4 .4 .4],'autoscalefactor',1.2);


% m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
% cmocean('delta');caxis([-1 1])


% m_pcolor(lon_sst,lat_sst,sst(:,:,i)');shading interp
% % colormap(flipud(brewermap(49,'Spectral')));
% colormap(jet);caxis([11 26])

%[ax,~]=m_contfbar([.18 .42],.865,[10 25],[10:1:25],'endpiece','no','axfrac',.05);
% cmocean('delta');caxis([-1 1])



% Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
% Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));

%  Cont_CEs=cat (3,Cyclonic_Cell(:,9), Cyclonic_Cell(:,10));
%  Cont_AEs=cat (3,Anticyclonic_Cell(:,9), Anticyclonic_Cell(:,10));
% 
% hold on
% 
% % 
% 
id_time=1;
date_num=date_of_interest;
% 
% id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
% id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
%     cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);
% 
% 
% id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;
% 
% id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
%     cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);
% 
% % All_eddy_AR=[AR_Traj{:,2}];
% 
% for loop_cyclo_id=id_in_anti'
%     m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'m','LineWidth',1)
% end
% 
% 
% 
% for loop_anti_id=id_in_cyclo'
%     m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'c','LineWidth',1)
% end
% 

indd=find(cell2mat(LC_VmaxContour(:,1))==date_of_interest);


 xx=LC_VmaxContour{indd,2};%xx(xx>180)=xx(xx>180)-360;

yy=LC_VmaxContour{indd,3};%ind=find(yy<21.25);yy(ind)=NaN;


m_plot(xx+360,yy,'Color','r','linewidth',2);
% m_plot(xx,yy,'Color',Colormap(id_min,:));

set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
set(findall(gcf,'-property','FontSize'),'FontSize',11);
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');
m_grid('linestyle','none');

hold on
m_line([Lon(1)+360 Lon(end)+360],[Lat(1) Lat(end)],'color','g','linewidth',2.5)

% ind=find(dates==dates(1)+i-1);

%  ind=find(dates==date_of_interest)

% m_plot(lon(ind),lat(ind),'*k')

%m_plot(lon,lat,'*k')
%title([ datestr(date_of_interest,'dd/mm/yy'),' LC ',num2str(round(lcnorm5(i))),' km '])

cmocean('-gray',12);caxis([0 1.2])

%cmocean('speed',12);caxis([0 1.2])

ax=m_contfbar([.455 .635],0.2,[0:.1:1.2],[0:.1:1.2],'endpiece','no','axfrac',.02);
ax=m_contfbar([.455 .635],0.2,[0:.1:1.2],[0:.1:1.2],'endpiece','no','axfrac',.02);

%ax.XLabel.String='Speed(ms^-^1)';

% print(gcf,'-painters','-depsc2','-r600','cluster4')

% print(gcf, '-dpng','-r200',[output_dir,datestr(time_step(i),'yyyymmdd')])
% 
% close all
% 
% end