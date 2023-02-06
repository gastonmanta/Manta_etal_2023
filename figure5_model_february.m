%% 

%STILL NEED TO CHECK THE DIFFERENT GRIDS
%THE BEST LONGITUDE TO PUT THE MOORING
%

%% THE MOST IMPORTANT PART. PICK THE LAG FOR CORRELATION AND SMOOTHING OF THE DATA
lag=30; %set the lag (in days)

smooth=21;

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





um=nmean(ugos,3);

vm=nmean(vgos,3);


adtm=nmean(adt,3);
adtt=nmean(nmean(adtm));
adtm=adtm-adtt;

spacing=3;


load('Anclajes.mat')

% 


ln5=movmedian(ln1,5);

load lc_rem
lcnorm5=movmedian(lc_norm,1);
lcnorm5=lc_norm;


%load YC observations
load('vmerged.mat');tdaily=time_merged;vdaily=vmerged;

%time intersection
[C,ia,ib] = intersect(time_remi,tdaily);

lcnormcut=lcnorm5(ia);


%in this way the moving average is actually the previous days
%vgos30=movmean(vgos,[0 smooth],3);

vgos30=vgos;

[C,ia,ib] = intersect(time,tdaily);

vgos30=vgos30(:,:,ia);

%ln30=movmean(lcnormcut,[0 smooth]);
ln30=lcnormcut;

%nr30=movmean(nr5,[0 lag]);

% searchs for regions of high correlation of Vgos from altimetry with llag


%ln30=ln30';


for i=1:length(longitude)
for j=1:length(latitude)

% vv=squeeze(adt(i,j,:));
%vv=squeeze(ugos(i,j,:));
%vv=squeeze(vgos(i,j,:));
%[rho(i,j),pval(i,j)]=corr(vv,ln5','rows','pairwise');
%[rho30(i,j),pval30(i,j)]=corr(vv(1:end-lag+1),ln5(lag:end)','rows','pairwise');

vv=squeeze(vgos30(i,j,:));
[rhoalt(i,j),pvalalt(i,j)]=corr(vv,ln30,'rows','pairwise');
[rho30alt(i,j),pval30alt(i,j)]=corr(vv(1:end-lag+1),ln30(lag:end),'rows','pairwise');

end
end


p05=prctile(rho30alt(:),.5)

p995=prctile(rho30alt(:),99.5)

ind=find(rho30alt>p995);

ind2=find(rho30alt<p05);

[ixa,iya]=find(rho30alt>p995);

[ix1,iy1]=find(rho30alt<p05);

v995=vgos(ixa,iya,:);v995=squeeze(nmean(nmean(v995)));

v5=vgos(ix1,iy1,:);v5=squeeze(nmean(nmean(v5)));


% %time intersection
% [C,ia,ib] = intersect(time,tdaily);
% 
% ln5cut=ln5(ia);

v995cut=v995;%(ia);

v5cut=v5;%(ia);

%v995=vgos(ix,iy,:);v995=squeeze(nmean(nmean(v995)));
%v05=vgos(ix1,iy1,:);v05=squeeze(nmean(nmean(v05)));

% [lag_corr,max_lag_corr,lag_,smooth_,lag_rmse,lag_pval,lagg_rmse,smooth_rmse]=lagged_smooth_correlation(v995(1:end-lag+1),ln30(lag:end),2,1920)



v995cut=v995(ia);
v5cut=v5(ia);


% %-86.0800
% predictorrslon1=squeeze(vdaily(:,25,:));
% % -85.3300
% predictorrslon2=squeeze(vdaily(:,50,:));
% 
% predictors=[predictorrslon1 ;predictorrslon2];
% vv=permute(vdaily,[3 1 2]);vv=vv(:,:);
% 
% v1=vv(:,maxrrho30);
% v2=vv(:,minrrho30);
% v3=vv(:,maxrrho30sub);



% 
% lag=21;
% 
% ln30=movmean(ln5cut,[0 lag]);
% 
% corr(v995cut(1:end-lag+1),ln30(lag:end)','rows','pairwise')
% 


% smooth=90;
% predictors=movmean(predictors,smooth,1);
% predictand=movmean(ln5cut_lag,smooth);


vdaily60=movmean(vdaily,[0 smooth],3);

% 
% load lc_rem
% lcnorm5=movmedian(lc_norm,5);

%%%%%  TO REPLACE LN5 FROM GANESH TO LCNORM5 FROM REMI %%%%% 
%ln5=lcnorm5;
%[C,ia,ib] = intersect(time_remi,tdaily);
%ln5cut=ln5(ia);

%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%lnnorm5cut60=movmean(lnnorm5,[0 ssmooth]);
lnnorm5cut60=movmean(lcnorm5,[0 smooth]);

lnnorm5cut60=lcnorm5;
%ln5cut30=movmean(ln5cut,[0 360]);



[C,ia,ib] = intersect(time_remi,tdaily);

lnnorm5cut60=lnnorm5cut60(ia);


lnnorm5cut60=lnnorm5cut60';

for j=1:length(X(1,:))
for i=1:length(Z(:,1))

vv=squeeze(vdaily60(i,j,:));
%vv=squeeze(vdailyanom_normalize(i,j,:));

[rho(i,j),pval(i,j)]=corr(vv,lnnorm5cut60','rows','pairwise');

[rho30(i,j),pval30(i,j)]=corr(vv(1:end-lag+1),lnnorm5cut60(lag:end)','rows','pairwise');
%[rho30(i,j),pval30(i,j)]=corr(vv(1:end-llag+1),ln5cut60(llag:end)','rows','pairwise');

end
end


load pgon_yuc

rrho30=rho30;

ind=find(pval30>0.001);rrho30(ind)=NaN;


rrho30s=rrho30;





% 
% ind8611=find(x_grid==-86.11)
% 
% zx2=find(Z(maxrrho30mid)==Z(:,1))
% 
% zx3=find(Z(minrrho301200)==Z(:,1))
% zx4=find(Z(maxrrho30deep)==Z(:,1))
% 
% x2=squeeze(vdaily60(zx2,ind8611,:));
% 
% 
% x2=squeeze(vdaily60(ind8611,zx2,:));
% 
% x3=squeeze(vdaily60(zx3,ind8611,:));
% 
% x4=squeeze(vdaily60(zx4,ind8611,:));


% predictors_mooring1=[vmax vmid  v1200min vdeep]
% 
% x1=squeeze(vdaily60(iy,iz,:))
% 
% predictors_mooring2=[vmin vmidmin v1200max vdeepmin]


%predictors=[v995cut vmax vmid vdeep vdeepmin];

%predictors=[v995cut predictors_mooring1 predictors_mooring2]

% model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% 
% predictors=[v995cut predictors_mooring1]


%v05cut=v05(ia);


%predictors=[v995cut vmid  v1130min vdeep]


%predictors=[v995cut vmid vdeep]

%predictors=[v995cut vmid v1130min vdeep]

% vdaily60=vdaily;

% 76 1510m 55 1090 38 750

%86 w es 28

%[28 44 EN Z 76 55 Y 38 EXPLICAN 87%)

tt1=26; %87


tt=44; %87

JJ=0;
 %87


predictors=[v995cut v5cut squeeze(vdaily60(76+JJ,tt1,:)) squeeze(vdaily60(76+JJ,tt,:))];

predictors=[v995cut v5cut squeeze(vdaily60(76+JJ,tt1,:)) squeeze(vdaily60(56+JJ,36,:))];

predictors=[v995cut v5cut squeeze(vdaily60(76+JJ,tt1,:)) squeeze(vdaily60(56+JJ,34,:))];


predictors=[v995cut v5cut squeeze(vdaily60(77+JJ,tt1,:)) squeeze(vdaily60(59+JJ,tt1,:)) squeeze(vdaily60(37+JJ,tt1,:))];


predictors=[v995cut v5cut squeeze(vdaily60(77+JJ,tt1,:)) squeeze(vdaily60(61+JJ,tt1,:))];


%predictors=[v995cut v5cut];

predictand=lnnorm5cut60;

predictand=lcnormcut;

 %predictors=[vmax vmin vmid vmidmin v1130min vdeep vdeepmin]

% model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))

% predictors=[v995cut v5cut]

train_period=round((length(predictors)/3)*2);

model1=fitlm(predictors(1:end-lag+1-train_period,:),predictand(lag:end-train_period))

ypred = predict(model1,predictors);



rmse(ypred,lnnorm5cut60')


corr(ypred,lnnorm5cut60')


rmse(ypred(1:length(ypred)/3*2),lnnorm5cut60(1:length(ypred)/3*2)')


rmse(ypred(1:length(ypred)/3*2),lnnorm5cut60(1:length(ypred)/3*2)')


rmse(ypred(1:length(ypred)/3*2),lcnormcut(1:length(ypred)/3*2))



rmse(ypred,lnnorm5cut60')


rmse(ypred,lnnorm5cut60')



figure('Renderer', 'painters', 'Position', [0 0 1100 700])
subplot(2,2,1)
m_proj('mercator','lon',[longitude(1) longitude(end)], 'lat',[latitude(1) latitude(end)]);hold on
m_pcolor(longitude,latitude,rho30alt');shading interp;

m_contour(longitude,latitude,rho30alt',[-.4 -.6 .4 .6],'k','showtext','on','Labelspacing',360)

% m_contour(longitude,latitude,rho30alt',[.6 .6],'k','showtext','on','Labelspacing',360)

m_contour(longitude,latitude,rho30alt',[.625 .625 ],'c','linewidth',1)

m_contour(longitude,latitude,rho30alt',[-.61 -.61],'c','linewidth',1)

%m_contour(longitude,latitude,rho30alt',[-.8:.1:.8 ],'k','showtext','on','Labelspacing',360)

cmocean('balance',16)
caxis([-.8 .8])
hold on
 %m_plot(longitude(ixa),latitude(iya),'m+','linewidth',1)

% m_plot(longitude(ix1),latitude(iy1),'m+','linewidth',1)

%colorbar
m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
m_line([Lon(1) Lon(end)],[Lat(1) Lat(end)],'color','g','linewidth',2)
set(gcf,'color','w');

m_grid('linestyle','none');

%title('Pearson correlation between Vgos and LC extension 1 month after (magenta >prctile99.5)')
set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
set(findall(gcf,'-property','FontSize'),'FontSize',12);

% ax=m_contfbar([.46 .62],0.22,[-.8:.1:.8],[-.8:.1:.8],'endpiece','no','axfrac',.02);
% ax.XLabel.String='Correlation';

subplot(2,2,2)
pcolor(X,Z,rrho30); hold on;

line([Lon(1) Lon(end)],[0 0],'color','g','linewidth',2)

cmocean('balance','pivot');%caxis([-.6 .6])
shading interp
%axis ij
%colorbar; axis ij
% for i=1:length(Lon)
% xline(Lon(i),'--','color',[0.81818 0.77647 0.70909])
% %plot(Lon(i),-20,'dk')
% end

contour(X,Z,rrho30,[-1:.1:1],'k','showtext','on','Labelspacing',360)

plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 0]);
box on

%title(' Correlation between observed Vdaily and LC REMI extension (only p<0.001)')
xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
%xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');

set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on


axx=m_contfbar([.2 .4],0.21,[-.8:.1:.8],[-.8:.1:.8],'endpiece','no','axfrac',.02);
axx.XLabel.String='Correlation';
% % 
% % plot(X(minrrho30),Z(minrrho30),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(maxrrho30),Z(maxrrho30),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(maxrrho30mid),Z(maxrrho30mid),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(maxrrho30deep),Z(maxrrho30deep),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(minrrho30mid),Z(minrrho30mid),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(minrrho30deep),Z(minrrho30deep),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(minrrho301200),Z(minrrho301200),'+m','markersize',10,'linewidth',2)
% % 
% % plot(X(maxrrho301200),Z(maxrrho301200),'+m','markersize',10,'linewidth',2)


%line([-86.11 -86.11],[-1650 -700],'color',[.7 .7 .7],'linewidth',2)

line([x_grid(tt1) x_grid(tt1)],[-1650 Z(59,1)],'color',[.7 .7 .7],'linewidth',2)

plot(x_grid(tt1) ,Z(77,1),'oc','markersize',10,'linewidth',2)

plot(x_grid(tt1) ,Z(61,1),'oc','markersize',10,'linewidth',2)

%plot(-86.11 ,Z(38,1),'oc','markersize',10,'linewidth',2)

caxis([-.8 .8])



% xline(x_grid(28),'--','color','c','linewidth',2)
% xline(x_grid(44),'--','color','c','linewidth',2)

subplot(2,2,3:4)
plot(time_merged,lnnorm5cut60,'color',[.6 .6 .6])
plot(time_merged,lcnormcut,'color',[.6 .6 .6])


hold on
plot(time_merged(round(1:round(length(time_merged)/3*2))),ypred(1:round(length(time_merged)/3*2)),'r')
plot(time_merged(round(round(length(time_merged)/3*2)):end),ypred(round(length(time_merged)/3*2):end),'--r')


title([' R^2=',num2str(round((model1.Rsquared.Ordinary),2)),'. RMSE=',num2str(round(rmse(lnnorm5cut60',ypred))),'km '])

legend('Observed','Predicted','orientation','horizontal','autoupdate','off','box','off')

datetick('x')
ylabel('Loop Current extension (km)')

axis tight
grid on


set(findall(gcf,'-property','FontSize'),'FontSize',12);
set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on


%xline(time_merged(round(length(time_merged)/3*2)),'y')


print(gcf,'-painters','-depsc2','-r600','figure5_model_february')
 
% savefig('figure5_model')








aa={};

for i=1:length(Detached)

aa{i}=Detached(i):Reattached(i)


end

detret=cell2mat(aa);


length(detret)/length(lc_norm)


[C,ia,ib] = intersect(detret,tdaily)



lcnormcutf=lcnormcut;lcnormcutf(ib)=0;ind=find(lcnormcutf==0);lcnormcutf(ind)=[];
ypredf=ypred;ypredf(ib)=0;ind=find(ypredf==0);ypredf(ind)=[];



rmse(lcnormcut,ypred)


rmse(lcnormcutf,ypredf)


% lcnorm5_sin=lcnorm5cut;lcnorm5_sin(ib)=NaN;
% 
% ypredf=ypred;ypredf(ib)=NaN;
% 
% 
% length(ib)/length(lcnorm5cut)
% 
% lcnorm5_sin_f=fillmissing(lcnorm5_sin,'linear')
% 
% 











