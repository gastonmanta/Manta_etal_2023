%% TOEddies/LC and LCE for el Maestro Gaston

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

load vmerged; vdaily=vmerged;tdaily=time_merged;

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


lag=30; %set the lag (in days)
llag=30;

smooth=21;
ssmooth=21;


%in this way the moving average is actually the previous days
vgos30=movmean(vgos,[0 smooth],3);

%ln30=movmean(ln5,[0 smooth]);
ln30=ln5;
%nr30=movmean(nr5,[0 lag]);


vdaily60=movmean(vdaily,[0 ssmooth],3);


load lc_rem
lcnorm5=movmedian(lc_norm,5);

%%%%%  TO REPLACE LN5 FROM GANESH TO LCNORM5 FROM REMI %%%%% 
%ln5=lcnorm5;
%[C,ia,ib] = intersect(time_remi,tdaily);
%ln5cut=ln5(ia);

%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%lnnorm5cut60=movmean(lnnorm5,[0 ssmooth]);
lnnorm5cut60=movmean(lcnorm5,[0 ssmooth]);
lnnorm5cut60=lc_norm;

%ln5cut30=movmean(ln5cut,[0 360]);



[C,ia,ib] = intersect(time_remi,tdaily);

lnnorm5cut60=lnnorm5cut60(ia);


lnnorm5cut60=lnnorm5cut60';



load pgon_yuc

%%%%%%%%%%%%%%%%       composites             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vdaily60anom=vdaily60-nmean(vdaily60,3);

lc95=prctile(lnnorm5cut60,80)


[pks,pksloc]=findpeaks(lnnorm5cut60)

lc5=prctile(lnnorm5cut60,10)



lnnorm5cut60diff=diff(lnnorm5cut60);


vdailyanom60=vdaily60-nmean(vdaily60,3);

% [eof_maps,pc,expvar]=eof(vdailyanom60);
% 
% for i=1:length(Shedding)
% xline(Shedding(i),'color', [0 .4 0],'linewidth',.65)
% end


vdaily60anom=vdaily60-nmean(vdaily60,3);

ind95=find(lnnorm5cut60>=lc95);
composite95=nmean(vdaily60anom(:,:,ind95),3);

ind5=find(lnnorm5cut60<lc5);
composite5=nmean(vdaily60anom(:,:,ind5),3);



[Csheds,iasheds,ibsheds]=intersect(Shedding,tdaily)




% [Composite95]=intersect(pksloc,ind95)



sheddingv=intersect(Shedding,tdaily)

load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat', 'Shedding')



[pks,pksloc]=findpeaks(lnnorm5cut60,'MinPeakDistance',365)



% figure;plot(tdaily,lnnorm5cut60)
% yyaxis right
% hold on;plot(tdaily,squeeze(nmean(nmean(vdaily60))))
% 
% 
% hold on
% 
% plot(tdaily(Composite95),lnnorm5cut60(Composite95),'r*')
% 


% figure;
% plot(tdaily,lnnorm5cut60,'k')
% hold on
% %plot(tdaily(Composite95),lnnorm5cut60(Composite95),'r*')
% 
% hold on
% plot(tdaily(pksloc),lnnorm5cut60(pksloc),'b*')
% 
% for i=1:length(sheddingv)
%  xline(sheddingv(i),'color', [0 .4 0],'linewidth',.65)
% end
% xlim([tdaily(1) tdaily(end)])
% datetick('x')
% 
% set(gcf,'color','w');



ff=vdaily60anom;

%ff=vdaily-nmean(vdaily,3);

%Si tenes una matriz ff de dim (lt,ly,lx) haces:
%Para hacer estadistica con Matlab sacando los NaNs y luego los pones.

%Si tenes una matriz ff de dim (lt,ly,lx) haces:
ff=(permute(ff,[3 1 2]));

ff=ff(:,:);   %pasa a tener dim  (lt,ly*lx)
mapf=ff(1,:);    % Save geographical dimension to use in mask.
ff(:,any(isnan(ff)))=[];%Remove NaNs
[lt,lxyf_sinNANs]=size(ff);   %dim (lt,lx*ly sin NaNs)



%V33 = clusterdata(ff,'linkage','ward','savememory','on','maxclust',3);

V = clusterdata(ff,'linkage','ward','savememory','on','maxclust',4);

% just to change the numeration of clusters
ind=find(V==3);V(ind)=12;
ind=find(V==4);V(ind)=11;
ind=find(V==2);V(ind)=13;
ind=find(V==1);V(ind)=14;
V=V-10;






%Vkmeans = kmeans(ff, 4);
%V=Vkmeans;

% to check the transitions between 
V1=V;ind=find(V1>1);V1(ind)=NaN;;V1=V1./V1;

for i=1:length(V1)-1;
if V1(i)>0 & V1(i+1)>0;
V1(i+1)=V1(i)+V1(i+1);
V1(i)=NaN;
end
end

V2=V;ind=find(V2>2 | V2<2);V2(ind)=NaN;V2=V2./V2;

for i=1:length(V2)-1;
if V2(i)>0 & V2(i+1)>0;
V2(i+1)=V2(i)+V2(i+1);
V2(i)=NaN;
end
end


V3=V;ind=find(V3>3 | V3<3);V3(ind)=NaN;V3=V3./V3;

for i=1:length(V3)-1;
if V3(i)>0 & V3(i+1)>0;
V3(i+1)=V3(i)+V3(i+1);
V3(i)=NaN;
end
end


V4=V;ind=find(V4>4 | V4<4);V4(ind)=NaN;V4=V4./V4;

for i=1:length(V4)-1;
if V4(i)>0 & V4(i+1)>0;
V4(i+1)=V4(i)+V4(i+1);
V4(i)=NaN;
end
end


% mean time
nmean(V1)

nmean(V2)

nmean(V3)

nmean(V4)




% mean time
nmax(V1)

nmax(V2)

nmax(V3)

nmax(V4)




% transitions
nsum(~isnan(V1))

nsum(~isnan(V2))

nsum(~isnan(V3))

nsum(~isnan(V4))


% 
% VconNaN=mask(V,mapf);
% VconNaN =reshape(VconNaN,301,121);

v1=find(V==1);
v2=find(V==2);
v3=find(V==3);
v4=find(V==4);


% percentage of time in each
length(v1)/length(V)

length(v2)/length(V)

length(v3)/length(V)

length(v4)/length(V)




lnnorm5cut60diff=diff(lnnorm5cut60);lnnorm5cut60diff=[lnnorm5cut60diff(1) lnnorm5cut60diff];


%removing the extreme values associated to sheddin reatach detach
lnnorm5cut60diffs=diff(lnnorm5cut60);lnnorm5cut60diffs=[lnnorm5cut60diffs(1) lnnorm5cut60diffs];

ind=find(lnnorm5cut60diffs>15);lnnorm5cut60diffs(ind)=NaN;

ind=find(lnnorm5cut60diffs<-15);lnnorm5cut60diffs(ind)=NaN;


nmean(lnnorm5cut60diff(v1))

nmean(lnnorm5cut60diff(v2))

nmean(lnnorm5cut60diff(v3))

nmean(lnnorm5cut60diff(v4))


nmean(lnnorm5cut60diffs(v1))

nmean(lnnorm5cut60diffs(v2))

nmean(lnnorm5cut60diffs(v3))

nmean(lnnorm5cut60diffs(v4))




% nmedian(lnnorm5cut60diff(v1))
% 
% nmedian(lnnorm5cut60diff(v2))
% 
% nmedian(lnnorm5cut60diff(v3))
% 
% nmedian(lnnorm5cut60diff(v4))






% % % v5=find(V==5);
% % % v6=find(V==6);
% 
% 
[Composite95,ia,ib]=intersect(sheddingv,tdaily);

Vib=V(ib)

length(v1)/length(V)

length(v2)/length(V)

length(v3)/length(V)

length(v4)/length(V)


%mean current extension asociated to each one
lc_v1=nmean(lnnorm5cut60(v1))

lc_v2=nmean(lnnorm5cut60(v2))

lc_v3=nmean(lnnorm5cut60(v3))

lc_v4=nmean(lnnorm5cut60(v4))



figure
subplot(1,2,2)
hist(V,4); title('all conditions')
subplot(1,2,1)
hist(V(ib),4);title('eddy shed') 
 set(gcf,'color','w');

load pgon_yuc



figure('Renderer', 'painters', 'Position', [100 100 780 640])

subplot(3,2,1:2)
hold on
% xline(Shedding(1),'r','linewidth',.65)
% xline(Detached(1),'color', [0 .4 0],'linewidth',.65)
% xline(Reattached(1),'color', [1.0000    0.8398 0],'linewidth',.65)
% 
% legend('Eddy Shedding','Detached','Reattached','orientation','horizontal','box','off','location','northoutside','Autoupdate','off')

 xline(Shedding(1),'r','linewidth',.65)
legend('Eddy Shedding','orientation','horizontal','box','off','location','north','Autoupdate','off')

for i=1:length(Shedding)
xline(Shedding(i),'r','linewidth',.65)
end


% for i=1:length(Detached)
% xline(Detached(i),'color', [0 .4 0],'linewidth',.65)
% end
% 
% 
% for i=1:length(Reattached)
% xline(Reattached(i),'color', [1.0000    0.8398 0],'linewidth',.65)
% end
xlim([time_merged(1) time_merged(end)])

scatter(tdaily,lnnorm5cut60,18,V,'filled');

colormap(parula(4));
aa=colorbar
ylabel(aa, 'Cluster')


plot(time_merged,lnnorm5cut60,'k'); hold on

xlim([time_merged(1) time_merged(end)])

datetick('x')
%xlim([datenum('1993-1-1') datenum('2022-1-1')])
xlim([time_merged(1) time_merged(end)])

ylabel('LC (km)')
box on; grid on


set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on




xlim([time_merged(1) time_merged(end)])



subplot(3,2,3)
pcolor(X,Z,nmean(vdaily60anom(:,:,v1),3)); hold on;
contour(X,Z,nmean(vdaily60anom(:,:,v1),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
cmocean('balance',22);caxis([-.22 .22]);shading interp;%colorbar
title('Cluster 1')
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});

subplot(3,2,4)
pcolor(X,Z,nmean(vdaily60anom(:,:,v2),3)); hold on;
contour(X,Z,nmean(vdaily60anom(:,:,v2),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
cmocean('balance',22);caxis([-.22 .22]);shading interp;%colorbar
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
%ylabel('Depth (m)');
title('Cluster 2')

subplot(3,2,5)
pcolor(X,Z,nmean(vdaily60anom(:,:,v3),3)); hold on;
contour(X,Z,nmean(vdaily60anom(:,:,v3),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
cmocean('balance',22);caxis([-.22 .22]);shading interp;%colorbar
title('Cluster 3')
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});

subplot(3,2,6)
pcolor(X,Z,nmean(vdaily60anom(:,:,v4),3)); hold on;
contour(X,Z,nmean(vdaily60anom(:,:,v4),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
cmocean('balance',22);caxis([-.22 .22]);shading interp;%colorbar
title('Cluster 4')
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','86.2ºW','','85.8ºW','','85.4ºW','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});

 set(gcf,'color','w');
set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on


savefig('figure4_cluster_february')

print(gcf,'-painters','-depsc2','-r600','figure4_cluster_february')
 print(gcf,'-dpng','-r600','figure4_cluster_february')



x=V;


xx=V;

for i=1:length(xx)-1
if xx(i)==xx(i+1);
xx(i)=NaN;
end
end

ind=isnan(xx);
ii=find(ind==1);xx(ii)=[]

x=xx

m = max(x);
n = numel(x);
y = zeros(m,1);
p = zeros(m,m);
for k=1:n-1
    y(x(k)) = y(x(k)) + 1;
    p(x(k),x(k+1)) = p(x(k),x(k+1)) + 1;
end
p = bsxfun(@rdivide,p,y); p(isnan(p)) = 0;



% 
% figure
% subplot(3,2,1:2)
% scatter(tdaily,lnnorm5cut60,10,V,'filled');cmocean('thermal',4);colorbar
% title('LC extension (km), colors are clusters basted on V at YC'); box on; grid on
% hold on
% for i=1:length(sheddingv)
%  xline(sheddingv(i),'color', [0 .4 0],'linewidth',.65)
% end
% xlim([tdaily(1) tdaily(end)])
% datetick('x')
% 
% 
% 
% subplot(3,2,3)
% pcolor(X,Z,nmean(vdaily60anom(:,:,v1),3)); hold on;
% cmocean('balance','pivot');shading interp;colorbar
% title('1')
% subplot(3,2,4)
% pcolor(X,Z,nmean(vdaily60anom(:,:,v2),3)); hold on;
% cmocean('balance','pivot');shading interp;colorbar
% title('2')
% subplot(3,2,5)
% pcolor(X,Z,nmean(vdaily60anom(:,:,v3),3)); hold on;
% cmocean('balance','pivot');shading interp;colorbar
% title('3')
% subplot(3,2,6)
% pcolor(X,Z,nmean(vdaily60anom(:,:,v4),3)); hold on;
% cmocean('balance','pivot');shading interp;colorbar
% title('4')
% 
%  set(gcf,'color','w');
% 





% rng(1); % For reproducibility
% numStates = 4;
% mc = V
% 
% 
% figure;
% graphplot(mc,'ColorEdges',true);
% 
% figure
% G = graph(V);
% plot(G,'-.dr','NodeLabel',{})
% 
% 
% G = graph(s,t);
% h = plot(G)
ind95=find(lnnorm5cut60>=lc95);
composite95=nmean(vdaily60anom(:,:,ind95),3);


figure;

subplot(2,2,1)

pcolor(X,Z,nmean(vdaily60anom(:,:,ind95),3)); hold on;

pcolor(X,Z,composite95); hold on;
cmocean('balance',11);caxis([-.2 .2])
shading interp
contour(X,Z,composite95,[-.2:.05:.2],'k','showtext','on','Labelspacing',360)
plot(pgon_yuc,'FaceColor',[.7 .7 .7])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');

title('p95 LC')

subplot(2,2,2)
pcolor(X,Z,composite5); hold on;
cmocean('balance',11);caxis([-.2 .2])
shading interp
contour(X,Z,composite5,[-.2:.05:.2],'k','showtext','on','Labelspacing',360)
plot(pgon_yuc,'FaceColor',[.7 .7 .7])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
title(' Correlation between observed Vdaily and LC REMI extension (only p<0.001)')
xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');

title('p5 LC')


subplot(2,2,3)
pcolor(X,Z,nmean(vdaily60anom(:,:,ibsheds),3)); hold on;

cmocean('balance',21);caxis([-.2 .2])
shading interp
hold on
contour(X,Z,nmean(vdaily60anom(:,:,ibsheds),3),[-.2:.05:.2],'k','showtext','on','Labelspacing',360)


%contour(X,Z,nmean(vdaily60anom(:,:,ibsheds),3),[-.2:.05:.2],'k','showtext','on','Labelspacing',360)
plot(pgon_yuc,'FaceColor',[.7 .7 .7])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
title(' Correlation between observed Vdaily and LC REMI extension (only p<0.001)')
xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');

title('Pre shade')



subplot(2,2,4)
pcolor(X,Z,nmean(vdaily60anom(:,:,ibsheds+30),3)); hold on;

cmocean('balance',21);caxis([-.2 .2])
shading interp

hold on
contour(X,Z,nmean(vdaily60anom(:,:,ibsheds+30),3),[-.2:.05:.2],'k','showtext','on','Labelspacing',360)
%contour(X,Z,nmean(vdaily60anom(:,:,ibsheds+30),3),[-.2:.05:.2],'k','showtext','on','Labelspacing',360)
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on
title(' Correlation between observed Vdaily and LC REMI extension (only p<0.001)')
xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
ylabel('Depth (m)');

title('Post shade')







set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on

set(findall(gcf,'-property','FontSize'),'FontSize',12);



% corr(vmid(1:end-llag+1),ln5cut60(llag:end)','rows','pairwise')
% 
% %corr(vmid,ln5cut60','rows','pairwise')
% 
% 
% %vdeep=squeeze(vdaily60(iy,iz,:));
% 
% 
% 
% 
% 
% predictors=[v995cut vmax vmin vmid vmidmin vdeep vdeepmin v1130max v1130min];
% 
% predictors=[v995cut vmax vmin vmid vmidmin vdeep vdeepmin];
% 
% 
% predictand=ln5cut60';
% 
% model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))



% predictors_mooring1=[vmax vmid  v1130min vdeep]
% 
% predictors_mooring2=[vmin vmidmin v1130max vdeepmin]
% 
% 
% %predictors=[v995cut vmax vmid vdeep vdeepmin];
% 
% %predictors=[v995cut predictors_mooring1 predictors_mooring2]
% 
% % model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% % 
% % predictors=[v995cut predictors_mooring1]
% 
% 
% v05cut=v05(ia);
% 
% 
% %predictors=[v995cut vmid  v1130min vdeep]
% 
% 
% predictors=[v995cut vmid vdeep]
% 
% predictors=[v995cut vmid v1130min vdeep]
% 
% predictand=lnnorm5cut60;
% 
% % predictors=[v995cut vmin vmid vmidmin v1130min vdeep vdeepmin]
% 
% model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% 
% train_period=round(length(predictors)/3);
% 
% 
% model1=fitlm(predictors(1:end-lag+1-train_period,:),predictand(lag:end-train_period))
% 
% 
% ypred = predict(model1,predictors);
% 
% 
% 
% figure;
% subplot(2,2,1)
% m_proj('mercator','lon',[longitude(1) longitude(end)], 'lat',[latitude(1) latitude(end)]);hold on
% m_pcolor(longitude,latitude,rho30');shading interp;
% m_contour(longitude,latitude,rho30',[-.5 .5 ],'k','showtext','on','Labelspacing',360)
% 
% cmocean('balance',16)
% caxis([-.8 .8])
% hold on
% m_plot(longitude(ix),latitude(iy),'m+')
% colorbar
% m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
% 
% set(gcf,'color','w');
% 
% title('Pearson correlation between Vgos and LC extension 1 month after (magenta >prctile99.5)')
% set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
% set(findall(gcf,'-property','FontSize'),'FontSize',12);
% m_grid
% 
% 
% 
% subplot(2,2,2)
% pcolor(X,Z,rrho30); hold on;
% 
% cmocean('balance','pivot');%caxis([-.6 .6])
% shading interp
% %axis ij
% %colorbar; axis ij
% contour(X,Z,rrho30,[-1:.1:1],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
% 
% xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
% box on
% 
% 
% %colormap(rednblue)
% %cmocean('balance');
%  %caxis([-120 120])
% %title(' Correlation between observed Vdaily and LC REMI extension (only p<0.001)')
% 
% xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
% yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
% ylabel('Depth (m)');
% 
% set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12);
% 
% 
% 
% minrrho30=find(rrho30==min(min(rrho30)))
% 
% [iy,iz]=find(rrho30==min(min(rrho30)));
% 
% vmin=squeeze(vdaily60(iy,iz,:));
% 
% 
% maxrrho30=find(rrho30==max(max(rrho30)))
% 
% [iy,iz]=find(rrho30==max(max(rrho30)));
% 
% maxrho30v=squeeze(vdaily60(iy,iz,:));
% 
% 
% vmax=squeeze(vdaily60(iy,iz,:));
% 
% 
% rrho30sub=rrho30;rrho30sub(1:72,:)=NaN;
% 
% maxrrho30deep=find(rrho30sub==max(max(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==max(max(rrho30sub)));
% 
% vdeep=squeeze(vdaily60(iy,iz,:));
% 
% 
% rrho30sub=rrho30;rrho30sub(1:30,:)=NaN;rrho30sub(70:end,:)=NaN;
% 
% maxrrho30mid=find(rrho30sub==max(max(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==max(max(rrho30sub)));
% 
% vmid=squeeze(vdaily60(iy,iz,:));
% 
% 
% rrho30sub=rrho30;rrho30sub(1:72,:)=NaN;
% 
% minrrho30deep=find(rrho30sub==min(min(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==min(min(rrho30sub)));
% 
% vdeepmin=squeeze(vdaily60(iy,iz,:));
% 
% 
% rrho30sub=rrho30;rrho30sub(1:30,:)=NaN;rrho30sub(70:end,:)=NaN;
% 
% minrrho30mid=find(rrho30sub==min(min(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==min(min(rrho30sub)));
% 
% vmidmin=squeeze(vdaily60(iy,iz,:));
% 
% 
% 
% rrho30sub=rrho30;rrho30sub(1:55,:)=NaN;rrho30sub(63:end,:)=NaN;
% 
% minrrho301130=find(rrho30sub==min(min(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==min(min(rrho30sub)));
% 
% v1130min=squeeze(vdaily60(iy,iz,:));
% 
% 
% 
% rrho30sub=rrho30;rrho30sub(1:55,:)=NaN;rrho30sub(63:end,:)=NaN;
% 
% maxrrho301130=find(rrho30sub==max(max(rrho30sub)));
% 
% [iy,iz]=find(rrho30sub==min(min(rrho30sub)));
% 
% v1130max=squeeze(vdaily60(iy,iz,:));
% 
% 
% plot(X(minrrho30),Z(minrrho30),'+m','markersize',10,'linewidth',2)
% 
% plot(X(maxrrho30),Z(maxrrho30),'+m','markersize',10,'linewidth',2)
% 
% plot(X(maxrrho30mid),Z(maxrrho30mid),'+m','markersize',10,'linewidth',2)
% 
% plot(X(maxrrho30deep),Z(maxrrho30deep),'+m','markersize',10,'linewidth',2)
% 
% plot(X(minrrho30mid),Z(minrrho30mid),'+m','markersize',10,'linewidth',2)
% 
% plot(X(minrrho30deep),Z(minrrho30deep),'+m','markersize',10,'linewidth',2)
% 
% plot(X(minrrho301130),Z(minrrho301130),'+m','markersize',10,'linewidth',2)
% 
% plot(X(maxrrho301130),Z(maxrrho301130),'+m','markersize',10,'linewidth',2)
% 
% 
% line([-86.1 -86.1],[-1650 -700],'color','g','linewidth',2)
% 
% 
% plot(-86.1 ,Z(maxrrho30deep),'or','markersize',10,'linewidth',2)
% 
% 
% plot(-86.1 ,Z(minrrho301130),'or','markersize',10,'linewidth',2)
% 
% plot(-86.1 ,Z(maxrrho30mid),'or','markersize',10,'linewidth',2)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% hold on
%  plot(time_merged,lnnorm5cut60)
% plot(time_merged,ypred)
% 
% ylabel('LC extension (km)')
% 
% legend('Observed','Predicted','orientation','horizontal','box','off')
% 
% title('Predicted and observed 1 month before')
% datetick('x')
% 
% 
% 
% 
% 
% 
% 
% 
% % oversmooth=1;lag=30;
% % predictors=movmean(predictors,oversmooth,1)
% % predictand=ln5cut60';
% % predictand=movmean(ln5cut,oversmooth)';
% 
% 
% % predictand=predictand(1:30:end);
% % 
% % predictors=predictors(1:30:end,:);
% %modelo=[predictors(1:end-lag+1,[1 4 5]),predictand(lag:end)];
% 
% % 
% % model=fitlm(predictors(1:end-lag+1,1),predictand(lag:end))
% % 
% % model1=fitlm((1:end-lag+1,:),predictand(lag:end))
% 
% 
% 
% 
% 
% % model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% % 
% % model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% % 
% % model1=fitlm(predictors(1:end-lag+1,[1 4 5]),predictand(lag:end))
% % 
% % model1=fitlm(predictors(1:end-lag+1,:),predictand(lag:end))
% % 
% % 
% % model1=fitlm(predictors(1:end-lag+1,[1]),predictand(lag:end))
% % 
% % model1=fitlm(predictors(1:end-lag+1,[1 4 5]),predictand(lag:end))
% % 
% % model1=fitlm(predictors(1:end-lag+1,[1 3 4 5]),predictand(lag:end))
% % 
% % model1=fitlm(predictors,predictand)
% % 
% % model1=fitlm(predictors(1:end-lag+1,2),predictand(lag:end))
% % 
% % model1=fitlm(predictors,predictand)
% % 
% % model1=fitlm(v995(1:end-lag+1),predictand)
% % 
% 
% 
% 
% % rrho30=rho30;
% % %ind=find(pval30>0.001);rrho30(ind)=NaN;
% % 
% % ind=find(rrho30>-.4 & rrho30<.4);rrho30(ind)=NaN;
% % 
% % rrho=rho;
% % ind=find(pval>0.001);rrho(ind)=NaN;
% 
% 
% 
% % 
% % 
% % figure
% % subplot(1,2,1)
% % pcolor(longitude,latitude,rho');shading interp;colorbar
% % 
% % 
% % hold on
% % 
% % 
% % contour(longitude,latitude,rho',[-.4 .4],'k')
% % 
% % 
% % cmocean('balance',13);caxis([-.6 .6])
% % subplot(1,2,2)
% % pcolor(longitude,latitude,rho30');shading interp;colorbar
% % cmocean('balance',13);caxis([-.6 .6])
% %  
% % 
% % 
% % 
% % for i=1:length(longitude)
% % for j=1:length(latitude)
% % 
% % for k=1:365
% % 
% % vv=squeeze(vgos(i,j,1:length(adt)-k+1));
% % lag_corr(i,j,k)=corr(vv,ln5(k:length(adt))');
% % 
% % 
% % end
% % end
% % end
% % 
% % 
% % 
% % 
% % cd scripps_videos
% % 
% % 
% % 
% % 
% % 
% % for i=1:365
% % 
% % figure(i);hold on
% %  set(gcf,'visible','off');
% % 
% % m_proj('mercator','lon',[longitude(1) longitude(end)], 'lat',[latitude(1) latitude(end)]);
% % m_pcolor(longitude,latitude,lag_corr(:,:,i)'); hold on;colorbar
% % 
% % cmocean('balance');caxis([-.5 .5])
% % shading interp
% % 
% % hold on
% % m_contour(longitude,latitude,lag_corr(:,:,i)',[-.5:.1:.5],'k','showtext','on','Labelspacing',360)
% % 
% % %axis ij
% % %colorbar; axis ij
% % %contour(X,Z,pposition,[0:5:50],'k','showtext','on','Labelspacing',360)
% % 
% % % plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% % % xlim([x_grid(1) x_grid(end)]);ylim([-2200 -20]);
% % % box on
% % 
% % title(['Lagged correlation between Vgos and LC extension after ',num2str(i-1), ' days'])
% % 
% % % xticks([-86.8:0.2:-85]);xticklabels({'86.8ºW','86.6ºW','86.4ºW','86.2ºW','85ºW','85.8ºW','85.6ºW','85.4ºW','85.2ºW','85ºW'})
% % % yticks([-2200:200:-200]);yticklabels({'2200','2000','1800','1600','1400','1200','1000','800','600','400','200'});
% % % ylabel('Depth (m)');
% % m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
% % m_grid
% % 
% % set(findall(gcf,'-property','FontSize'),'FontSize',12);set(gcf,'color','w');
% % print(gcf, '-dpng','-r200',[num2str(i)])
% % 
% % 
% % close all
% % end
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% 
% % addpath(genpath('/Users/rlaxenaire/Library/CloudStorage/OneDrive-FloridaStateUniversity/On_Work/Eddy_Detection/TOEddies/Scripts'))
% % addpath(genpath('/Users/rlaxenaire/Library/CloudStorage/OneDrive-FloridaStateUniversity/Boulot/Script_generaux/Matlab'))
% % 
% % %% Directories
% % 
% % Main_Input_dir='/Users/rlaxenaire/Library/CloudStorage/OneDrive-FloridaStateUniversity/Boulot/Data/AVISO/Global/Delayed_Time/DT_2018/SymbolLink';
% % Main_TOEddiesOutput_dir='/Users/gaston/Desktop/Scripps';
% % Main_Dir_Paper='/Users/gaston/Desktop/Scripps/2023_Manta';
% % Main_Paper_Output_dir=fullfile(Main_Dir_Paper,'MAT');
% % Main_Paper_Figure_Output_dir=fullfile(Main_Dir_Paper,'Figure');
% % 
% % 
% % addpath(genpath(Main_Dir_Paper))
% 
% 
% % load mask_HYCOM_25.mat
% % load LC_VmaxContour_V2.mat
% % load LC_mat_VmaxContour_V2_Shed_TOEddies.mat
% % load LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat
% % 
% % 
% % 
% % 
% % 
% % lc=LC_VmaxContour(:,8);
% % aa=cell2mat(lc);
% % lcaa=aa(:,2);
% % aa=cell2mat(lc);
% % aa=aa(:,1);
% % 
% % %lc=(lc(:,8));
% % %% Parameters TOEddies
% % 
% % %%%% Step by Step
% % % Name the output files according to their type (sla or roms etc) and
% % % the date
% % type_SSH='ADT';
% % output_name = [lower(type_SSH),'_'];
% % 
% % % Degree of ghost point if needed (if WL=EL)
% % ghost_point_deg= 3;
% % 
% % % Define area of interest. Both -180:180 and 0-360 longitudes
% % % convention are possible as the best is used according to the
% % % area of interest
% % NL = 31 ;
% % SL = 17;
% % WL = 261;
% % EL = 284;
% % 
% % % Minimum amplitude requested to detec eddies [m]
% % min_amp=1e-3;
% % 
% % % Number maximum of eddies
% % max_eddies=1e2;
% % 
% % % Parameters trajectories
% % min_intersect_all=[0,50];% Minimum percent of area superposition neededed to follow a trajectory
% % t_max_asso=5;% Max extension of map tracking
% % type_contour_all={'_max','_out'};% Type of contour used in the superposition
% % type_eddy_all={'Anti','Cyclo'};
% % 
% % % Parameters filtering trajectories
% % min_lifetime=7;  % Minimum lifetime for trajectories that are not filtered out
% % min_med_ampl=0; % Minimum mediane amplitude for trajectories that are not filtered out
% % 
% % %%%% Figures
% % type_eddy_full={'Anticyclonic','Cyclonic'};
% % %date_interest=datenum():datenum();
% % 
% % %% Parameters LC
% % 
% % %%%%%%% Parameters LC
% % YC=[-87.25, 21.25; -84, 22.15]; % X and Y cordinates of the Yucatan Channel
% % FS=[-80.6, 25.5; -80.6,22.75];% X and Y cordinates of the Straits of Florida
% % SSH_resolution = 5e-3; % SSH resolution to define a contour [m]
% % Velocity_resolution=[40,50]; % Maximum number of contours tested
% % IsRegularGrid=1; % Assume irregular grid
% % 
% % % Sections
% % Name_Sections={'Yucatan','FlStr'};
% % Vec_Sections={YC,FS};
% % 
% % 
% % %%%%%%% Phases
% % Mat_Values_Phases=[2.2,3];
% % 
% % % Create TOpo
% % 
% % dir_etopo='/Users/gaston/Documents/MATLAB/toolboxes_ocean/m_map/etopo1_ice_c_i2';
% % samplefactor=1;
% % Bathy_file='/Users/gaston/Documents/MATLAB/toolboxes_ocean/m_map/etopo1_ice_c_i2/etopo1_ice_c_i2.bin';
% % 
% % Rad_Earth= 6371000;    %  [m] Mean radius of earth  A.E.Gill
% % 
% % 
% % % Define Mask GoM
% % X_Lim=[262.120002746582,262.160003662109,270.879997253418,271.959999084473,275.9604,276.0101,277.7521,279.6500,279.040000915527,277.839996337891,262.120002746582];
% % Y_Lim=[30.6964759826660,18.1296672821045,18.3576030731201,21.234476931611,21.234476931611,22.3293979432848,22.8565821340907,22.8565821340907,26.0168933868408,31.1768054962158,30.6964759826660]';
% % 
% % 
% % %%%% Cuba
% % 
% % SL_topo=min([Vec_Sections{1}(:,2);Vec_Sections{2}(:,2)])-.5;
% % NL_topo=max([Vec_Sections{1}(:,2);Vec_Sections{2}(:,2)])+.5;
% % 
% % WL_topo=min([Vec_Sections{1}(:,1);Vec_Sections{2}(:,1)])-.5;
% % EL_topo= max([Vec_Sections{1}(:,1);Vec_Sections{2}(:,1)])+.5;
% % 
% % [Bathy, ~] = etopo([dir_etopo,'/etopo1_ice_c_i2.bin'], samplefactor ...
% %     ,[SL_topo NL_topo] ...
% %     ,[WL_topo EL_topo]);
% % 
% % Lon_topo=linspace(WL_topo,EL_topo,size(Bathy,2))';
% % Lat_topo=linspace(SL_topo,NL_topo,size(Bathy,1))';
% % 
% % % Find limit of Cuba
% % cout = contourc(Lon_topo,Lat_topo,Bathy,[0 0]);
% % cstruct_SecCuba = struct('x', cell(1,1), 'y', cell(1,1));
% % count_tmp=0;idxlabel = 1;
% % while idxlabel <= size(cout,2)
% %     n = cout(2,idxlabel);
% %     if n>count_tmp
% %         cstruct_SecCuba(1).x = cout(1,idxlabel+(1:n));
% %         cstruct_SecCuba(1).y = cout(2,idxlabel+(1:n));
% %         count_tmp=n;
% %     end
% %     idxlabel = idxlabel+n+1;
% % end
% % 
% % % Create a vector containing the 2 sections and Cuba
% % [X0_Yuc,Y0_Yuc,I,~] = my_intersections(cstruct_SecCuba(1).x,cstruct_SecCuba(1).y,Vec_Sections{1}(:,1),Vec_Sections{1}(:,2));
% % [~,id_intrest]=min(X0_Yuc);
% % I_Yuc=floor(I(id_intrest));
% % X0_Yuc=X0_Yuc(id_intrest);
% % Y0_Yuc=Y0_Yuc(id_intrest);
% % 
% % [X0_Str,Y0_Str,I,~] = my_intersections(cstruct_SecCuba(1).x,cstruct_SecCuba(1).y,Vec_Sections{2}(:,1),Vec_Sections{2}(:,2));
% % [~,id_intrest]=min(Y0_Str);
% % I_Str=ceil(I(id_intrest));
% % X0_Str=X0_Str(id_intrest);
% % Y0_Str=Y0_Str(id_intrest);
% % 
% % cstruct_SecCuba_all=cstruct_SecCuba;
% % 
% % cstruct_SecCuba(1).x=[X0_Str,cstruct_SecCuba(1).x(I_Str:I_Yuc),X0_Yuc];
% % cstruct_SecCuba(1).y=[Y0_Str,cstruct_SecCuba(1).y(I_Str:I_Yuc),Y0_Yuc];
% % 
% % %         % Verify
% % %         figure
% % %         hold on
% % %         pcm(Lon_topo,Lat_topo,Bathy);shading flat
% % %         plot(cstruct_SecCuba_all(1).x,cstruct_SecCuba_all(1).y,'k')
% % %         plot(Vec_Sections{1}(:,1),Vec_Sections{1}(:,2),'r')
% % %         plot(Vec_Sections{2}(:,1),Vec_Sections{2}(:,2),'r')
% % %
% % %         plot(cstruct_SecCuba(1).x,cstruct_SecCuba(1).y,'g')
% % 
% % X_Cuba=[285,285,275,275,285];
% % Y_Cuba=[17,24,24,17,17]';
% % 
% % %% Parameters LCE
% % 
% % min_time_LCE_in_GoM=7;
% % 
% % %% Flags
% % 
% % id_TOEddies_recompute=0;
% % id_LCFront_recompute=0;
% % id_LCE_recompute=0;
% % 
% % id_plot_map_LC=0;
% % id_plot_hist_LC=0;
% % 
% % 
% % %% Figure 1 : L_{LC} overlap and hist
% % 
% % Dl_cells=.1;
% % Lim_cont=1:Dl_cells:5;
% % 
% % if id_plot_map_LC || id_plot_hist_LC
% %     Dir_figure_TMP=[Main_Paper_Figure_Output_dir,'/LLC_Map'];
% %     if exist(Dir_figure_TMP,'dir')==0;
% %         mkdir(Dir_figure_TMP);
% %     end;
% %     Colormap=jet(length(Lim_cont)-1);
% %     Colormap(1,:)=[0 0 0];
% % end
% % 
% % BD_LC=cellfun(@(X) X(2),LC_VmaxContour(:,8));
% % Ratio_dist_LC=cellfun(@(X) X(1)./X(2),LC_VmaxContour(:,8));
% % Time_LC=[LC_VmaxContour{:,1}];
% % Time_LC_diff_year=(Time_LC-Time_LC(1)+1)/365.25;
% % Area_VmaxContour=cellfun(@(X) X(1),LC_VmaxContour(:,7));
% % 
% % Rcorr=corr(Area_VmaxContour(:,1),Ratio_dist_LC);
% % p = polyfit(Ratio_dist_LC,Area_VmaxContour(:,1),1);
% % yfit =  p(1) * Ratio_dist_LC + p(2);
% % yresid = Area_VmaxContour(:,1) - yfit;
% % SSresid = sum(yresid.^2);
% % SStotal = (length(Area_VmaxContour(:,1))-1) * var(Area_VmaxContour(:,1));
% % R2= 1 - SSresid/SStotal;
% % disp([' Correlation coefficient of Norm length LC with area is ',num2str(round(1e2*Rcorr)/1e2),' and R^2 ',num2str(round(1e2*R2)/1e2)])
% % 
% % disp('%%%%')
% % disp([' Mean (Median) LC is', num2str(mean(Ratio_dist_LC)), '(',num2str(median(Ratio_dist_LC)),')  '...
% %     'Ratio LC perc<=',num2str(Mat_Values_Phases(1)),'=',num2str(round(sum(Ratio_dist_LC<=Mat_Values_Phases(1))*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(Mat_Values_Phases(2)),'=',num2str(round(sum(Ratio_dist_LC>=Mat_Values_Phases(2))*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(2.9),'=',num2str(round(sum(Ratio_dist_LC>=2.9)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3),'=',num2str(round(sum(Ratio_dist_LC>=3)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3.2),'=',num2str(round(sum(Ratio_dist_LC>=3.2)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3.25),'=',num2str(round(sum(Ratio_dist_LC>=3.25)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3.3),'=',num2str(round(sum(Ratio_dist_LC>=3.3)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3.5),'=',num2str(round(sum(Ratio_dist_LC>=3.5)*100/length(Ratio_dist_LC))),...
% %     ' perc>=',num2str(3.6),'=',num2str(round(sum(Ratio_dist_LC>=3.6)*100/length(Ratio_dist_LC))),...
% %     ' Mean +/- Std BD>=',num2str(mean(BD_LC)),'+/-',num2str(std(BD_LC)),' km',...
% %     ])
% % 
% % disp([' Not connected =',num2str(round(sum(Id_non_continues)./size(LC_VmaxContour,1)*10000)/100),' %'])
% % 
% % 
% %  Colormap=parula(length(Lim_cont));
% % 
% % % Colormap=cmocean('parula',length(Lim_cont));
% % 
% % %if id_plot_map_LC
% % 
% %     [~,id_sort_LC]=sort(Ratio_dist_LC,'descend');
% % 
% %     [~,id_sort_LC]=sort(aa,'descend');
% % 
% % 
% %     % Map
% %     figure%('DefaultAxesFontSize',14,'DefaultAxesFontWeight','bold')
% %     hold on
% %     %Colormap=colormap(parula(21))
% % m_proj('mercator','lon',[longitude(1) longitude(end)], 'lat',[latitude(1) latitude(end)]);hold on
% %   
% % 
% % hold on
% % % m_proj('mercator','lon',[WL-360 EL-360], 'lat',[SL NL]);
% %     for loop_tmp=id_sort_LC'
% %         %             if isnan(LC_VmaxContour{loop_tmp,6})
% %         %                 continue
% %         %             end
% %         id_min=find(Lim_cont(2:end)>Ratio_dist_LC(loop_tmp),1);
% %         if isempty(id_min);id_min=size(Colormap,1);end;
% % 
% %         xx=LC_VmaxContour{loop_tmp,2};xx(xx>180)=xx(xx>180)-360;
% % 
% % yy=LC_VmaxContour{loop_tmp,3};ind=find(yy<21.25);yy(ind)=NaN;
% % 
% % 
% %         m_plot(xx,yy,'Color',Colormap(id_min,:));
% %        m_plot(xx,yy,'Color',Colormap(id_min,:));
% %     end
% % 
% % hold on 
% % [CS,CH]=m_etopo2('contour',[-7000:1000:-1000 -200 0],'color',[.7 .7 .7]);
% % 
% % m_contour(longitude,latitude,adtm',[.17 .17],'r','linewidth',2)
% % m_line([Lon(1) Lon(end)],[Lat(1) Lat(end)],'color','g','linewidth',2)
% % 
% % 
% % 
% %  set(gcf,'color','w');
% %  set(findall(gcf,'-property','FontSize'),'FontSize',12);
% % 
% % ax=m_contfbar([.455 .62],0.21,[650:200:2800],[650:200:2800],'endpiece','no','axfrac',.02);
% % ax.XLabel.String='Extension (km)';
% % 
% % m_gshhs_l('patch',[0.81818 0.77647 0.70909]);
% %  m_grid('linestyle','none');
% % 
% % 
% % print(gcf,'-painters','-depsc2','-r600','fig1c_yuc')
% % 
% % 
% % 
% % 










