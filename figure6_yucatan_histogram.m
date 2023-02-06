


load vort_merged.mat

load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat', 'Shedding')



[Csheds,iasheds,ibsheds]=intersect(Shedding,timee);

load vort_merged.mat

vort_anom=vort_merged-nmean(vort_merged,3);

load('PresenceEddies4Gaston.mat')

cyclones=Mat_Presence_Vmax;ind=find(cyclones>0);cyclones(ind)=0;

anticyclones=Mat_Presence_Vmax;ind=find(anticyclones<0);anticyclones(ind)=0;


load('/Users/gaston/Desktop/Scripps/Remi_scripts_LC/MAT/LC_mat_VmaxContour_V2_Shed_TOEddies_Full.mat')

load lc_rem
lcnorm5=movmedian(lc_norm,5);

[Csheds,iasheds,ibsheds]=intersect(Shedding,time_remi)

lc_sheddings=lc_norm(ibsheds)

% figure;
% subplot(1,3,1)
% pcolor(X,Y,nmean(Mat_Presence_Vmax,3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% subplot(1,3,2)
% pcolor(X,Y,nmean(anticyclones,3)');shading interp;cmocean('balance',8,'pivot'); colorbar
% 
% subplot(1,3,3)
% pcolor(X,Y,nmean(cyclones,3)');shading interp;cmocean('balance',8,'pivot'); colorbar
% 
% figure;pcolor(X,Y,nmean(vort_merged,3)');shading interp
% cmocean('balance','pivot')
% 

% 
% xx = 650:50:2800;
% lcround = interp1(xx, xx, lcnorm5, 'nearest', 'extrap')
% 
% lcsheddings = interp1(xx, xx, lc_sheddings, 'nearest', 'extrap')
% 
% 

[C,ia,ibb] = intersect(Detached,time_remi)

nmean(lc_norm(ibb-1))


[C,ia,ib] = intersect(Shedding,time_remi)

nmean(lc_norm(ib-1))



xx = 650:100:2800;
lcround = interp1(xx, xx, lc_norm, 'nearest', 'extrap')
lcsheddings = interp1(xx, xx, lc_norm(ib-1), 'nearest', 'extrap')

lcdetachments = interp1(xx, xx, lc_norm(ibb-1), 'nearest', 'extrap')


lcroundd=lcround;



numIntervals = 23;
intervalWidth = 100;
x = 550:100:2750;
y=lcround
y1=lc_norm(ib-1)
y1d=lc_norm(ibb-1)


%Next, use the HISTC function to find the frequency of each data range "x" in the given data set "y". This function returns the histogram count for a data set and range.
ncount = histc(y ,x);
%Calculate the relative frequency of each data range by dividing the frequency by the total number of data points:
relativefreq = ncount/length(y);
%Finally plot the relative frequency versus the data ranges as a bar chart. On this chart, the bars will be adjoining, and the tick marks on the x-axis will label the extent of each bar's data range.



%Next, use the HISTC function to find the frequency of each data range "x" in the given data set "y". This function returns the histogram count for a data set and range.
ncount1 = histc(y1 ,x);
relativefreq1 = ncount1/length(y1);


%Next, use the HISTC function to find the frequency of each data range "x" in the given data set "y". This function returns the histogram count for a data set and range.
ncount1d = histc(y1d ,x);
relativefreq1d = ncount1d/length(y1d);










load vmerged


vdaily=vmerged; tdaily=time_merged;

ssmooth=21;

llag=30;

vdaily60=movmean(vdaily,[0 ssmooth],3);

vdailyanom60=vdaily60-nmean(vdaily60,3);


[C,ia,ib] = intersect(tdaily,time_merged);

cyclones=-cyclones(:,:,ia);timee

cyclones_anom=cyclones-nmean(cyclones,3);

[Ca,iaa,iba] = intersect(Shedding,timee);

cyclones_anomm=movmean(cyclones_anom,[0 ssmooth],3);

[Cd,iad,ibd] = intersect(Detached,timee);


figure
% subplot(2,3,1)

hold on

bar(x-intervalWidth/2, relativefreq1d,1,'r','FaceAlpha', 0.5)

bar(x-intervalWidth/2, relativefreq1,1,'facecolor', [0 .4 0],'FaceAlpha', 0.5)
bar(x-intervalWidth/2, relativefreq,1,'b','FaceAlpha', 0.5)

xlim([min(x) max(x)])
set(gca, 'xtick', x)

xticks([600:400:2600])

ylabel('Relative frecuency')
xlabel('Loop Current extension (km)')
grid on; box on

legend('Sheding','Dettachment','All days','location','northwest')
set(gcf,'color','w');

% 
% vanom=vmerged-nmean(vmerged,3);
% 
% subplot(2,3,2)
% [Cs,ias,ibs] = intersect(Shedding,time_merged);
% 
% pcolor(x_grid,Z(:,1),nmean(vdailyanom60(:,:,ibs-1),3));shading interp;cmocean('balance',21,'pivot')
% colorbar
% 
% % pcolor(x_grid,Z(:,1),nmean(vanom(:,:,ibs-1),3));shading interp;cmocean('balance',21,'pivot')
% % colorbar
% 
% %caxis([-0.5 0.5])
% 
% 
% subplot(2,3,3)
% 
% [Cs,ias,ibs] = intersect(Detached,time_merged);
% 
% pcolor(x_grid,Z(:,1),nmean(vdailyanom60(:,:,ibs-1),3));shading interp;cmocean('balance',21,'pivot')
% colorbar
% 
% % pcolor(x_grid,Z(:,1),nmean(vanom(:,:,ibs-1),3));shading interp;cmocean('balance',21,'pivot')
% % colorbar
% % 
% 
% 
% subplot(2,3,4)
% pcolor(X,Y,nmean(cyclones,3)');shading interp;cmocean('balance',16)
% caxis([-0.8 0.8])
% 
% subplot(2,3,5)
% pcolor(X,Y,nmean(cyclones_anomm(:,:,[iba]),3)');shading interp;cmocean('balance',11)
% caxis([-0.2 0.2])
% 
% 
% subplot(2,3,6)
% pcolor(X,Y,nmean(cyclones_anomm(:,:,[ibd]),3)');shading interp;cmocean('balance',11)
% caxis([-0.2 0.2])
% 
% 
% 
% 
% 
% 
% 


% 
% ws1(i)=ws;
% nr1(i)=nr;
% ln1(i)=ln;
% ar1(i)=ar;
% 
% 
% 
% plot(ln1,nr1,'.')
% hold on
% yyaxis right
% plot(ln1,ws1,'.')
% 
% xlim([500 2800])






load pgon_yuc

[C,ia,ib] = intersect(Shedding,time_merged);


% pcolor(X,Z,nmean(vdailyanom60(:,:,ib-1),3)); shading interp
% caxis([-.1 .1])
% cmocean('balance',11)
% hold on
% contour(X,Z,nmean(vdailyanom60(:,:,ib-1),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
% 
% 
% 
% 
% figure
% plot(time_merged,squeeze(nmean(nmean(vmerged))))
% hold on
% 
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
% 
% xlim([time_merged(1) time_merged(end)])
% 
% datetick('x')


% figure
% 
% hold on
% 
% bar(x-intervalWidth/2, relativefreq1d,1,'r','FaceAlpha', 0.5)
% 
% bar(x-intervalWidth/2, relativefreq1,1,'facecolor', [0 .4 0],'FaceAlpha', 0.5)
% bar(x-intervalWidth/2, relativefreq,1,'b','FaceAlpha', 0.5)
% 
% xlim([min(x) max(x)])
% set(gca, 'xtick', x)
% 
% xticks([600:200:2600])
% ylabel('Relative frecuency')
% xlabel('Loop Current extension (km)')
% grid on; box on
% 
% legend('Sheding','Dettached','All days','location','northwest')
% set(gcf,'color','w');
% 
% 
% subplot(1,3,2)
% load pgon_yuc
% 
% [C,ia,ib] = intersect(Shedding,time_merged);
% 
% 
% pcolor(X,Z,nmean(vdailyanom60(:,:,ib-1),3)); shading interp
% caxis([-.1 .1])
% cmocean('balance',11)
% hold on
% contour(X,Z,nmean(vdailyanom60(:,:,ib-1),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
% 
% subplot(1,3,3)
% load pgon_yuc
% 
% [C,ia,ibb] = intersect(Detached,time_merged);
% 
% 
% pcolor(X,Z,nmean(vdailyanom60(:,:,ibb-1),3)); shading interp
% caxis([-.1 .1])
% cmocean('balance',11)
% hold on
% contour(X,Z,nmean(vdailyanom60(:,:,ibb-1),3),[-1:.1:1],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Mat_Presence_Vmax=movmean(Mat_Presence_Vmax,[0 30],3);
% 
% Mat_Presence_Vmax_anom=Mat_Presence_Vmax-nmean(Mat_Presence_Vmax,3);
% 
% Mat_Presence_Vmax_anom=Mat_Presence_Vmax-nmean(Mat_Presence_Vmax,3);
% 
% cyclones_anom=cyclones-nmean(cyclones,3);
% 
% figure;
% subplot(1,2,1)
% pcolor(X,Y,nmean(Mat_Presence_Vmax_anom(:,:,ib-1),3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% subplot(1,2,2)
% pcolor(X,Y,nmean(Mat_Presence_Vmax_anom(:,:,ibb-1),3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% 
% 
% figure;
% 
% subplot(1,3,1)
% pcolor(X,Y,nmean(cyclones,3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% subplot(1,3,2)
% pcolor(X,Y,nmean(cyclones(:,:,ib-1),3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% subplot(1,3,3)
% pcolor(X,Y,nmean(cyclones(:,:,ibb-1),3)');shading interp;cmocean('balance',16,'pivot'); colorbar
% 
% 
% 
% 
% 
% 
% pcolor(X,Y,nmean(anticyclones,3)');shading interp;cmocean('balance',8,'pivot'); colorbar
% 
% subplot(1,3,3)
% pcolor(X,Y,nmean(cyclones,3)');shading interp;cmocean('balance',8,'pivot'); colorbar
% 
% figure;pcolor(X,Y,nmean(vort_merged,3)');shading interp
% cmocean('balance','pivot')
% 
% 
% ind=find(time_merged==(datenum('2015-9-25')))
% 
% load pgon_yuc
% 
% figure
% 
% subplot(1,2,1)
% pcolor(X,Z,nmean(vmerged(:,:,ind),3)); shading interp
% caxis([-.1 1.4])
% cmocean('balance')
% hold on
% contour(X,Z,nmean(vmerged(:,:,ind),3),[-1:.1:1.4],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
% subplot(1,2,2)
% pcolor(X,Z,nmean(vmerged,3)); shading interp
% caxis([-.1 1.4])
% cmocean('balance')
% hold on
% contour(X,Z,nmean(vmerged,3),[-1:.1:1.4],'k','showtext','on','Labelspacing',360)
% 
% plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])
% 
