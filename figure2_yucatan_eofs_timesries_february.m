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


nmean(ln5)
nstd(ln5)

nmean(lcnorm5)
nstd(lcnorm5)





corr(ln5(1:length(lcnorm5))',lcnorm5)


load vmerged.mat
vdaily=vmerged; tdaily=time_merged;

ssmooth=30;

llag=30;

vdaily60=movmean(vdaily,[0 ssmooth],3);

vdailyanom60=vdaily60-nmean(vdaily60,3);

[eof_maps,pc,expvar]=eof(vdailyanom60);
% 
% for i=1:length(Shedding)
% xline(Shedding(i),'color', [0 .4 0],'linewidth',.65)
% end
% 

%[eof_maps,pc,expv] = eof(vdaily60);

% scalinig eof pc -1 to 1
for k = 1:size(pc,1)

   % Find the index of the maximum value in the time series:
   [maxval,ind] = max(abs(pc(k,:)));

   % Divide the time series by its maximum value:
   pc(k,:) = pc(k,:)/maxval;

   % Multiply the corresponding EOF map:
   eof_maps(:,:,k) = eof_maps(:,:,k)*maxval;

end

load pgon_yuc






figure
ylabel('Depth (m)');
hold on
for i=1:3%length(eof_maps(1,1,:))

%subplot(2,round(length(eof_maps(1,1,:))/2),i)

subplot(4,3,i+6)

pcolor(X,Z,-eof_maps(:,:,i));shading interp
hold on
contour(X,Z,-eof_maps(:,:,i),[-1:.1:1],'k','showtext','on','Labelspacing',360)
%cmocean('balance','pivot') 
cmocean('balance',12*2) 

caxis([-.6 .6])

title(['EOF ',num2str(i),' (',num2str(round(expvar(i)),'%0.1f'),'%)'])
%colorbar
hold on

%plot(pgon_yuc,'FaceColor',[.7 .7 .7])
plot(pgon_yuc,'FaceColor',[0.81818 0.77647 0.70909])


xlim([x_grid(1) x_grid(end)]);ylim([-2200 -10]);
box on

xticks([-86.8:0.2:-85]);xticklabels({'','86.6ºW','','','','85.8ºW','','','','85ºW'})
yticks([-2200:200:-200]);yticklabels({'2200','','1800','','1400','','1000','','600','','200'});
end


subplot(4,3,4:6)

plot(tdaily,-pc(1:3,:));axis tight%shading interp
datetick('x')
legend('PC1','PC2','PC3','orientation','horizontal','box','off')
grid on
axis tight
ylabel('PCs')

subplot(4,3,1:3)

plot(time,movmedian(ln1,5),'k'); hold on
%ylim([500 3000]);yyaxis right
plot(time_remi,movmedian(lc_norm,5),'b'); hold on
ylim([500 3000])
hold on

xline(Shedding(1),'r','linewidth',.65)
xline(Detached(1),'color', [0 .4 0],'linewidth',.65)
xline(Reattached(1),'color', [1.0000    0.8398 0],'linewidth',.65)

legend('0.17m','< V >','Eddy Shedding','Detached','Reattached','orientation','horizontal','box','off','location','northoutside','Autoupdate','off')

xticks([]);xticklabels({})



for i=1:length(Shedding)
xline(Shedding(i),'r','linewidth',.65)
end


for i=1:length(Detached)
xline(Detached(i),'color', [0 .4 0],'linewidth',.65)
end


for i=1:length(Reattached)
xline(Reattached(i),'color', [1.0000    0.8398 0],'linewidth',.65)
end


datetick('x')
%xlim([datenum('1993-1-1') datenum('2022-1-1')])
xlim([time_merged(1) time_merged(end)])

ylabel('LC extension (km)')
box on; grid on


set(gcf,'color','w');set(findall(gcf,'-property','FontWeight'),'FontWeight','Normal');box on
% 
% savefig('figure2_eofs_february')
% print(gcf,'-painters','-depsc2','-r600','figure2_eofs_february')

