
load 00NewNe.mat
[Nlat, Nlon, Nh]=size(ne_background);
Ne=reshape(ne_background,Nlat*Nlon*Nh,1);
ne0=reshape(ne_ini,Nlat*Nlon*Nh,1);
Ne_mart0=reshape(ne_mart0,Nlat*Nlon*Nh,1);
Ne_mart1=reshape(ne_mart1,Nlat*Nlon*Nh,1);

figure
plot(Ne)
hold on
plot(ne0)
hold on
plot(Ne_mart0)
hold on
plot(Ne_mart1)
% set(gca,'Ylim')
xlabel('n_{voxel}','fontsize',14,'fontname','Times New Roman','fontweight','bold')
ylabel('n_{e}/el/m^3','fontsize',14,'fontname','Times New Roman','fontweight','bold')
% title('MART','fontsize',12,'fontname','Times New Roman','fontweight','bold')
legend('Background', 'initial', 'MART0', 'MART1')
set(gca,'fontsize',14,'fontname','Times New Roman','fontweight','bold')
grid on

neb=reshape(ne_background(5,5,:), [],1);
neb=smooth(neb);
nem0=reshape(ne_mart0(5,5,:), [],1);
nem0=smooth(nem0);
nem1=reshape(ne_mart1(5,5,:), [],1);
nem1=smooth(nem1);
nei=reshape(ne_ini(5,5,:), [],1);
nei=smooth(nei);

hh=100:25:1000-25;
% Re_Error=mean(abs(Ne_mart(2000:12000)-Ne(2000:12000))./Ne(2000:12000));
figure
% plot(neb,hh, nei,hh,nem0,hh, nem1,hh,  'linewidth', 1.5)
plot(neb,hh, nem0,hh, nem1,hh,  'linewidth', 1.5)
xlabel('n_{e}el/m^3','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Altitude/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
legend('Background', 'initial', 'MART0', 'MART1')
grid on