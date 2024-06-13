%%--plot difference--%
load ./Figs12/12NewNe
rne1=ne_background-ne_mart0;
rne2=ne_background-ne_mart1;

figure
xx=lat(1:Nlat);
yy=lon(1:Nlon);
zz=h(1:Nh);

[x,y,z]=meshgrid(xx,yy,zz);
% xslice=[xx(1),xx(9),xx(11),xx(13),xx(end)];
xslice=[xx(1),xx(10),xx(end)];
yslice=yy(end);
zslice=h(1);
slice(x,y,z,rne1,xslice,yslice,zslice)
xlim([lat(1) lat(end)])
ylim([lon(1) lon(end)])
zlim([h(1) h(end)])
xlabel('Lat/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Lon/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
zlabel('H/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
set(gca,'Ydir','Normal','fontsize',12,'fontname','Times New Roman','fontweight','bold')
shading interp
colorbar

% slice 画法
figure
% xslice=[xx(1),xx(9),xx(11),xx(13),xx(end)];
xslice=[xx(1),xx(10),xx(end)];
yslice=yy(end);
zslice=h(1);
slice(x,y,z,rne2,xslice,yslice,zslice)
xlim([lat(1) lat(end)])
ylim([lon(1) lon(end)])
zlim([h(1) h(end)])
xlabel('Lat/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Lon/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
zlabel('H/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
set(gca,'Ydir','Normal','fontsize',12,'fontname','Times New Roman','fontweight','bold')
shading interp
% title(['lamda=0.5 k=500 MART','交点覆盖率=',num2str(cover_rate)],'fontsize',12,'fontname','Times New Roman','fontweight','bold')

colorbar
