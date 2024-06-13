clear 
close all
clc
%%%%%%%%%%%%%%%%%%       AIS+GNSS CT迭代  AIS辅助反演 20201030
%%用IRI数据作输入 精简版
%%乘法代数重建算法%%

tic
%%%%%%%%%%%%%文件载入
load Neb2022031512             %背景电离层电子浓度
load STEC
load net

Xg=lat;
Yg=lon;
Zg=h;        %单位km

Nlat=length(lat)-1;
Nlon=length(lon)-1;
Nh=length(h)-1;
dh=h(2)-h(1);
ne_background=ne_background(1:Nlon, 1:Nlat, 1:Nh);
Ne=reshape(ne_background,Nlat*Nlon*Nh,1);

Re=6378135; %地球半径 单位 m%

%IRI作初值
load Ne003
ne_ini=ne0(1:Nlon, 1:Nlat, 1:Nh);
ne00=reshape(ne_ini,Nlat*Nlon*Nh,1);

% ne0=ones(Nlat*Nlon*Nh,1);
% ne_ini=reshape(ne0, Nlon, Nlat, Nh);
% %随机初值，全10^11
% ne_ini=ones(14, 14, 14)*10^11;
% ne0=reshape(ne_ini,size(ne_ini,1)*size(ne_ini,2)*size(ne_ini,3),1);

load RaysS
%%--MART 系数矩阵--
load MART_coefficient_matrix_sp
% [Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg);

% D1 = Grid_raysinfo.MART_coefficient_matrix_3index;

D=MART_coefficient_matrix;
% D(all(D==0,2),:)=[];

%%%%%%% method 1 %%%%%%%%%%%
% dis1=sum(D(10,:));   %求斜射线长度
% dis2=h(end)-h(1);           %垂直高度
% VTEC=STEC(1)*(cos(dis2/dis1));      %求测量值的VTEC
% VTEC0=sum(ne_ini(1,1,:)*dh); %均匀分布的网格求初始值的VTEC
% ne0=ne0*delta;      %修正初值

%%%%%% method 2 %%%%%%%%%%%
STEC0=D*ne00;
delta=mean(STEC./STEC0);     %求修正因子
ne01=delta*ne00;      %修正初值

%乘法代数重建算法
tic
lamdak=1.5;       %驰豫因子
Ne_mart0=ne00;
[Nray,Nnet]=size(D);
kn=10;
NUM=sum(D,1);
NonZeroCross=find(NUM~=0);       %找到被覆盖的网格
cover_rate=length(NonZeroCross)/Nnet;

tic
for j=1:kn            %迭代循环     
    j
    xx=Ne_mart0;
    for i=1:Nray
        DD=norm(D(i,:));
        ss=STEC(i)/dot(Ne_mart0,D(i,:));
        for l=1:Nnet
            Ne_mart0(l)=Ne_mart0(l)*(ss^(lamdak*D(i,l)/DD));
        end 
    end
    xxx=Ne_mart0;
    SD(j)=rms(abs(xxx-xx));    
    
end 
%--计算误差等--
ne_mart0=reshape(Ne_mart0,Nlat,Nlon,Nh);

Re_Error0=mean(abs(Ne_mart0-Ne)./Ne);
Re_Error2=mean(abs(ne00-Ne)./Ne);
ad0=mean(abs(Ne_mart0-Ne));
ad2=mean(abs(ne00-Ne));
rmse0=sqrt(mean((Ne_mart0-Ne).^2));
rmse2=sqrt(mean((ne00-Ne).^2));


%--新的初始值--%
Ne_mart1=ne01;
tic
for j=1:kn            %迭代循环     
    j
    xx=Ne_mart1;
    for i=1:Nray
        DD=norm(D(i,:));
        ss=STEC(i)/dot(Ne_mart1,D(i,:));
        for l=1:Nnet
            Ne_mart1(l)=Ne_mart1(l)*(ss^(lamdak*D(i,l)/DD));
        end 
    end
    xxx=Ne_mart1;
    SD(j)=rms(abs(xxx-xx));    
    
end 
%--计算误差等--%
ne_mart1=reshape(Ne_mart1,Nlat,Nlon,Nh);
toc
Re_Error1=mean(abs(Ne_mart1-Ne)./Ne);
ad1=mean(abs(Ne_mart1-Ne));
rmse1=sqrt(mean((Ne_mart1-Ne).^2));

neb=reshape(ne_background(11,5,:), [],1);
neb=smooth(neb);
nei=reshape(ne_ini(11,5,:), [],1);
nei=smooth(nei);
nem0=reshape(ne_mart0(11,5,:), [],1);
nem0=smooth(nem0);
nem1=reshape(ne_mart1(11,5,:), [],1);
nem1=smooth(nem1);

hh=h(1:end-1);
% Re_Error=mean(abs(Ne_mart(2000:12000)-Ne(2000:12000))./Ne(2000:12000));
figure
plot(neb,hh, nem0,hh, nem1,hh, nei,hh, 'linewidth', 1.5)
xlabel('ne/m^3','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Altitude/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
legend('Background','MART0','MART1','Initial')
grid on
savefig('./Figs12/12NewLine.fig')
%%保存计算结果

% save 'ne_mart_500.mat' 'ne_mart'
% save 'ne0.mat' 'ne_background' 'lat' 'lon' 'h'
% slice 画法
figure
xx=lat(1:Nlat);
yy=lon(1:Nlon);
zz=h(1:Nh);

[x,y,z]=meshgrid(xx,yy,zz);
% xslice=[xx(1),xx(9),xx(11),xx(13),xx(end)];
xslice=[xx(1),xx(10),xx(end)];
yslice=yy(end);
zslice=h(1);
slice(x,y,z,ne_background,xslice,yslice,zslice)
xlim([lat(1) lat(end)])
ylim([lon(1) lon(end)])
zlim([h(1) h(end)])
xlabel('Lat/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Lon/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
zlabel('Altitude/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
set(gca,'Ydir','Normal','fontsize',12,'fontname','Times New Roman','fontweight','bold')
title('12:00:00','fontsize',12,'fontname','Times New Roman','fontweight','bold')
shading interp
% colormap jet
colorbar
savefig('./Figs12/12NewNeb.fig')

% slice 画法
figure
% xslice=[xx(1),xx(9),xx(11),xx(13),xx(end)];
xslice=[xx(1),xx(10),xx(end)];
yslice=yy(end);
zslice=h(1);
slice(x,y,z,ne_mart0,xslice,yslice,zslice)
xlim([lat(1) lat(end)])
ylim([lon(1) lon(end)])
zlim([h(1) h(end)])
xlabel('Lat/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Lon/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
zlabel('Altitude/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
set(gca,'Ydir','Normal','fontsize',12,'fontname','Times New Roman','fontweight','bold')
title('MART0','fontsize',12,'fontname','Times New Roman','fontweight','bold')
shading interp
% title(['lamda=0.5 k=500 MART','交点覆盖率=',num2str(cover_rate)],'fontsize',12,'fontname','Times New Roman','fontweight','bold')
% colormap jet
colorbar
savefig('./Figs12/12NewNem0.fig')

figure
% xslice=[xx(1),xx(9),xx(11),xx(13),xx(end)];
xslice=[xx(1),xx(10),xx(end)];
yslice=yy(end);
zslice=h(1);
slice(x,y,z,ne_mart1,xslice,yslice,zslice)
xlim([lat(1) lat(end)])
ylim([lon(1) lon(end)])
zlim([h(1) h(end)])
xlabel('Lat/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('Lon/\circ','fontsize',12,'fontname','Times New Roman','fontweight','bold')
zlabel('Altitude/km','fontsize',12,'fontname','Times New Roman','fontweight','bold')
set(gca,'Ydir','Normal','fontsize',12,'fontname','Times New Roman','fontweight','bold')
title('MART1','fontsize',12,'fontname','Times New Roman','fontweight','bold')
shading interp
% title(['lamda=0.5 k=500 MART','交点覆盖率=',num2str(cover_rate)],'fontsize',12,'fontname','Times New Roman','fontweight','bold')
% colormap jet
colorbar
savefig('./Figs12/12NewNem1.fig')

figure
plot(Ne)
hold on
plot(Ne_mart0)
hold on
plot(Ne_mart1)
hold on
plot(ne00,'.')
hold on
% set(gca,'Ylim')
xlabel('n_{cell}','fontsize',12,'fontname','Times New Roman','fontweight','bold')
ylabel('ne_{mart}','fontsize',12,'fontname','Times New Roman','fontweight','bold')
title('MART','fontsize',12,'fontname','Times New Roman','fontweight','bold')
legend('Background','MART0','MART1','Initial')
savefig('./Figs12/New12Voxel.fig')
% title(['MART','交点覆盖率=',num2str(cover_rate)])
save ./Figs12/12NewNe ne_mart0 ne_mart1 ne_background ne_ini
save ./Figs12/12NewData delta ad0 ad1 ad2 rmse0 rmse1 rmse2 Re_Error0 Re_Error1 Re_Error2