
%%%%%%%%%% Hall ladder
% 能带测量和布居数分布

clear;
clc;
close all

t1=xlsread('tek0010.csv','A22:A100021');
y1=xlsread('tek0010.csv','C22:C100021');
y2=xlsread('tek0010.csv','C22:C100021');


figure
plot(t1,y1,'k','LineWidth',0.2)
xlabel('t');
ylabel('I');
set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold','linewidth',1.5);
axis tight


figure
plot(t1,y2,'k','LineWidth',0.2)
xlabel('t');
ylabel('I');
set(gca,'FontSize',16,'FontName','Times New Roman','FontWeight','bold','linewidth',1.5);
axis tight


C=2*142;
s=length(t1)/C;

for i=1:s
A(i,:)=y2(1+(i-1)*C:1+(i-1)*C+C-1);

end

[T,Ome]=meshgrid(1:C,1:s);

AA=-A-min(min(-A));
Au=max(max(AA))-0.1*max(max(AA));
T_out=AA./Au;

figure
surf(T,Ome,T_out)
shading interp%装饰，变平滑
colorbar
xlabel('t','FontName','Times New Roman'); 
ylabel('\omega','FontName','Times New Roman'); 
zlabel('T_{out}','FontName','Times New Roman'); 
set(gca,'FontSize',15','FontName','Times New Roman');
view(0,90)
axis tight
%print('-dpng','-r300',['C:\Users\acer\Desktop\',num2str(1)])

%%
% 归一化导频率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f=2./s;
oavg=152;
omin=oavg-round(0.5./f);
omax=oavg+round(0.5./f);
D=1;
oo=[omin:D:omax];

ft=4./C;
tovg=108;
tmin=tovg-round(1./ft);
tmax=tovg+round(1./ft);
tt=[tmin:tmax];

CC=1:C;
SS=1:s;

[T1,Ome1]=meshgrid(CC(tt)*ft-CC(tovg)*ft,SS(oo)*f-SS(oavg)*f);

figure
surf(T1,Ome1,T_out(oo,tt))
shading interp%装饰，变平滑
colorbar
xlabel('t','FontName','Times New Roman'); 
ylabel('\omega','FontName','Times New Roman'); 
zlabel('T_{out}','FontName','Times New Roman'); 
set(gca,'FontSize',15','FontName','Times New Roman');
view(0,90)
axis tight

axis off
% print('-dpng','-r300',['C:\Users\acer\Desktop\',num2str(1)])


%%
% 拟合实部虚部%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 双峰拟合

K=1:2:C./2;
T_out2=T_out(oo,tt);

xi_pu=SS(oo)*f-SS(oavg)*f;

Ren_min=-0.4; Ren_max=0.1; Ren_sta=-0.25;
Rep_min=-0.1; Rep_max=0.3; Rep_sta=0.25;
Imp_min=0.0002; Imp_max=0.2; Imp_sta=0.1;
Imn_min=0.0002; Imn_max=0.2; Imn_sta=0.1;
Rp_min=0.0002; Rp_max=1; Rp_sta=0.2;
Rn_min=0.0002; Rn_max=1; Rn_sta=0.2;

ff=fittype('-imag(Rp/(Ep_RE-delta+i*Ep_Im)+(Rn)/(En_RE-delta+i*En_Im))','independent','delta',...
    'coefficients',{'Ep_RE','Ep_Im','En_RE','En_Im','Rp','Rn'});

% 
% ff=fittype('abs((Rp/(Ep_RE-delta+i*Ep_Im)+(Rn)/(En_RE-delta+i*En_Im)))^2','independent','delta',...
%     'coefficients',{'Ep_RE','Ep_Im','En_RE','En_Im','Rp','Rn'});


options = fitoptions(ff);
options.StartPoint =[Rep_sta Imp_sta Ren_sta Imn_sta Rp_sta Rn_sta];
options.Lower = [Rep_min Imp_min Ren_min Imn_min Rp_min Rn_min];
options.upper = [Rep_max Imp_max Ren_max Imn_max Rp_max Rn_max];


for k=1:length(K)
    kk=K(k);
    T_out3=T_out2(:,kk)./(max(max(T_out2(:,kk))));
   
% enhanced_T = peak_enhancement_filter(T_out3, 3 , 3);
% T_out3=enhanced_T;
    cfun_pu=fit(xi_pu',T_out3,ff,options);

Ep_RE=cfun_pu.Ep_RE;
Ep_Im=cfun_pu.Ep_Im;
En_RE=cfun_pu.En_RE;
En_Im=cfun_pu.En_Im;

Rp=cfun_pu.Rp;
Rn=cfun_pu.Rn;

Y=cfun_pu(xi_pu);

figure
plot(xi_pu,T_out3,'ko','LineWidth',2)
hold on
plot(xi_pu,Y,'b','Linewidth',2)
xlabel('\delta\omega/\Omega','FontName','Times New Roman'); 
ylabel('T_{out}','FontName','Times New Roman'); 
title(['t=',num2str(kk),' k=',num2str(k)])
set(gca,'FontSize',15','FontName','Times New Roman');
axis tight
%axis([-50 -49 -2.4 -2.35 0 inf])


EP_RE(k)=cfun_pu.Ep_RE;
EP_IM(k)=cfun_pu.Ep_Im;
EN_RE(k)=cfun_pu.En_RE;
EN_IM(k)=cfun_pu.En_Im;

k
end
%

figure

subplot(221)
plot(K./C*4,EP_RE,'ro','Linewidth',2)
hold on
plot(K./C*4,EN_RE,'bo','Linewidth',2)
xlabel('\it{k_f}\rm{\Omega}/\pi','FontName','Times New Roman'); 
ylabel('ReE/\Omega','FontName','Times New Roman'); 
set(gca,'FontSize',15','FontName','Times New Roman');
axis([0 2 min(xi_pu) max(xi_pu)])
% axis tight



subplot(222)
plot(K./C*4,-abs(EP_IM),'r*','Linewidth',2)
hold on
plot(K./C*4,-abs(EN_IM),'b*','Linewidth',2)
xlabel('\it{k_f}\rm{\Omega}/\pi','FontName','Times New Roman'); 
ylabel('ImE/\Omega','FontName','Times New Roman'); 
set(gca,'FontSize',15','FontName','Times New Roman');
axis([0 2 0.02 0.12])
axis tight



subplot(223)
plot(EP_RE,-abs(EP_IM),'rd','Linewidth',2)
hold on
plot(EN_RE,-abs(EN_IM),'bd','Linewidth',2)
xlabel('ReE/\Omega','FontName','Times New Roman'); 
ylabel('ImE/\Omega','FontName','Times New Roman'); 
set(gca,'FontSize',15','FontName','Times New Roman');
axis([min(xi_pu) max(xi_pu) -0.2 -0.02])
axis tight

subplot(224)
scatter(EP_RE,-abs(EP_IM),20,K./C*4,'filled')
hold on
scatter(EN_RE,-abs(EN_IM),20,K./C*4,'filled')
colormap(gca,'jet')
colorbar
%grid on
box on
grid on
%axis equal
xlabel('Re(E)/\Omega'); 
ylabel('Im(E)/\Omega'); 
set(gca,'FontSize',18,'FontName','Times New Roman','linewidth',1.5);
axis([-0.3 0.3 -0.2 0])

set(gcf,'unit','normalized','position',[0.1,0.3,0.8,0.5]);

figure
plot(xi_pu,T_out3,'.')
