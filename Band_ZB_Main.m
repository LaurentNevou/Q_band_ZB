%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% last update 8Nov2021, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Here, you have to choose your material among the following %%%%%%%%%% 
%%%%%% Take care, not all the models are available for all the materials! %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material='AlAs';
Material='GaAs';
%Material='InAs';
%Material='AlSb';
%Material='GaSb';
%Material='InSb';
%Material='AlP';
%Material='GaP';
%Material='InP';
%Material='C';
%Material='Si';
%Material='Ge';
%Material='Sn';
%Material='3C-SiC';
%Material='ZnSe';
%Material='ZnTe';
%Material='CdS';
%Material='ZnS';
%Material='CdSe';
%Material='CdTe';
%Material='Void';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=300;                  % Temperature in Kelvin, Has an impact only on k.p model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library
ExtractParameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Model activation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on

epm_model = 0;               %% Empirical Pseudo Potential
TB_model  = 1; SpinOrbit=1;  %% Tight Binding; SpinOrbit coupling = 0 or 1
kp8_model = 1;               %% kp 8 bands model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nk=50;          %% number of k-points for the dispersion
Nband=12;        %% number of band to plot

k=kZB_f(Nk,a);   %% function to compute the k-vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if epm_model==1
  if strcmp(M{1},'CdS')==1 || strcmp(M{1},'CdSe')==1
    display(strcat('Error: No parameter in the Pseudo Potential model for: ',M{1}))
    return
  else
    Ek1=epm_f(M{3},a,k,Nband);
  end
end

if TB_model==1
  if SpinOrbit==0
    if strcmp(M{1},'CdTe')==1 || strcmp(M{1},'3C-SiC')==1
      display(strcat('Error: No parameter in the Tight Binding (without Spin) model for: ',M{1}))
      return
    else 
      Ek2=TightBinding_f(M{4},a,k,SpinOrbit);
    end
  elseif SpinOrbit==1
    if strcmp(M{1},'C')==1 || strcmp(M{1},'Ge')==1 || strcmp(M{1},'Sn')==1 || strcmp(M{1},'3C-SiC')==1 || strcmp(M{1},'ZnTe')==1
      display(strcat('Error: No parameter in the Tight Binding (with Spin) model for: ',M{1}))
      return
    else 
      Ek2=TightBinding_f(M{5},a,k,SpinOrbit);
    end
  end
end

if kp8_model==1
  if strcmp(M{1},'C')==1 || strcmp(M{1},'Sn')==1 || strcmp(M{1},'3C-SiC')==1
    display(strcat('Error: No parameter in the k.p 8bands model for: ',M{1}))
    return
  else
    Ek3=kp_8bands_DKK_f(k, Eg, EP_L, Dso, F, g1, g2, g3);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X0fig=-1800; Y0fig=100;
X0fig=100; Y0fig=100;
Wfig=800;Hfig=800;

figure('Name','Results','position',[X0fig Y0fig Wfig Hfig],'color','w')

FS=18;
LW=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xscale=[0 2*Nk];
yscale=[-2 4];
%xscale=[0 3.5*Nk];
%yscale=[-13 7];

subplot(1,1,1,'fontsize',FS)
hold on;grid on;box on;
s={};
c=0;
if epm_model==1
  c=c+1;
  plot(1:length(k(:,1)),Ek1-max(Ek1(4,:)),'b-', 'linewidth',LW)
  s{c}=strcat('\fontsize{',num2str(FS),'}\color{blue}PseudoPotential');
end

if TB_model==1
  c=c+1;
  plot(1:length(k(:,1)),Ek2,'g-', 'linewidth',LW)
  s{c}=strcat('\fontsize{',num2str(FS),'}\color{green}TightBinding sp3s*');
end

if kp8_model==1
  c=c+1;
  plot(1:length(k(:,1)),Ek3,'r-', 'linewidth',LW)
  s{c}=strcat('\fontsize{',num2str(FS),'}\color{red}k.p 8bands');
end

title(strcat(M{1},' bandstructure'));
xlabel('Wavevector k');
ylabel('Energy (eV)');

plot( 0*[Nk Nk] , yscale , 'k' )
plot( [Nk Nk] , yscale , 'k' )
plot( 2*[Nk Nk] , yscale , 'k' )
plot( 2.5*[Nk Nk] , yscale , 'k' )
plot( 3.5*[Nk Nk] , yscale , 'k' )

text(xscale(2)*0.75,yscale(2)*1.0,s,'edgecolor','k','backgroundcolor','w')
set(gca,'xtick',[0 Nk 2*Nk 2.5*Nk 3.5*Nk ],'ytick',yscale(1):1:yscale(2), 'xticklabel',{'L' '\Gamma' 'X' 'K' '\Gamma'})

xlim(xscale)
ylim(yscale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%