function[k_ZB]=kZB_f(Nk,a)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Brillouin zone vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% L-valley
kL=linspace(-0.5,0,Nk);
k1 = [ kL'  kL'  kL' ];

%%% X-valley
kx=linspace(0,1,Nk);
k2 = [ kx'  0*kx'  0*kx' ];

%%% X point to U point
%ky=linspace(0,0.25,floor(Nk/2));
%kz=linspace(0,0.25,floor(Nk/2));
%k3 = [ ky'*0+1  ky'  kz' ];

ky=linspace(0,0.125,floor(Nk/2));
kz=linspace(0,0.125,floor(Nk/2));
k3 = [ ky'*0+1  ky'  kz' ];

%%% U point to Gamma point
kx=linspace(0.875,0,Nk);
ky=linspace(0.875,0,Nk);
k4 = [ kx'  ky'  kx'*0 ];

%kx=linspace(1,0,Nk);
%ky=linspace(0.125,0,Nk);
%k4 = [ kx'  ky'  ky' ];


%k_ZB=[ k1 ; k2  ]*2*pi/a;
k_ZB=[ k1 ; k2 ; k3 ; k4 ]*2*pi/a;

