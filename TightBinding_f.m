%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% last update 3Jan2022, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Ek]=TightBinding_f(Vpara,a,k,SpinOrbit)
%sp3s*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Semiconductor parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the "_" is for the "*"
Esa     = Vpara(1);
Epa     = Vpara(2);
Esc     = Vpara(3);
Epc     = Vpara(4);
Es_a    = Vpara(5);
Es_c    = Vpara(6);
Vss     = Vpara(7); Vs_s=0; % hopefully...
Vxx     = Vpara(8);
Vxy     = Vpara(9);
Vsa_pc  = Vpara(10);
Vsc_pa  = Vpara(11); Vpa_sc = Vsc_pa; % hopefully...
Vs_a_pc = Vpara(12);
Vpa_s_c = Vpara(13);

if SpinOrbit==1
  Da      = Vpara(14);
  Dc      = Vpara(15);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj=1:length(k(:,1))

kx = k(jj,1)*a/4;
ky = k(jj,2)*a/4;
kz = k(jj,3)*a/4;

g0=+cos(kx)*cos(ky)*cos(kz) -1i*sin(kx)*sin(ky)*sin(kz);
g1=-cos(kx)*sin(ky)*sin(kz) +1i*sin(kx)*cos(ky)*cos(kz);
g2=-sin(kx)*cos(ky)*sin(kz) +1i*cos(kx)*sin(ky)*cos(kz);
g3=-sin(kx)*sin(ky)*cos(kz) +1i*cos(kx)*cos(ky)*sin(kz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here is the Vogl notation
% J. Phys. Chem. Solids Vol 44, No. 5, pp365-378, 1983

%Hdiag = [ Esa Esc   Epa Epa Epa   Epc Epc Epc   Es_a Es_c ];
%
%% sa    sc      pa           pa           pa         pc         pc         pc         sa*       sc*
%
%H=[
%  0   Vss*g0     0            0            0       Vsa_pc*g1  Vsa_pc*g2  Vsa_pc*g3     0         0          % sa
%  0      0   -Vpa_sc*g1'  -Vpa_sc*g2'  -Vpa_sc*g3'    0          0          0          0         0          % sc
%  0      0       0            0            0       Vxx*g0     Vxy*g3     Vxy*g2        0     -Vpa_s_c*g1    % pa
%  0      0       0            0            0       Vxy*g3     Vxx*g0     Vxy*g1        0     -Vpa_s_c*g2    % pa
%  0      0       0            0            0       Vxy*g2     Vxy*g1     Vxx*g0        0     -Vpa_s_c*g3    % pa
%  0      0       0            0            0          0          0          0      Vs_a_pc*g1    0          % pc
%  0      0       0            0            0          0          0          0      Vs_a_pc*g2    0          % pc
%  0      0       0            0            0          0          0          0      Vs_a_pc*g3    0          % pc
%  0      0       0            0            0          0          0          0          0       Vs_s*g0      % sa*
%  0      0       0            0            0          0          0          0          0         0          % sc*
%];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here is the DiCarlo / Klimeck notation
% Superlattices and Microstructures, Vol. 27, No. 5/6, 2000
% Semicond. Sci. Technol. 18 (2003) R1–R31

Hdiag = [ Esa Epa Epa Epa Es_a Esc Epc Epc Epc Es_c ];

% sa      pa       pa       pa       sa*          sc          pc          pc           pc          sc*

H=[
  0       0        0        0        0         Vss*g0   Vsa_pc*g1   Vsa_pc*g2    Vsa_pc*g3         0          % sa
  0       0        0        0        0     -Vpa_sc*g1      Vxx*g0      Vxy*g3       Vxy*g2     -Vpa_s_c*g1    % pa
  0       0        0        0        0     -Vpa_sc*g2      Vxy*g3      Vxx*g0       Vxy*g1     -Vpa_s_c*g2    % pa
  0       0        0        0        0     -Vpa_sc*g3      Vxy*g2      Vxy*g1       Vxx*g0     -Vpa_s_c*g3    % pa
  0       0        0        0        0            0    Vs_a_pc*g1  Vs_a_pc*g2   Vs_a_pc*g3         0          % sa*
  0       0        0        0        0            0           0           0            0           0          % sc
  0       0        0        0        0            0           0           0            0           0          % pc
  0       0        0        0        0            0           0           0            0           0          % pc
  0       0        0        0        0            0           0           0            0           0          % pc
  0       0        0        0        0            0           0           0            0           0          % sc*
];



H=H'+H+diag(Hdiag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SpinOrbit==1

H0=zeros(10,10);
H=[H H0 ; H0 H'];

lambdaA=Da/3;
lambdaC=Dc/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hso11=[ 0     0     0     0     0
        0     0    -1i    0     0
        0    1i     0     0     0
        0     0     0     0     0
        0     0     0     0     0
];

Hso12=[ 0     0     0     0     0
        0     0     0     1     0
        0     0     0    -1i    0
        0    -1    1i     0     0
        0     0     0     0     0
];

Ha11=lambdaA*Hso11;       %SO Hamiltonian (up,up) or (1,1) -anion
Hc11=lambdaC*Hso11;       %SO Hamiltonian (up,up) or (1,1) -cation

Ha12=lambdaA*Hso12;       %SO Hamiltonian (up,down) or (1,2) -anion
Hc12=lambdaC*Hso12;       %SO Hamiltonian (up,down) or (1,2) -cation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0=zeros(5,5);

HSO=[ Ha11    H0      Ha12        H0
       H0    Hc11      H0        Hc12
      Ha12'   H0    conj(Ha11)    H0
       H0    Hc12'     H0      conj(Hc11) 
];


H=H+HSO;

end


Ek(:,jj) = eig(H);

end

end


