function[Ek]=epm_f(Vpseudo,a,k,Nband)

% https://nanohub.org/resources/4883/download/tutorial_semiempirical_bandstructure_methods.pdf
% https://web.northeastern.edu/afeiguin/phys5870/phys5870/node47.html
% https://www.ece.nus.edu.sg/stfpage/eleadj/pseudopotential.htm
% Yu and Cardona, "Fundamentals of Semiconductors", p58, 2.5/Band Structure Calculations by Pseudopotential Methods
% Cohen and Bergstresser, Phys. Rev, 141,p789, (1966)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NG=137; %% number of plane waves (size of the matrix to diagonalize max 137 from the epm_Gvec.txt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

me=9.10938188E-31;              %% electron mass kg
e=1.602176487E-19;              %% charge de l electron Coulomb
h=6.62606896E-34;               %% Planck constant J.s
hbar=h/(2*pi);
Ry=13.6056925330;               %% Rydberg energy in eV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Semiconductor parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vs3 =Vpseudo(1);
Vs8 =Vpseudo(2);
Vs11=Vpseudo(3);
Va3 =Vpseudo(4);
Va4 =Vpseudo(5);
Va11=Vpseudo(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Building Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gvec=load('epm_Gvec.txt');       %% the G vector are save in a text file


for i=1:3
    Gveci = repmat( Gvec(1:NG,i)  , 1  , NG );
    Gvecj = Gveci' ;
    K(:,:,i) = Gveci - Gvecj;
end



KK=sum(K.*K,3);
idx3  = find(KK==3);
idx4  = find(KK==4);
idx8  = find(KK==8);
idx11 = find(KK==11);


TT=zeros(NG,NG,3)+a*1/8;

Vs=zeros(NG,NG);
Vs(idx3)  = Vs3  * Ry ;
Vs(idx8)  = Vs8  * Ry ;
Vs(idx11) = Vs11 * Ry ;

Va=zeros(NG,NG);
Va(idx3)  = Va3  * Ry ;
Va(idx4)  = Va4  * Ry ;
Va(idx11) = Va11 * Ry ;

IDX=1-diag(ones(1,NG));
HV=Vs .* cos(2*pi/a*  sum(K.*TT,3)   ) + 1j*Va .* sin(2*pi/a*  sum(K.*TT,3)   ) ;
%HV= HV .* IDX ;

Gvec=Gvec*2*pi/a;


for i=1:length(k(:,1))
    
    kk=repmat(k(i,:),NG,1);
    Gk = hbar^2/(2*me*e) * sum( (Gvec(1:NG,:) + kk) .* (Gvec(1:NG,:) + kk) , 2 ) ;
    Gkdiag=diag(Gk);
    
    H = Gkdiag + HV ;
    % Take care... 
	% (sparse+eig) does not work on Matlab BUT does work with Octave
	% (eig) works on Matlab AND Octave
	% (sparse+eigs) makes weird band diagram
    % H=sparse(H);
    [Vk,Emat]=eig(H);
    E(:,i)=diag(Emat);  
    
    if k(i,:)==[0  0  0]
        V=Vk;
    end

end;

Ek=E(1:Nband,:);        %% cut the solution to have just the first ones

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
