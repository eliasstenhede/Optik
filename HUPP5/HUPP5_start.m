close all
clear

N=256;
sidlaengd_Plan1=100e-6;
%sidlaengd_Plan1=30e-6;
omega_in=6e-6;
lambda_noll=1550e-9;
k_noll=2*pi/lambda_noll;
a=sidlaengd_Plan1/N;

xvekt=-N/2*a:a:(N/2-1)*a;
yvekt=xvekt;
[xmat,ymat]=meshgrid(xvekt,yvekt);
rmat=sqrt(xmat.^2+ymat.^2);

% brytningsindexvariation i (x,y)-led
n_core=1.51; % i k�rnan
n_clad=1.50; % i h�ljet
D_core=60e-6;
%nmat=(rmat<=D_core/2)*n_core+(rmat>D_core/2)*n_clad;
nmat = nmat_GRIN(n_core,n_clad,D_core,xmat,ymat);
figure(1)
imagesc(xvekt*1e6,yvekt*1e6,nmat)
xlabel('x [�m]')
ylabel('y [�m]')
colormap(jet)
colorbar
title('Brytningsindexvariation')
drawnow


% d�mpmatris (beh�ver normalt inte �ndras)
r_daemp_start=0.8*N/2*a; % ut till denna radie sker ingen d�mpning
kantvaerde=0.8; % v�rdet p� daempmat vid kanten av ber�kningsf�nstret (l�ngs xy-axlarna)
daempmat=(rmat<=r_daemp_start)*1+(rmat>r_daemp_start).*(1-(1-kantvaerde)/(N/2*a-r_daemp_start)^2.*(rmat-r_daemp_start).^2);
deampmat(rmat<D_core/2)=1;
daempmat(daempmat<0)=0;
figure(2)
mesh(xvekt*1e6,yvekt*1e6,daempmat)
xlabel('x [�m]')
ylabel('y [�m]')
title('Daempmat')
colormap(jet)
drawnow
disp('Tryck valfri tangent f�r att forts�tta!')
pause

% startf�lt "l�ngst till v�nster"
E_start=exp(-(xmat.^2+ymat.^2)/omega_in^2);
E_start_offset = exp(-((xmat-10e-6).^2+(ymat).^2)/omega_in^2);
E_start_ny_mod = E_start.*ymat;
% total propagationsstr�cka och BPM-steg-storlek
L=1000e-6*12;
delta_z=lambda_noll;
Lvekt=delta_z:delta_z:L;

%E1=(E_start_ny_mod + max(max(E_start_ny_mod))*E_start/(max(max(E_start))));
E1=E_start_offset.*exp(-1j*k_noll*sin(0.0698).*ymat);
E_sida=zeros(N,length(Lvekt));
I_sida_norm=zeros(N,length(Lvekt));
steg_nummer=0;
for akt_L=Lvekt
    steg_nummer=steg_nummer+1;
    
    E2=BPM_steg(E1,delta_z,N,a,lambda_noll,nmat,daempmat); % F�ltf�rdelning i Plan 2
    
    I2=abs(E2).^2; % Intensitetsf�rdelning i Plan 2
    
    % Vektorer/matriser f�r plottning sett "fr�n sidan"
    E2_laengs_yaxeln=E2(:,N/2+1);
    E_sida(:,steg_nummer)=E2_laengs_yaxeln;
    I2_laengs_yaxeln=abs(E2_laengs_yaxeln).^2;
    I2_laengs_yaxeln_norm=I2_laengs_yaxeln/max(I2_laengs_yaxeln);
    I_sida_norm(:,steg_nummer)=I2_laengs_yaxeln_norm;
    
    if rem(steg_nummer,500)==0 % f�r att snabba p� simuleringen plottas inte alla BPM-steg
        
        figure(10)
        imagesc(xvekt*1e6,yvekt*1e6,I2)
        hold on
        plot(D_core/2*cos(linspace(0,2*pi,50))*1e6,D_core/2*sin(linspace(0,2*pi,50))*1e6,'Color',[1 1 1]*0.6,'LineWidth',2)
        title(['Intensitet i tv�rsnitt efter ' num2str(akt_L*1e3) ' mm propagation' ])
        axis('square')
        xlabel('x [�m]')
        ylabel('y [�m]')
        colormap(jet)
        hold off
        drawnow
        %pause
        
        
        figure(11)
        imagesc(Lvekt*1e3,yvekt*1e6,I_sida_norm)
        title(['I_-sida: Normerad intensitet efter ' num2str(akt_L*1e3) ' mm propagation' ])
        xlabel('z [mm]')
        ylabel('y [�m]')
        colormap(jet)
        drawnow
        %pause
        
    end
    
    E1=E2;
end
disp("Origo: " + I2( round(length(I2)/2), round(length(I2)/2) ))

function E2=BPM_steg(E1,delta_z,N,a,lambda_noll,nmat,daempmat)
    k_noll=2*pi/lambda_noll;
    n_PAS=mean(mean(nmat));
    E2_PAS=PAS(E1,delta_z,N,a,lambda_noll,n_PAS);
    faskorrektion=k_noll*delta_z*(nmat-n_PAS);
    if max(max(abs(faskorrektion)))>(2*pi*0.02) % max till�ten faskorrektion �r 2% av 2pi
         disp('Stegl�ngd f�r stor')
         E2=0;
    else
        E2=E2_PAS.*exp(1i*faskorrektion).*daempmat;
    end
end
