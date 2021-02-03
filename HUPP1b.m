clear; close all; full_white_value=255;

T_DOE = load("T_DOE_gen2").T_DOE_gen2;
N=1024; % NxN är matrisstorleken
sidlaengd_Plan1=4e-3;
a=sidlaengd_Plan1/N; % samplingsavstånd i Plan 1 (och Plan 2 eftersom vi använder PAS)
L=20e-3; % propagationssträcka (dvs avstånd mellan Plan 1 och 2)

lambda_noll=633e-9; % vakuumvåglängd
n_medium=1; % brytningsindex för medium mellan Plan 1 och 2
k=2*pi*n_medium/lambda_noll; 

xvekt=-N/2*a:a:(N/2-1)*a; % vektor med sampelpositioner i x-led
yvekt=xvekt;
[xmat,ymat]=meshgrid(xvekt,yvekt);
rmat=sqrt(xmat.^2+ymat.^2);

f_dvl=-0.15; % fokallängd på de vises lins före Plan 1
T_dvl=exp(-1i*k*rmat.^2/(2*f_dvl)); % Transmissionsfunktion för en lins (linsen är TOK)
f_eye=20e-3; % fokallängd på linsen före Plan 1
T_eye=exp(-1i*k*rmat.^2/(2*f_eye)); % Transmissionsfunktion för en lins (linsen är TOK)

E1 = ones(N,N).*T_eye.*T_dvl.*T_DOE;

I1=abs(E1).^2;

E2=PAS(E1,L,N,a,lambda_noll,n_medium);

I2=abs(E2).^2;

mattnadsfaktor_plot=max(50);
figure(3)
image(xvekt*1e3,yvekt*1e3,I2/max(max(I2))*full_white_value*mattnadsfaktor_plot)
title(['Intensitet efter ' num2str(L*1e3) ' mm propagation (mattnadsfaktor=' num2str(mattnadsfaktor_plot) '). Verkar OK, eller?'])
xlabel('x[mm]')
ylabel('y[mm]')
colormap(gray) 
drawnow
axis('equal')
