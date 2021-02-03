clear; close all; full_white_value=255;

T_DOE = load("T_DOE_gen2").T_DOE_gen2;
N=1024; % NxN �r matrisstorleken
sidlaengd_Plan1=4e-3;
a=sidlaengd_Plan1/N; % samplingsavst�nd i Plan 1 (och Plan 2 eftersom vi anv�nder PAS)
L=20e-3; % propagationsstr�cka (dvs avst�nd mellan Plan 1 och 2)

lambda_noll=633e-9; % vakuumv�gl�ngd
n_medium=1; % brytningsindex f�r medium mellan Plan 1 och 2
k=2*pi*n_medium/lambda_noll; 

xvekt=-N/2*a:a:(N/2-1)*a; % vektor med sampelpositioner i x-led
yvekt=xvekt;
[xmat,ymat]=meshgrid(xvekt,yvekt);
rmat=sqrt(xmat.^2+ymat.^2);

f_dvl=-0.15; % fokall�ngd p� de vises lins f�re Plan 1
T_dvl=exp(-1i*k*rmat.^2/(2*f_dvl)); % Transmissionsfunktion f�r en lins (linsen �r TOK)
f_eye=20e-3; % fokall�ngd p� linsen f�re Plan 1
T_eye=exp(-1i*k*rmat.^2/(2*f_eye)); % Transmissionsfunktion f�r en lins (linsen �r TOK)

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
