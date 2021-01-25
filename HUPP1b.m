clear
close all
full_white_value=64; % äldre matlabversion - detta värde plottas som vitt (max) med image-kommandot
%full_white_value=255; % nyare matlabversion - prova denna om din plot verkar mörk!

T_DOE = load("T_DOE_gen2").T_DOE_gen2;
N=1024; % NxN är matrisstorleken (rekommenderad storlek N=1024)
sidlaengd_Plan1=4e-3; % det samplade områdets storlek (i x- eller y-led) i Plan 1 (rekommenderad storlek 4 mm)
a=sidlaengd_Plan1/N; % samplingsavstånd i Plan 1 (och Plan 2 eftersom vi använder PAS)
L=20e-3; % propagationssträcka (dvs avstånd mellan Plan 1 och 2)

lambda_noll=633e-9; % vakuumvåglängd för rött ljus från en HeNe-laser
n_medium=1; % brytningsindex för medium mellan Plan 1 och 2
k=2*pi*n_medium/lambda_noll; 

xvekt=-N/2*a:a:(N/2-1)*a; % vektor med sampelpositioner i x-led
yvekt=xvekt; % och y-led
[xmat,ymat]=meshgrid(xvekt,yvekt); % koordinatmatriser med x- och y-värdet i varje sampelposition
rmat=sqrt(xmat.^2+ymat.^2); % avståndet till origo i varje sampelpunkt. Observera att alla operationer är elementvisa!

%******* Fält i Plan 1
f_dvl=20e-3; % fokallängd på linsen före Plan 1
T_dvl=exp(-1i*k*rmat.^2/(2*f_dvl)); % Transmissionsfunktion för en lins (linsen är TOK)

f_eye=20e-3; % fokallängd på linsen före Plan 1
T_eye=exp(-1i*k*rmat.^2/(2*f_eye)); % Transmissionsfunktion för en lins (linsen är TOK)

E1 = ones(N,N).*T_eye.*T_dvl.*T_DOE;

I1=abs(E1).^2;

E2=PAS(E1,L,N,a,lambda_noll,n_medium);

I2=abs(E2).^2;

mattnadsfaktor_plot=50; % anger hur många gånger maxvärdet ska vara mättat i plotten (>1, kan vara bra om man vill se svagare detaljer)
figure(3)
image(xvekt*1e3,yvekt*1e3,I2/max(max(I2))*full_white_value*mattnadsfaktor_plot)
title(['Intensitet efter ' num2str(L*1e3) ' mm propagation (mattnadsfaktor=' num2str(mattnadsfaktor_plot) '). Verkar OK, eller?'])
xlabel('x[mm]')
ylabel('y[mm]')
colormap(gray) 
drawnow
axis('equal')
