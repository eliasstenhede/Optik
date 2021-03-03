% 
function [x,y,M]=xy_source(N,D_star,separation)

% Utdata
% x: matris med punktk�llornas x-positioner [m] 
% y: matris med punktk�llornas y-positioner [m]
% M: punktk�llornas antal

% Indata
% N: Antalet observationspunkter (punkter l�ngs u-axeln d�r f�ltet ber�knas)
% D_star: Stj�rnans diameter [m]
% separation: avst�nd mellan punktk�llorna p� stj�rnan [m],
%             kan l�mpligen anges som br�kdel av D_star, 
%             t.ex. separation=D_star/30


xpos=-D_star/2:separation:D_star/2; 
ypos=xpos;
[xposmat,yposmat]=meshgrid(xpos,ypos); 
rposmat=sqrt(xposmat.^2+yposmat.^2);

indexvekt=find(rposmat<(D_star/2));
M=length(indexvekt);

x_vekt_source=xposmat(indexvekt);
y_vekt_source=yposmat(indexvekt);

alfa_plot=linspace(0,2*pi,100);
x_circumference_plot=D_star/2*cos(alfa_plot);
y_circumference_plot=D_star/2*sin(alfa_plot);

%figure(1000)
%close(1000)
%figure(1000)
%plot(x_vekt_source/1e6/1e3,y_vekt_source/1e6/1e3,'r.')
%hold on
%plot(x_circumference_plot/1e6/1e3,y_circumference_plot/1e6/1e3,'k')
%xlabel('x [miljoner km]')
%ylabel('y [miljoner km]')
%title(['K�llpositioner p� stj�rnan. Antal k�llor M=' num2str(M) ])
%axis('equal')
%hold off

x=repmat(x_vekt_source,1,N);
y=repmat(y_vekt_source,1,N);


