function E2=PAS(E1,L,N,a,lambda_noll,n_medium)

% Varje sampelpunkt i k-planet motsvarar en plan våg med en viss riktning (kx,ky,kz)
delta_k=2*pi/(N*a); % samplingsavstånd i k-planet
kxvekt=-N/2*delta_k:delta_k:(N/2-1)*delta_k; % vektor med sampelpositioner i kx-led
kyvekt=kxvekt;
[kxmat,kymat]=meshgrid(kxvekt,kyvekt);

k=2*pi*n_medium/lambda_noll; 
kzmat=sqrt(k^2-kxmat.^2-kymat.^2);
fasfaktor_propagation=exp(1i.*kzmat*L); 
A= a^2/(2*pi)^2*fft2c(E1);
B=A.*fasfaktor_propagation;
E2= delta_k^2*N^2*ifft2c(B);
