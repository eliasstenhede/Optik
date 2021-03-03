n=10000;
instantan_summa = 0; I_obs_summa = 0; Gamma_I_summa = 0;
N = 50; %Värden över 25 verkar ge okej resultat.
D_star = 13927*70e5;
separation = D_star / N;
k_noll = 2*pi/650e-9;
L_star  = 6.62237e17;
[x,y,M]=xy_source(N,D_star,separation);

for i = 1:n
    fas_vekt=rand(M,1)*2*pi;
    fas=repmat(fas_vekt,1,N);

    uvekt = linspace(0, 20, length(x(1,:)));
    u=repmat(uvekt,M,1);

    r = u.*x/L_star;

    E_k_obs = exp(1i*(fas+k_noll*r));
    E_obs=sum(E_k_obs,1);

    instantan_produkt=E_obs(1)*conj(E_obs);
    instantan_summa = instantan_summa + instantan_produkt;
    I_obs_inst=abs(E_obs).^2;
    I_obs_summa = I_obs_summa + I_obs_inst;
    Gamma_I = I_obs_inst(1)*I_obs_inst ; %Uppgift f
    Gamma_I_summa = Gamma_I_summa + Gamma_I; %Uppgift f
    
    if mod(i, 1000) == 0
        disp("steg: " + i)
        
        figure(1000)
        plot(uvekt, abs(instantan_summa/max(instantan_summa)))
        xlabel('u [m]')
        ylabel('Koherens')
        title(["Normerad \Gamma_{AB}, N = " + num2str(N)])
        
        figure(2000)
        plot(uvekt, I_obs_summa/max(I_obs_summa))
        xlabel('u [m]')
        ylabel('Intensitet')
        title(["Normerad intensitet, N = " + num2str(N)])
        ylim([0,1])
        
        figure(3000) %Uppgift f
        plot(uvekt, Gamma_I_summa/max(Gamma_I_summa))
        xlabel('u [m]')
        ylabel('Koherens')
        title(["Normerad \Gamma_{I}, N = " + num2str(N)])
    end
end

