clc, clf
d = 10e-6;
n_o = 1.5;
N_skiv = 100;
lambda_r = 650e-9;
lambda_g = 500e-9;
lambda_b = 450e-9;
thetas = linspace(0, pi/2, 100);
E_0 = [1;0];

hold on;
plot(thetas, g(thetas, N_skiv, d, n_o, lambda_r, E_0), 'r')
plot(thetas, g(thetas, N_skiv, d, n_o, lambda_g, E_0), 'g')
plot(thetas, g(thetas, N_skiv, d, n_o, lambda_b, E_0), 'b')

function y = g(thetas, N_skiv, d, n_o, lambda, E_0)
    for i = 1:100
        y(i) = norm(prop(thetas(i), N_skiv, d, n_o, lambda, E_0));
    end
end

function E = prop(theta, N_skiv, d, n_o, lambda, E_0)
    E = E_0;
    delta = d/N_skiv;
    for n = 1:N_skiv
        alfa_n = pi*n/(2*N_skiv);
        E = J_ret(alfa_n, delta, f(theta, n_o), n_o, lambda)*E;
    end
    E = J_pol(pi/2)*E;
end

function matris = J_proj(alfa)
    matris = [cos(alfa), sin(alfa); -sin(alfa), cos(alfa)]; %eller andra hållet?
end

function matris = J_pol(alfa)
    matris = J_proj(-alfa)*[1, 0; 0, 0]*J_proj(alfa);
end

function matris = J_ret(alfa_n, delta, n_eo, n_o, lambda)
    phi = 2*pi/lambda*(n_eo-n_o)*delta;
    matris = J_proj(-alfa_n)*[exp(1i*phi), 0; 0, 1]*J_proj(alfa_n);
end

function out = f(theta, n_o)
    n_eo = 1.6;
    out  = n_eo*n_o/norm([n_o*cos(theta), n_eo*sin(theta)]);
end