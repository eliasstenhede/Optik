%G�ger glas �r polariserad i x-led och v�nster glas i y-led.

%c)

E_in = [1;1];

E_ut_i = J_pol(0)*J_ret(-pi/4, pi/2)*J_ret(pi/4, pi/2)*J_pol(0)*E_in;
E_ut_ii = J_pol(0)*J_ret(pi/4, pi/2)*J_ret(pi/4, pi/2)*J_pol(0)*E_in;
E_ut_iii = J_pol(0)*J_ret(-pi/4, pi/2)*J_ret(-pi/4, pi/2)*J_pol(0)*E_in;
E_ut_iv = J_pol(0)*J_ret(pi/4, pi/2)*J_ret(-pi/4, pi/2)*J_pol(0)*E_in;

I_ut_i = norm(E_ut_i)^2;
I_ut_ii = norm(E_ut_ii)^2;
I_ut_iii = norm(E_ut_iii)^2;
I_ut_iv = norm(E_ut_iv)^2;

%d)
E_ut_i = J_pol(0)*J_ret_icke_ideal(-pi/4, pi/2)*J_ret_icke_ideal(pi/4, pi/2)*J_pol(0)*E_in;
E_ut_ii = J_pol(0)*J_ret_icke_ideal(pi/4, pi/2)*J_ret_icke_ideal(pi/4, pi/2)*J_pol(0)*E_in;
E_ut_iii = J_pol(0)*J_ret_icke_ideal(-pi/4, pi/2)*J_ret_icke_ideal(-pi/4, pi/2)*J_pol(0)*E_in;
E_ut_iv = J_pol(0)*J_ret_icke_ideal(pi/4, pi/2)*J_ret_icke_ideal(-pi/4, pi/2)*J_pol(0)*E_in;

I_ut_i = norm(E_ut_i)^2;
I_ut_ii = norm(E_ut_ii)^2;
I_ut_iii = norm(E_ut_iii)^2;
I_ut_iv = norm(E_ut_iv)^2;

%e)
E_ut_i = J_pol(0)*J_ret_icke_ideal(-pi/4, pi/2)*J_ret_icke_ideal(3*pi/4, pi/2)*J_pol(pi/2)*E_in;
E_ut_ii = J_pol(0)*J_ret_icke_ideal(pi/4, pi/2)*J_ret_icke_ideal(3*pi/4, pi/2)*J_pol(pi/2)*E_in;
E_ut_iii = J_pol(0)*J_ret_icke_ideal(-pi/4, pi/2)*J_ret_icke_ideal(pi/4, pi/2)*J_pol(pi/2)*E_in;
E_ut_iv = J_pol(0)*J_ret_icke_ideal(pi/4, pi/2)*J_ret_icke_ideal(pi/4, pi/2)*J_pol(pi/2)*E_in;

I_ut_i = norm(E_ut_i)^2
I_ut_ii = norm(E_ut_ii)^2
I_ut_iii = norm(E_ut_iii)^2
I_ut_iv = norm(E_ut_iv)^2

%pris: intensitet f�rloras �ven i den "korrekta" riktningen.

function matris = J_proj(alfa)
    matris = [cos(alfa), sin(alfa); -sin(alfa), cos(alfa)]; %eller andra h�llet?
end

function matris = J_pol(alfa)
    matris = J_proj(-alfa)*[1, 0; 0, 0]*J_proj(alfa);
end

function matris = J_ret(alfa, phi)
    matris = J_proj(-alfa)*[exp(1i*phi), 0; 0, 1]*J_proj(alfa);
end

function matris = J_ret_icke_ideal(alfa, phi)
    matris = J_proj(-alfa)*[exp(1i*phi*1.25), 0; 0, 1]*J_proj(alfa);
end

