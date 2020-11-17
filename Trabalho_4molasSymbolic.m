
syms ixy integer;
syms iy ix positive;
syms k1 k2 k3 k4 m1 a L b c exit rfa rfb rfc positive;
M = [ix,ixy,0;ixy,iy,0;0,0,m1];
syms z_p z_pp;
vel_z = [0,0,z_p];
acel_z = [0,0,z_pp];
syms teta beta teta_pp beta_pp t Amp z
rotx = [teta,0,0];
roty = [0,beta,0];
alphax = [teta_pp,0,0];
alphay = [0,beta_pp,0];
MM= M*(conj(alphax') + conj(alphay') + conj(acel_z'))
rk1 = [a,b,c];
rf = [rfa,rfb,rfc];
d1 = [0,0,z-(b*teta)+(a*beta)];
Fk1 = -k1*d1;
Fexit = [0,0,exit];
K = cross(rk1,Fk1) +Fk1
SK_s = collect(K, {'teta' 'beta' 'z'})