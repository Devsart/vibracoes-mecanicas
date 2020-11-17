syms m1 ixx ixy ixz iyy iyz izz a b c af bf cf L 
i= [ixx,ixy,ixz;ixy,iyy,iyz;ixz,iyz,izz];
syms k1 k2 k3 k4;
syms z_p z_pp;
vel_z = [0,0,z_p];
acel_z = [0,0,z_pp];
rk1 = [-a,-b,-c];
rk2 = [L-a,-b,-c];
rk3 = [L-a,L-b,-c];
rk4 = [-a,L-b,-c];
rf = [-a+af,bf-b,-c+cf];
syms teta beta teta_pp beta_pp t Amp z
rotx = [teta,0,0];
roty = [0,beta,0];
alphax = [teta_pp,0,0];
alphay = [0,beta_pp,0];
Fexit = [0,0,Amp*cos(34*pi*t)];
d1 = [0,0,z-(b*teta)+(a*beta)];
d2 = [0,0,z-(b*teta)-((L-a)*beta)];
d3 = [0,0,z+((L-b)*teta)-((L-a)*beta)];
d4 = [0,0,z+((L-b)*teta)+(a*beta)];
Fk1 = -k1*d1;
Fk2 = -k2*d2;
Fk3 = -k3*d3;
Fk4 = -k4*d4;
M = [ixx,ixy,0;ixy,iyy,0;0,0,m1];
SF = cross(rk1,Fk1)+cross(rk2,Fk2)+cross(rk3,Fk3)+cross(rk4,Fk4) +Fk1+Fk2+Fk3+Fk4;
SF_s = collect(SF, {'teta' 'beta' 'z'})
%Avet = vpa(Avet,5);
%Ava = vpa(Ava, 5);
%Ava1 = vpa(Ava(1,1),5);