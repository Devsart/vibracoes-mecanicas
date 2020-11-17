m1=5737;
ixx=8359.73;
iyy=8946.76;
izz=6735.07;
ixy=-257.68;
ixz=460.34;
iyz=-1211.23;
i= [ixx,ixy,ixz;ixy,iyy,iyz;ixz,iyz,izz];
a=2;
b=1.76;
c=1.31;
af=1.45;
bf=1.90;
cf=1;
Lx=3;
Ly=3;
L=3;
k1=1117958;
k2=1117958;
k3=1117958;
k4=1117958;
syms z_p z_pp;
vel_z = [0,0,z_p];
acel_z = [0,0,z_pp];
rk1 = [-a,-b,-c];
rk2 = [Lx-a,-b,-c];
rk3 = [Lx-a,Ly-b,-c];
rk4 = [-a,Ly-b,-c];
rf = [-a+af,bf-b,-c+cf];
syms teta beta teta_pp beta_pp t Amp z
rotx = [teta,0,0];
roty = [0,beta,0];
alphax = [teta_pp,0,0];
alphay = [0,beta_pp,0];
Fexit = [0,0,Amp*cos(34*pi*t)];
d1 = [0,0,z-(b*teta)+(a*beta)];
d2 = [0,0,z-(b*teta)-((Lx-a)*beta)];
d3 = [0,0,z+((Ly-b)*teta)-((Lx-a)*beta)];
d4 = [0,0,z+((Ly-b)*teta)+(a*beta)];
Fk1 = -k1*d1;
Fk2 = -k2*d2;
Fk3 = -k3*d3;
Fk4 = -k4*d4;
M = [ixx,ixy,0;ixy,iyy,0;0,0,m1]
SK = cross(rk1,Fk1)+cross(rk2,Fk2)+cross(rk3,Fk3)+cross(rk4,Fk4) +Fk1+Fk2+Fk3+Fk4;
SK_s = collect(SK, {'teta' 'beta' 'z'});
K =[ -(- b^2*k1 - b^2*k2 - k3*(L - b)^2 - k4*(L - b)^2),-(a*b*k1 - b*k2*(L - a) - a*k4*(L - b) + k3*(L - a)*(L - b)),-(b*k1 + b*k2 - k3*(L - b) - k4*(L - b))
-(a*b*k1 - b*k2*(L - a) - a*k4*(L - b) + k3*(L - a)*(L - b)),-(- a^2*k1 - a^2*k4 - k2*(L - a)^2 - k3*(L - a)^2),- (k2*(L - a) - a*k4 - a*k1 + k3*(L - a))
- (b*k1 + b*k2 - k3*(L - b) - k4*(L - b)),- (k2*(L - a) - a*k4 - a*k1 + k3*(L - a)),- (- k1 - k2 - k3 - k4)]
SF = cross(rf,Fexit) + Fexit;
SF_ = SF';
Mi = inv(M)
D = Mi*K
[Avet,Ava]=eig(D)
Wn=sqrt(Ava)
Fn=Wn/(2*pi)
Km=Avet'*K*Avet
Mm=Avet'*M*Avet
Untitled5

M =

   1.0e+03 *

    8.3597   -0.2577         0
   -0.2577    8.9468         0
         0         0    5.7370


K =

   1.0e+07 *

    1.0364   -0.0581   -0.1163
   -0.0581    1.1180    0.2236
   -0.1163    0.2236    0.4472


Mi =

   1.0e-03 *

    0.1197    0.0034         0
    0.0034    0.1119         0
         0         0    0.1743


D =

   1.0e+03 *

    1.2388   -0.0311   -0.1315
   -0.0293    1.2487    0.2461
   -0.2027    0.3897    0.7795


Avet =

    0.1729    0.8690    0.4366
   -0.3414    0.4934   -0.7159
    0.9239    0.0371   -0.5449


Ava =

   1.0e+03 *

    0.5975         0         0
         0    1.2156         0
         0         0    1.4539


Wn =

   24.4441         0         0
         0   34.8653         0
         0         0   38.1298


Fn =

    3.8904         0         0
         0    5.5490         0
         0         0    6.0685


Km =

   1.0e+07 *

    0.3717    0.0000   -0.0000
    0.0000    1.0063   -0.0000
   -0.0000   -0.0000    1.1693


Mm =

   1.0e+03 *

    6.2200    0.0000   -0.0000
    0.0000    8.2781   -0.0000
   -0.0000   -0.0000    8.0429



