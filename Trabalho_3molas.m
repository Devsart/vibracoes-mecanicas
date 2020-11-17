%% Trabalho de Vibrações Mecânicas ( Não sei nem como começar a fazer, Deus me salve)
%% Atribuição de variaveis
theta = 0;
gamma = 0;

%%
F0 = 60000; %% Amplitude da força, dependerá da rotação da máquina
m = 5737; %%Massa do sistema
cm = [2 1.76 1.31]; %% centro de massa
I = [0.22 -0.40 0.89; -0.97 -0.05 0.22; -0.05 -0.92 -0.40]; %% Tensor de Inércia no CM
J = [6079.79 8478.34 9483.43]; %% Momentos de Inércia no CM
P_F = [1.45 1.90 1]; %% Ponto de Aplicação da Força
F = [0 0 56280]; %% Intensidade da Força na vertical (somente a amplitude)
n = 3 %% Numero de Molas

%% EDOS para cada grau de liberdade (rolagem_x, guinada_y, deslocamento_z)
%% theta = angulo de rolagem_x;
%% gamma = angulo de arfagem_y;
%% z = deslocamento vertical;
%%rot = (gamma theta 0);

%% Problemas ao usar o pacote Symbolic do Octave
%% sym z(t);
%% sym Z_pp(t);

%%Posicao das molas
  
pm1 = [3 1.76 0];
pm2 = [0 0 0];
pm3 = [0 3 0];

%% Podemos resolver as equações separadamente e somá-las, encontrar uma solução homogênea para depois encontrar a particular. Modelo W2BV-484.
k1 = 1958076;
k2 = 1958076;
k3 = 1958076;
K = [3.0976*k2+1.5376*k3, 3.52*k2-2.48*k3, -1.76*k2+1.24*k3;...
     3.52*k2-2.48*k3, 4*k2+4*k3+k1, -2*k2-2*k3+k1;...
     -1.76*k2+1.24*k3, -2*k2-2*k3+k1, k1+k2+k3];

 
 Lxx = 8359.73;
 Lxy = -257.68;
 Lyx = -257.68;
 Lyy = 8946.76;
 
 Ixx = 35914.44;
 Ixy = 19886.91;
 Iyy = 41754.68;
 
 M = [Lxx Lxy 0; Lxy Lyy 0; 0 0 m]
 
%% Analise Modal

Minv = pinv(M); %%inverso da matriz de massa;
A = Minv*K;
[U, lambda] = eig(A) %% U é a matriz de autovetores e lambda é a diagonal com os autovalores;
wn = sqrt(diag(lambda)) %% Vetor com as frequências naturais associadas a cada mola, em coluna;
Mmodal = (U'*M*U)
Mmlin = diag(Mmodal)./diag(Mmodal);
Kmodal = (U'*K*U)
Kmlin = diag(Kmodal)./diag(Mmodal);
Wnlin = sqrt(abs(Kmlin));
Fmodal = U'*F'
Flin = Fmodal./diag(Mmodal);




for i = 1:3
  Fo = Flin(i);
  wnm = Wnlin(i);
  mm = Mmlin(i);
  wf = 34*pi;
  
  x0 = 56280/Kmodal(i,i); %% ponto de equilibrio
  x0ponto = 0.1;
  f_0 = Fo/mm; 
  
  for j = 1:101
    t(j) = 5*(j-1)./100;
    w(j) = 200*(j-1)./100;
    x(j) = x0ponto*sin(wnm*t(j))/wnm + (x0 - f_0/(wnm.^2 - wf.^2))*cos(wnm*t(j)) + (f_0/(wnm.^2-wf.^2))*cos(wf*t(j)); %% eq de movimento p forcamento sem amortecimento
    T(j) = abs(Kmodal(i,i))./((-w(j).^2).*Mmodal(i,i)+abs(Kmodal(i,i)));
    XF(j) = 1./((-w(j).^2).*Mmodal(i,i)+abs(Kmodal(i,i)));
  endfor
  G = 1 ./ (((1-(wf./wnm).^2).^2).^(1/2))
  subplot(2,3,i);
  plot(t,x);
  xlabel ('t');
  ylabel ('x(t)'); %% posicao
  subplot(2,3,3+i);
  plot(w,T);
  xlabel ('w');
  ylabel ('T(w)'); %% Transmitância

endfor
  




