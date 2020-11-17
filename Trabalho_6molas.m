%% Trabalho de Vibrações Mecânicas ( Não sei nem como começar a fazer, Deus me salve)
%% Atribuição de variaveis
%%
F0 = 56280; %% Amplitude da força, dependerá da rotação da máquina
m = 5737; %%Massa do sistema
cm = [2 1.76 1.31]; %% centro de massa
I = [0.22 -0.40 0.89; -0.97 -0.05 0.22; -0.05 -0.92 -0.40]; %% Tensor de Inércia no CM
J = [6079.79 8478.34 9483.43]; %% Momentos de Inércia no CM
P_F = [1.45 1.90 1]; %% Ponto de Aplicação da Força
F = [0 0 F0]; %% Intensidade da Força na vertical (somente a amplitude)
n = 6 %% Numero de Molas

%% EDOS para cada grau de liberdade (rolagem_x, guinada_y, deslocamento_z)
%% theta = angulo de rolagem_x;
%% gamma = angulo de arfagem_y;
%% z = deslocamento vertical;

%%Posicao das molas
  
pm1 = [3 3 0];
pm2 = [3 0 0];
pm3 = [1.5 0 0];
pm4 = [0 0 0];
pm5 = [0 3 0];
pm6 = [1.5 3 0];

%% Podemos resolver as equações separadamente e somá-las, encontrar uma solução homogênea para depois encontrar a particular. Modelo W2BV-481.
k1 = 945684;
k2 = 945684;
k3 = 945684;
k4 = 945684;
k5 = 945684;
k6 = 945684;



K = [k1+k2+k3+k4+k5+k6, k1+k2-0.5*k3-2*k4-2*k5-0.5*k6,1.24*k1-1.76*k2-1.76*k3-1.76*k4+1.24*k5+1.24*k6;...
      k1+k2-0.5*k3-2*k4-2*k5-0.5*k6, k1+k2+0.25*k3+4*k4+4*k5+0.25*k6, 1.24*k1-1.76*k2+0.88*k3+3.52*k4-2.48*k5-0.62*k6;...
      1.24*k1-1.76*k2-1.76*k3-1.76*k4+1.24*k5+1.24*k6, 1.24*k1-1.76*k2+0.88*k3+3.52*k4-2.48*k5-0.62*k6, 1.5376*k1+2.1824*k2+2.1824*k3+2.1824*k4+1.5375*k5+1.5376*k6]
      
      

 
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
wn = sqrt(abs(diag(lambda))) %% Vetor com as frequências naturais associadas a cada mola, em coluna;
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