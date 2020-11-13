%% Trabalho de Vibrações Mecânicas ( Não sei nem como começar a fazer, Deus me salve)
%% Atribuição de variaveis
theta = 0;
gamma = 0;

m = 5737; %%Massa do sistema
cm = [2 1.76 1.31]; %% centro de massa
I = [0.22 -0.40 0.89; -0.97 -0.05 0.22; -0.05 -0.92 -0.40]; %% Tensor de Inércia no CM
J = [6079.79 8478.34 9483.43]; %% Momentos de Inércia no CM
P_F = [1.45 1.90 1]; %% Ponto de Aplicação da Força
F = [0 0 0]; %% Intensidade da Força na vertical
n = input("Quantas molas deseja? (Escolha entre 3 e 6):");
%% Analise Modal

M = 0;  %%matriz de massa;
K = 0;  %%matriz de Rigidezas;
Minv = pinv(M) %%inverso da matriz de massa;

[U, lambda] = eig(Minv*K); %% U é a matriz de autovetores e lambda é a diagonal com os autovalores;
wn = sqrt(diag(lambda)) %% Vetor com as frequências naturais associadas a cada mola, em coluna;

%% EDOS para cada grau de liberdade (rolagem_x, guinada_y, deslocamento_z)
%% theta = angulo de rolagem_x;
%% gamma = angulo de guinada_y;
%% z = deslocamento vertical;
rot = (gamma theta 0);


%% Podemos corrigir o Tensor de Inércia utilizando a diagonalização 
if(n = 4)
  %%Posicao das molas
  m1 = [0 0 0];
  m2 = [3 0 0];
  m3 = [0 3 0];
  m4 = [3 3 0];
elseif(n=5)
  m1 = [0 0 0];
  m2 = [3 0 0];
  m3 = [0 3 0];
  m4 = [3 3 0];
  m5 = [1.5 1.5 0];
  z5 = z -(m5-cm)*rot';
elseif(n=6)
  m1 = [0 0 0];
  m2 = [3 0 0];
  m3 = [0 3 0];
  m4 = [3 3 0];
  m5 = [1.5 0 0];
  m6 = [1.5 3 0];
  z5 = z -(m5-cm)*rot';
  z6 = z -(m6-cm)*rot';
endif

  z1 = z -(m1-cm)*rot'
  z2 = z -(m2-cm)*rot'
  z3 = z -(m3-cm)*rot'
  z4 = z -(m4-cm)*rot'
  pf = P_F - cm;

%% Podemos resolver as equações separadamente e somá-las, encontrar uma solução homogênea para depois encontrar a particular
 m*z_pp - F(3) + k1*z1+ k2*z2+ k2*z2+ k4*z4+ k5*z5 +k6*z6 = 0
