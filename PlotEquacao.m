V0=0.5; %%Velocidade Inicial (m/s)
P=25627; %%Peso da massa do sistema (Kgf)
k=64824; %%Rigidez vertical do sistema (N/m)
c=0; %%Característica do amortecedoor (N*s*m^-1)
g=9.8; %%Gravidade (m/s^2)
omega_n=sqrt(k/P); %%Frequencia Natural (rad/s)
cc=2*P*omega_n; %%Coeficiente de Amortecimento Critico 
xi=c/cc; %%Fator de Amortecimento (N*s*m^-1)
omega_d=sqrt(1-(xi.^2))*omega_n; %%Frequencia Natural Amortecida (rad/s)

A1=-(e.^(-xi*omega_n*0)*0*sin(omega_d*0))/(e.^(-xi*omega_n*0)*cos(omega_d));
A2=(A1*(-xi*omega_n*(e.^(-xi*omega_n*0))*cos(omega_d*0)-omega_d*(e.^(-xi*omega_n*0))*sin(omega_d*0))-V0)/(-xi*omega_n*(e.^(-xi*omega_n*0))*sin(omega_d*0)-omega_d*(e.^(-xi*omega_n*0))*cos(omega_d*0));
t=0:.01:10;
u=e.^(-xi*omega_n*t).*(A1*cos(omega_d*t)+A2*sin(omega_d*t));
plot(t,u,':b')