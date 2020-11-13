t=0:.01:20;
global omegan omegaf mo;
omegan = 6;
omegaf=17;
theta0= 0.7;
omega0=0.5;
mo=10;

function dy=forcedundampeddot(t,x)
  global omegan omegaf mo;
  dy = zeros(2,1);
  dy(1)= x(2);
  dy(2)= mo*cos(omegaf*t)-omegan^2*x(1);
end

theta=omega0/omegan*sin(omegan*t)+(theta0-(mo/(omegan^2-omegaf^2)))*cos(omegan*t)+(mo/(omegan^2-omegaf^2))*cos(omegaf*t);
thetaf=(mo/(omegan^2-omegaf^2))*cos(omegaf*t);

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4]);
[TL,YL] = ode45(@forcedundampeddot,[0 20],[theta0 omega0],options);

plot(TL,YL(:,1),'o',t,theta,'*')

