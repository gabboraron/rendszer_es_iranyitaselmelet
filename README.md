# rendszer és irányításelmélet
Linkek: 
- https://regi.tankonyvtar.hu/hu/tartalom/tamop412A/2011-0042_iranyitaselmelet/ch01.html
- http://docplayer.hu/46320484-Bevezetes-rendszer-es-iranyitaselmelet.html
- 
# EA1 - bevezetés
- dirac impulzus dirac delta
![dirac impulzus és delta](https://slideplayer.hu/slide/12423963/74/images/6/Szabv%C3%A1nyos+vizsg%C3%A1l%C3%B3+jelek.jpg)
bővebben: https://slideplayer.hu/slide/12423963/
- laplace transzformáció: https://www.mateking.hu/analizis-3/laplace-transzformacio/a-laplace-transzformalt
- inverz laplace transzfomráció: https://www.mateking.hu/analizis-3/laplace-transzformacio/inverz-laplace-transzformalt

linearitásról: http://sysbook.sztaki.hu/sysbook6.php?page=54&left=intro&right=edu

# EA2 - inga - pendullum kimozdulása
> Az inga kimozdulását az alábbi egyenlet határozza meg, a [Feedback Systems - An Introduction for Scientists and Engineers könyv 85. oldalán](http://www.cds.caltech.edu/~murray/books/AM08/pdf/fbs-public_24Jul2020.pdf):
>
> ![képlet](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/ingadiffegynelet.png)

Matlab kód: [Simulation.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/Simulation.m), [balancingSystem.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/balancingSystem.m)
```Matlab
%Balancing System
M = 1;
m = 0.5; %kg
J = 0.3;
l = 0.5;
c = 0.3;
g = 9.81;
gamma = 0.2;

A = 0.1;
omega = 0.5;

tGrid = 0:0.1:50; %s

ut = tGrid;
%u = A*sin(omega*ut); %Forcing term
u = zeros(1,length(ut));
u(1) = 0.01;

x0 = [0,0.01,0,0]; %Initial condition

[t,x] = ode45(@(t,x)balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma),[0 10],x0);

figure
hold on
subplot(2,2,1);
plot(t,x(:,1:2)) %- inga kimozdulása idő függvényében
xlabel('ido')
ylabel('inga elmozdulása')
title('State Space Models - Ordinary Differential Equations')
subplot(1,2,2);
plot(t,x(:,1),t,rad2deg(x(:,2))) %- inga kimozdulása szögértékben
xlabel('ido')
ylabel('inga elmozdulása szögekben')
title('State Space Models - Ordinary Differential Equations')
hold off
grid on
```
![kimozdulás](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/inga-ido%2Bszog.jpg)

