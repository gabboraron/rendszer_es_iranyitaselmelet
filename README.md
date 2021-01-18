# rendszer és irányításelmélet
Linkek: 
- https://regi.tankonyvtar.hu/hu/tartalom/tamop412A/2011-0042_iranyitaselmelet/ch01.html
- http://docplayer.hu/46320484-Bevezetes-rendszer-es-iranyitaselmelet.html
- http://www.cds.caltech.edu/~murray/books/AM08/pdf/fbs-public_24Jul2020.pdf

youtube anyagok:
- https://www.youtube.com/playlist?list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m
- https://www.youtube.com/user/ControlLectures/videos
- https://www.youtube.com/channel/UCYO_jab_esuFRV4b17AJtAw/playlists?view=50&sort=dd&shelf_id=20

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
>
> [![3BLUE1BROWN SERIES  S4 • E1 Differential equations, studying the unsolvable | DE1](https://img.youtube.com/vi/p_di4Zn4wz4/0.jpg)](https://www.youtube.com/watch?v=p_di4Zn4wz4)

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


# EA3 - autó tempomat emelkedőn
> ![kocsi ingával](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/kocsi-pendulummal.png)
>
> **Alap rendszer**
>
> ha a bemenetünk *(külső kényszerhatás)* valami *szinuszos gerjesztés* (`u(t) = A*sin(omega*t)`)
> 
> ideális eset, amikor az első(`y`) és a második(`x`) állapot egybe esik:
>
> `y(t) = x(t)`
>
> általános eset:
>
> `y(t) = h(x(t))`


> **Változtatott rendszer, hogy a kimenet nullába tartson**
>
> *Az az eset érdekel minekt mikor az origója a rendszernek egybeesik azzal, hogy  kocsi a 0 pozícióban van és a pendullum is 0 pozícióban van. Azaz mi történik ha a rendszert kitérítjük a nyugalmi állapotból és hogy választunk olyan `u(t)`-t, hogy visszatérjen a nyugalmi állapotba.*
>
> ![visszacsatolásos modell, forrás: https://mersz.hu/dokumentum/m31ireet__55](https://mersz.hu/mod/object.php?objazonosito=m31ireet_chap07_level2_sec7.2.5_i1_idx)
>
> Ha a bemenetünkön *(külső kényszerhatás)* változtatunk, visszacsatolást alkalmazva (`u(t) = -K*x`), ez egy skalár szorzás ahol `K` és `x` vektorok.
> 
> ideális eset, amikor az első(`y`) és a második(`x`) állapot egybe esik:
>
> `y(t) = x(t)`
>
> *Ezen a ponton kimérünk mindent amit lehet, hogy optimális bemenetet határozzunk meg, majd visszacsatolunk.* 
>
> Ehhez egy zárt rendszert alakítunk ki ami nem függ `u`-tól, hogy ne legyen külső beavatkozásunk. Ekkor `dx/dt=F(x,u)`-hoz `u(t) = -K*x`-t hozzávéve `dx/dt=F(x,-K*x)=G(x)` 

fájlok: [zart_rendszeru_inga/zartrendszer.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/zart_rendszeru_inga/zartrendszer.m), [zart_rendszeru_inga/balancingSystem.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/zart_rendszeru_inga/balancingSystem.m)
```Matlab
function dxdt = balancingSystem(t,x,ut,u,M,m,J,l,c,g,gamma)

K = [-1 120 -4 20]; %tankönyvi adatok, csak ezekhez a paraméterekhez működik! Minden más, különböző paraméterű rendszerhez más értékek valók!

u= -K*x;

Mt = M + m;
Jt = J + m*l^2;

dxdt = zeros(4,1);

dxdt(1) = x(3);
dxdt(2) = x(4);
dxdt(3) = (-m*l*sin(x(2))*x(4) + m*g*((m*l^2)/Jt)*sin(x(2))*cos(x(2))...
    -c*x(3)+u)/(Mt - m*((m*l^2)/Jt)*cos(x(2))^2);
dxdt(4) = (-m*l^2*sin(x(2))*cos(x(2))*x(4)^2 + Mt*g*l*sin(x(2)) +...
    c*l*cos(x(2))*x(3) + gamma*x(4) + l*cos(x(2))*u)/(Jt*(Mt/m)-m*(l*cos(x(2)))^2);

end
```
![végeredménye](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/zart_rendszeru_inga/zart_rendszeru_inga.jpg)
