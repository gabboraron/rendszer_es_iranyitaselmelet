# rendszer és irányításelmélet
Linkek: 
- https://regi.tankonyvtar.hu/hu/tartalom/tamop412A/2011-0042_iranyitaselmelet/ch01.html
- http://docplayer.hu/46320484-Bevezetes-rendszer-es-iranyitaselmelet.html
- http://www.cds.caltech.edu/~murray/books/AM08/pdf/fbs-public_24Jul2020.pdf

youtube anyagok:
- https://www.youtube.com/playlist?list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m
- https://www.youtube.com/user/ControlLectures/videos
- https://www.youtube.com/channel/UCYO_jab_esuFRV4b17AJtAw/playlists?view=50&sort=dd&shelf_id=20

# Bevezetés
- dirac impulzus dirac delta: egy olyan függvény ami egy végtelen magas, de egység területű végtelenükl kicsire összenyomott sézlességű téglalappot ad.
![dirac impulzus és delta](https://slideplayer.hu/slide/12423963/74/images/6/Szabv%C3%A1nyos+vizsg%C3%A1l%C3%B3+jelek.jpg)
bővebben: https://slideplayer.hu/slide/12423963/
- laplace transzformáció: https://www.mateking.hu/analizis-3/laplace-transzformacio/a-laplace-transzformalt
- inverz laplace transzfomráció: https://www.mateking.hu/analizis-3/laplace-transzformacio/inverz-laplace-transzformalt

linearitásról: http://sysbook.sztaki.hu/sysbook6.php?page=54&left=intro&right=edu

# Inga - pendullum kimozdulása
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


# Autó tempomat emelkedőn
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

## Tempomat rendszer
![blokkdiagram](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/cruise_control.png)
- `v`  - kimenetei sebesség
- `vr` - reefeerncia sebesség
- `F` - a kocsit mozgató erő
- `Fd` - zavaró erő, pl egy rámpa
- `T` - nyomaték a kerekeknél
- `phi` - a gázkar változása szögben, ez hat a motorra 

> Ennek modellezéséhez használjuk [Newton II-t](https://hu.wikipedia.org/wiki/Newton_t%C3%B6rv%C3%A9nyei#Newton_II._t%C3%B6rv%C3%A9nye_%E2%80%93_a_dinamika_alapt%C3%B6rv%C3%A9nye): `tomeg * gyorsulás = a testre ható erők összegével` => `m*dv/dt = F-Fd`, kérdés mik az `F`, `Fd` erők? Aaz hogyan változik a nyomaték a szögsebesség függvényében?
>
> Ehhez használjuk a motor forgatónyomatékát![forgatónyomaték és szögsebesség](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/forgat%C3%B3nyomat%C3%A9k_es_szogsebesseg.png)
>
> a kerék surlódását a földdel ![KERÉK SURLÓDÁSA FÖLDDEL](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/ker%C3%A9k_sul%C3%B3d%C3%A1s.png) és az aerodinamikai hatást, erőt ![aerodinamikai hatás](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/aerodinamikai_ero.png) valamint a gravitáció sem elhanyagolható: ![gravitációs erő](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/gavitacios_ero.png)

Teljes forrásfájlok: [cruise_control/imulation-Cruisecontrol.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/simulation-Cruisecontrol.m), [cruise_control/cruiseControl.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/cruise_control/cruiseControl.m)
```Matlab
%CSAK RÉSZLET
ut = tGrid;
u = ones(1,length(ut))*0.5;

thetat = tGrid;
theta = zeros(1,length(thetat));

iStart = find(tGrid == 20);  %bukkanó megadása
iEnd = find(tGrid == 30);  %

theta(iStart:end) = ones(1,length(theta(iStart:end)))*deg2rad(angle); %
theta(iStart:iEnd) = interp1([thetat(iStart),thetat(iEnd)],[0,deg2rad(angle)],thetat(iStart:iEnd)); % lineáris interpoláció, hogy a bukkanó ne éles szögben váltson
                                                                                                    % ismerjük a két pont koordinátáit és az x koordinátákat a kettő között

theta = zeros(1,length(thetat));


x0 = [v0,0]; %Initial condition

[t,x] = ode15s(@(t,x)cruiseControl(t,x,Tm,beta,alphan,omegam,m,g,Cr,rho,Cd,A,ut,u,thetat,theta,vref),[0 tEnd],x0); %integrálja a korábbi hibákat, hogy összeadja őket, hogy emlékezzen rájuk.

figure
plot(t,x)
grid on

figure
plot(thetat,theta)
grid on
```

## Fászisportré - fázisdiagram
> ### Mi a [skalármező](https://en.wikipedia.org/wiki/Scalar_field)?
> 
> A `z=f(x,y)` függvény a két dimenziós `R` számokhoz képez `R` számokat. Ekkor ha egy háromdimenziós térben egy síkon levő ponthoz `(x1,y1)` társítok egy harmadik pontot. Pl egy koordinátához egy hőmérsékletet. Pl:
> ![skalármező](https://www.researchgate.net/profile/Miguel_Negrao/publication/229038491/figure/fig4/AS:300840888881156@1448737465473/An-example-of-a-scalar-field-used-in-the-Sine-Field-system.png)

> ### [Vektormező](https://en.wikipedia.org/wiki/Vector_field#:~:text=In%20vector%20calculus%20and%20physics,a%20point%20in%20the%20plane.)
> 
> Ha van egy két dimenziós `R`-ből két dimenziós `R`-be képező `(z,w) = F(x,y)` függvény. Ez az előbbitől eltérően minden egyes ponthoz egy vektor értéket tesztel, nem egy skalár értéket, azaz, minden egyes ponthoz egy vektorértéket társít. Pl, egy [széltérkép](https://www.idokep.eu/ceu/szel). Ezek a vektorok nem csak síkbeliek, de térbeliek is lehetnek!
>
> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/b9/VectorField.svg/1024px-VectorField.svg.png"  width="400" height="400" /> <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Cessna_182_model-wingtip-vortex.jpg/1920px-Cessna_182_model-wingtip-vortex.jpg" width="400" height="400" />
> 
> Tehát ebben az esetben egy olyan függvényről beszélünk aminek két bemenete van és két kimenete.
> ```Matlab
> function dxdt = massSpringDamper(t,x,ut,u,m,k,c)
> u = interp1(ut,u,t);
> dxdt = zeros(2,1);
> dxdt(1) = x(2);
> dxdt(2) = -(c/m)*x(2) - (k/m)*x(1) + u/m;
> end
> ```
> Ezzel egy egyensúlyi állpot alakul ki, csakhogy ez egy instabil egyensúlyi álllapot lesz.

fájlok: [fazisdiagram/Simulation-fazis.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/fazisdiagram/Simulation-fazis.m), [fazisdiagram/massSpringDamper.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/fazisdiagram/massSpringDamper.m)
```Matlab
% %Mass spring damper
m = 0.5; %kg
c = 0.1;
k = 0.3;

A = 0;
omega = 0.5;

tGrid = 0:0.01:7; %s

ut = tGrid;
u = A*sin(omega*ut); %Forcing term

x0 = [0.01,0]; %Initial condition

[t,x] = ode45(@(t,x)massSpringDamper(t,x,ut,u,m,k,c),tGrid,x0);


xRange = -1:0.1:1;
yRange = -1:0.1:1;

[X,Y] = meshgrid(xRange,yRange); % az összes xRange és yRange közti értékre elvégzi a Descartes szorzást

dxdt = Y;
dydt = -(c/m)*Y - (k/m)*X;

figure
plot(t,x)
grid on

figure
plot(x(:,1),x(:,2)) %fázisállapot, ahol x tengelyen az első állapotváltozó, y tengelyen a második
title('fázisállaopotok')
xlabel("elso állapot változás")
ylabel("második állapot változása")
grid on

figure
quiver(X,Y,dxdt,dydt) %minden egyes pontpárhoz hozzárendel vektort
title('vektor rendelése minden egyes pontpárhoz')
grid on
hold on
plot(x(:,1),x(:,2))
hold off
```

> ### invertált pendulum
>
> <img src="https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/invertalt_pendulum.png"  width="200" height="400" /> 
> Az előző probléma részben leegyszerűsítve, itt nincs nincs kis kocsi, , rögzített az inga, ha hatunk a csúcsára akkor változtatja a pozícióját, és az elmozduás szge az ami érdekes számunkra. Ehhez ezt a képletet kell  használnunk: <img src="https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/keplet.png"  width="200" height="50" /> Ebből jól látszik, hogy az ingának két egynesúlyi állapota van, az első instabil, ebből megyy át a msáodikba ami stabil *(miután megbökjük és leesik)*. Ezt látjuk alább:
>
> <img src="https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/inga_allapotvaltozasa.jpg"  width="400" height="400" />  <img src="https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/pos.jpg"  width="400" height="400" /> 
>
> Ezek alapján felmerül mi az egyensúlyi állapot?
>
> **Egyensúlyi állapot: Van egy rendszer ami beáll egy konstans értékre amiből nem változik a rendszer állapot akkor ha a konstanst deriváljuk 0-át kapunk. Tehát a rendszer megváltozása `dx/dt = 0` amikor egyensúlyi állapotban vagyunk.**

fájlok: [invert_pendulum/invertedPendulum.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/invertedPendulum.m), [invert_pendulum/simulation-inverted_pendulum.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/invert_pendulum/simulation-inverted_pendulum.m)

> ### [Lorenz model](https://en.wikipedia.org/wiki/Lorenz_system)
>
> Amikor egy érdekes felületen kerül egyensúlyi állapotba a rendszer.
```Matlab
% Solve over time interval [0,100] with initial conditions [1,1,1]
% ''f'' is set of differential equations
% ''a'' is array containing x, y, and z variables
% ''t'' is time variable

sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,a] = ode45(f,[0 100],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
plot3(a(:,1),a(:,2),a(:,3))
```

# Egyensúlyi állapotok
> **Equilibrium pont: egy olyan kitüntett pont ahol a derivált az  konstans 0 lesz. Ez az egyensúlyi állapot.**
>
> **stabil megoldás: `a` és `b`kezdeti feltételek, `x(t,a)` stabil minden `E>0`-ra, ha van olyan `d>0`, hogy a kettő közötti távolság(`||b-a||`) kisebb mint `d` akkor `||x(t,b)-x(t,a)||<E` minden `t>0`-ra azaz mindörökre.** Tehát gyakorlatilag ha van egy tengelyünk ami körül mozog egy légy, ha tudunk egy olyan hengert adni aminél messszebb nem megy a légy a tengelytől akkor ez egy stabil rendszer.
>
> ![stabilitás meghatásrozása](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-4/stabil_keplet.png)
>
> **a renszer instabil ha nincs ilyen `d`**
>
> **lokálisan stabil: az a része a rendszernek ami egy adott területen stabil**, például van egy pontunk, és körülötte egy adott körben minden abba az egy pontba tart, de a körön kívűl eső pontok más irányba tartanak. Ilyen lehet mondjuk egy örvény, turbulencia. Ha egy ilyen kör radiusza végtelen nagy akkor globálisan stabil.


## Lineáris rendszerek
> Ebben az esetben az origó mindig equilibrium pont lesz. Ehhez a `dx/dt = Ax` egyenletben levő `A` mátrix saját értékeit használjuk fel. Ezek komlex sázmok lesznek. Ha a saját érték valós része `negtív` akkor a rendszer stabil, ha `nulla` akkor nem biztos, hogy stabil,  ha pozitív akkor biztos, hogy nem stabil.
>
> **Egy lináris rendszer a `dx/dt = Ax`formában aszimptotikusan stabil akkor és csak akkor, ha az összes sajátértéke az `A`-nak szigorúan negatív valós részű és instabil ha bármely sajátértéke az `A`-nak szigorúan pozitív valós részű.

> ## [DC motor](https://ctms.engin.umich.edu/CTMS/index.php?example=MotorSpeed&section=SystemModeling) modellezése
>
> ![motor modell](https://ctms.engin.umich.edu/CTMS/Content/MotorSpeed/System/Modeling/figures/motor.png)
> 
> ahol `T ~ theta` nyomaték. `dtheta/dt = szögsebesség = thetav` és `ddtheta/dt = szöggyorsulás = dthetav/dt`
>
> Erre használhatjuk a [Kisrchoff](https://hu.wikipedia.org/wiki/Kirchhoff-t%C3%B6rv%C3%A9nyek) szabályt ésa Newton II szabály, és így <img src="https://ctms.engin.umich.edu/CTMS/Content/MotorSpeed/System/Modeling/html/MotorSpeed_SystemModeling_eq17969768335097513228.png"  width="100" height="20" />  és <img src="https://ctms.engin.umich.edu/CTMS/Content/MotorSpeed/System/Modeling/html/MotorSpeed_SystemModeling_eq13914453465500333701.png"  width="100" height="20" /> azaz `dthetav = (K/J)*i - (b/J)*thetav` és `di/dt = (1/L)*V-(K/L)*thetav - (R/L)*i` ahol `V` feszültség az input, tehát nem zárt a rendszer. Ez általános alakban: `dx/dt = F(x,u)` azaz `F` az állapot *(`x`)* és input *(`u`)* fügvénye. Ez egy lineáris rendszerre *(`A*x+B*u` alakra)* is átírható => 
> ```
> ( - (b/J)   (K/J) )  ( thetav )  +  (  0  ) v 
> ( - (K/L)  -(R/L) )  (   i    )     ( 1/L )
>
>        A                x       +      B    u 
> ```
>
> ```Matlab
> function dxdt = DCmotor(t,x,ut,u,J,b,Ke,Kt,R,L)
> u = interp1(ut, u, t);
> dxdt = zeros(2,1);
> A = [-b/J, Ke/J; -Ke/L, -R/L];
> B = [0; 1/L];
> dxdt = A*x + B*u;
> end
> ```

fájlok: [Órai anyagok-4/DCmotor.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-4/Live_simulation.m), [Órai anyagok-4/Live_simulation.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-4/Live_simulation.m)


# Nem lineáris rendszerek stabilitása
> Ehhez [Lajpunov függvényeket](http://phys.ubbcluj.ro/~zneda/nemlin-math/c6.pdf) fogunk használni, amik pozitív definit függvények. A lineráis rendszernek csak egyetlen egyensúlyi pontja van ,az origóban *(kivéve ha többszörös sajátértékei vannak, akkor végtelen sok)*, míg a nem lineáris rendszerekhez lineáris approximációra van szükség. Ehhez [Taylosr sorokat](https://en.wikipedia.org/wiki/Taylor_series) fogunk használni, így a be tudunk helyetesíteni `sin(x1)` helyére `x1`-et, és ugyanígy `cos(x2)`-vel is. Ez viszont már egy lineáris rendszer! Inentől ismét egy lineársi rendszer stabilitása a kérdés.
>
> A kérdés a stabil egyensúlyi pontban `(pi,0)` pontban kérdéses. Tehát, ha `z1 = x1-pi` => `x1 = z1 + pi` és `z2 = x2-0` => `x2 = z2` akkor `d(z1 + pi)/dt = z2` azaz `dz1/dt + dpi/dt = dz1/dt` és `dz2/dt = sin(z1 + pi) -c * z2 = -z1 - c*z2`
>
> Általános esetbenben a kimenenet egy lineáris rendszer. Az alábbi képlettel adható meg:
>
> ![egyenlet](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-5/egyenlet.png)

fájlok: [Órai anyagok-5/Live_simulation.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-5/Live_simulation.m), [Órai anyagok-5/DCmotor.m](https://github.com/gabboraron/rendszer_es_iranyitaselmelet/blob/main/%C3%93rai%20anyagok-5/DCmotor.m)
```Matlab
A = [-b/J, Ke/J; -Ke/L, -R/L];
B = [0; 1/L];
C = [1,0];
D = 0;

sys = ss(A,B,C,D); % folytonos idejű állapotteres modellt ad
step(sys); % step response kiplottolva
impulse(sys); % impulzus válasza a rendszernek
```


