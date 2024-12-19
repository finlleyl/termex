import numpy as n
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#задается интервал
step = 500
t = n.linspace(0, 10, step)
#задается функция
phi = n.sin(2*t)
theta = n.linspace(0, -2*n.pi, step)

O = 0
r2 = 0.5#внешний радиус пружины
r1 = 0.2#внутренний радиус пружины

R1 = 3 #радиус кольца

#координаты точки А
Xa = (R1 - 1) * n.cos(n.pi/2+phi)
Ya = (R1 - 1) * n.sin(n.pi/2+phi)
#координаты точки Б
Xb = - R1 * n.cos(n.pi/2+theta)
Yb = R1 * n.sin(n.pi/2+theta)

#окно и грaфик
fgr = plt.figure()
grf = fgr.add_subplot(1, 1, 1)
grf.axis('equal')
grf.set(xlim = [-5, 5], ylim = [-5, 5])
grf.set_aspect( 1 )

#треугольник
grf.plot([0,0.5,-0.5,0], [0,-1,-1,0], color = 'black')

#точки на графике
p1 = grf.plot(O, O, marker = 'o', color = 'black')[0]
pA = grf.plot(Xa[0], Ya[0], marker = 'o', color = 'black')[0]
pB = grf.plot(Xb[0], Yb[0], marker = 'o', color = 'black')[0]

#прямая ОА
O1 = grf.plot([Xa[0], O],[Ya[0], O],color = 'black')[0]

#большая окружность  
circle = plt.Circle(( O, O ), R1 , fill = False)
grf.add_artist(circle)
               
#маленькая окружность
circleA = plt.Circle(( Xa[0], Ya[0]), 1 , fill = False)
grf.add_artist(circleA)



#спиральная пружина
Ns = 2
numpnts = n.linspace(0, 1, 50*Ns+1)
Betas = numpnts*(Ns * 2*n.pi - phi[0])
Xs = ((r2-r1)*numpnts)*n.cos(Betas + n.pi/2)
Ys = ((r2-r1)*numpnts)*n.sin(Betas + n.pi/2)

SpPruzh = grf.plot(Xs, Ys, color = 'black')[0]


#анимация
def run(i):
    pA.set_data([Xa[i]], [Ya[i]])
    pB.set_data([Xb[i]], [Yb[i]])
    O1.set_data([Xa[i], O], [Ya[i], O])
    circleA.center = (Xa[i],Ya[i])

    Betas = numpnts*(Ns * 2*n.pi - phi[i])
    Xs = -((r2-r1)*numpnts)*n.cos(Betas+n.pi/2)
    Ys = ((r2-r1)*numpnts)*n.sin(Betas+n.pi/2)

    SpPruzh.set_data(Xs, Ys)
    return

anim = FuncAnimation(fgr, run, interval = 1, frames = step)

fgr.show()