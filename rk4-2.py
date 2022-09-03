# 4次ルンゲ・クッタ法
#   連立1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

from fnmatch import fnmatch
import numpy as np
import math
import matplotlib.pyplot as plt

def grk(n, x, y):
    #計算する関数
    # x: 変数
    # y: 解
    # n: 方程式数
    f = np.zeros(n)
    ganma = 0.15

    f[0] = y[1]                         #dy/dx = v
    f[1] = -2 * ganma * y[1] - y[0]     #dv/dx = -2 * ganma * v - y
    return f

def rk4(N, x, y, h):
    #4次ルンゲクッタ法によって解を求める
    # N: 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
    # x:  変数
    # y:  解
    # h:  xの増加ステップ
    """
    #4次4段
    s = 4                       #次数
    c = (0, 0.5, 0.5, 1)        #4次        c1〜cs     
    b = (1/6, 1/3, 1/3, 1/6)    #重み付け   b1〜bs
    """
    #5次6段
    S = 6
    c = (0, 1/4, 3/8, 12/13, 1, 1/2)
    b = np.zeros((2, S))
    b[0] = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0)           #4次
    b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)    #5次
    #Rb = b[1] - b[0]
    Rb = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)
    
    ty = np.zeros(N)
    fn = np.zeros(N)
    k = np.zeros((N, S))
    y5 = np.zeros(N)
    #y5 = y とだけするとyに影響が出て計算結果がかわる？？？？
    for n in range(N):
        y5[n] = y[n]    
    
    #係数Kの計算
    for s in range(S):
        tx = x + c[s] * h
        for n in range(N):
            if s == 0:
                ty[n] = y[n]
            else:
                ty[n] = y[n] + c[s] * k[n][s - 1]
        fn = grk(N, tx, ty)
        for n in range(N):
            k[n][s] = h * fn[n]

    #解yの計算
    x += h
    for n in range(N):
        for s in range(S):
            y5[n] += b[1][s] * k[n][s]     #5次
            y[n] += b[0][s] * k[n][s]      #4次
        #print(y[n],y5[n])
    
    return x, y


# MAIN ------------------------------------------------------------

x = 0
xe = 20
h = 1.0e-3              #刻み
Nmax = int(xe / h)      #計算回数
Ns = 100                #表示点数
step = int(Nmax / Ns)   #表示間隔
N = 2
y = np.zeros(N)

#初期値
x = 0
y[0] = 1
y[1] = -0.15

gx = []
gy = []

for i in range(0, Nmax + 1):
    x, y = rk4(N, x, y, h)
    if i % step == 0:
        #表示
        print('{:8d} x:{:9.5f} y:{:13.9f} v:{:13.9f}   h:{:9.6f}'.format(i, x, y[0], y[1], h))
        gx.append(x)
        gy.append(y[0])

#グラフ表示
ganma = 0.15
v0 = -0.15
y0 = 1
A = 1       #####
alfa = 0    #####

gx2 = []
gy2 = []
for i in range(0, 20000):
    x = float(i / 1000)
    y = A * math.exp(-ganma * x) * math.cos(math.sqrt(1 - ganma ** 2) * x - alfa)
    gx2.append(x)
    gy2.append(y)

plt.plot(gx, gy, linestyle = 'None', marker = 'x', markeredgecolor = 'r')
plt.plot(gx2, gy2, color = 'g')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#  stop
#end program main

