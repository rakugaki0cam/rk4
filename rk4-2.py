# 4次ルンゲ・クッタ法
#   連立1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

from fnmatch import fnmatch
import numpy as np
import math
import matplotlib.pyplot as plt

def grk(x, y, n):
    #計算する関数
    # x: 変数
    # y: 解
    # n: 方程式数
    f = np.zeros(n)

    f[0] = y[1]
    f[1] = -0.3 * y[1] - y[0]
    return f

def rk4(x, y, N, h):
    #4次ルンゲクッタ法によって解を求める
    # x:  変数
    # y:  解
    # nn: 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
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
    b[0] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)    #5次
    b[1] = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0)           #4次
    Rc = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)
    
    k = np.zeros((N, S))
    ty = np.zeros(N)
    fn = np.zeros(N)

    #係数の計算
    for s in range(S):
        tx = x + c[s] * h
        for n in range(N):
            if(s == 0):
                ty[n] = y[n]
            else:
                ty[n] = y[n] + c[s] * k[n][s - 1]
        fn = grk(tx, ty, N)
        for n in range(N):
            k[n][s] = h * fn[n]

    #解yの計算
    for n in range(N):
        for s in range(S):
            y[n] += b[0][s] * k[n][s]       #5次
            #y[n] += b[1][s] * k[n][s]      #4次

    x += h
    
    return x, y


# MAIN ------------------------------------------------------------

h = 1.0e-3     #刻み
Nmax = 20000    #計算回数
step = 200      #表示間隔
n = 2
y = np.zeros(n)

#初期値
x = 0
y[0] = 1
y[1] = -0.15

gx = []
gy = []

for i in range(Nmax):
    x, y = rk4(x, y, n, h)
    if (i % step == 0):
        #表示
        print(x, y[0:n])
        gx.append(x)
        gy.append(y[0])

#グラフ表示
ganma = 0.15
v0 = -0.15
y0 = 1
gx2 = []
gy2 = []
for i in range(0, 2000):
    x = float(i / 100)
    y = math.exp(-ganma * x) * math.cos(math.sqrt(1 - ganma ** 2) * x)
    gx2.append(x)
    gy2.append(y)

plt.plot(gx, gy, linestyle = 'None', marker = 'x', markeredgecolor = 'r')
plt.plot(gx2, gy2, color = 'g')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#  stop
#end program main

