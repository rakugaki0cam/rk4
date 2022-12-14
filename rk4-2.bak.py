# 4次ルンゲ・クッタ法
#   連立1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

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

def rk4(x, y, n, h):
    #4次ルンゲクッタ法によって解を求める
    # x: 変数
    # y: 解
    # n: 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
    # h: xの増加ステップ
    
    s = 4                       #次数
    c = (0, 0.5, 0.5, 1)        #4次
    b = (1/6, 1/3, 1/3, 1/6)    #重み付け　b1〜s
    k = np.zeros((n, s))
    ty = np.zeros(n)
    """
    #係数の計算
    for j in range(s):
        tx = x + c[j] * h
        for i in range(n):
            if(j == 0):
                ty[i] = y[i]
            else:
                ty[i] = y[i] + c[j] * k[i][j - 1] 

        fn = grk(tx, ty, n)
        for i in range(n):
            k[i][j] = h * fn[i]

    x += h
    for i in range(n):
        for j in range(s):
            y[i] += b[j] * k[i][j]
    """
    for i in range(n):
        yb = y[i]  
        #係数の計算
        for j in range(s):
            tx = x + c[j] * h
            if(j != 0):
                ty[i] = yb +  c[j] * k[i][j - 1]
            k[i][j] = h * grk(tx, ty, n)
            #解yの計算
            y[i] += b[j] * k[i][j]
    x += h

    return x, y


# MAIN ------------------------------------------------------------

h = 1.0e-3      #刻み
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
for i in range(0, 200):
    x = float(i / 10)
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

