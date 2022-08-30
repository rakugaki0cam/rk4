# 4次ルンゲ・クッタ法
#   1階微分方程式
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
    f[0] = -y[0] * math.sin(x)
    return f

def rk4(x, y, n, h):
    #4次ルンゲクッタ法によって解を求める
    # x: 変数
    # y: 解
    # n: 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
    # h: xの増加ステップ
    s = 4                       #次数
    c = (0, 0.5, 0.5, 1)        #4次　c1〜cs
    b = (1/6, 1/3, 1/3, 1/6)    #重み付け　b1〜s
    k = np.zeros((n, s))
    fn = np.zeros(n)
    ty = np.zeros(n)

    #係数の計算
    for j in range(s):
        tx = x + c[j] * h
        for i in range(n):
            if(j == 0):
                ty[i] = y[i]
            else:
                ty[i] = y[i] +  c[j] * k[i][j - 1]

        fn = grk(tx, ty, n)
        for i in range(n):
            k[i][j] = h * fn[i]

    x += h
    #解yの計算
    for i in range(n):
        for j in range(s):
            y[i] += b[j] * k[i][j]

    return x, y


# MAIN ------------------------------------------------------------

h = 1.0e-3      #刻み
Nmax = 20000    #計算回数
step = 200      #表示間隔
n = 1
y = np.zeros(n)

#初期値
x = 0
y[0] = 2

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
plt.plot(gx, gy)
plt.show()


#  stop
#end program main

