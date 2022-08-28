# 4次ルンゲ・クッタ法
#   連立1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

import numpy as np
import math
import matplotlib.pyplot as plt

def grk(x, y, f):
    #計算する関数
    f[0] = y[1]
    f[1] = -0.3 * y[1] - y[0]
    return x, y, f

def rk4(x, y, n, h):
    #4次ルンゲクッタ法によって解を求める
    # x  --> 変数
    # y  --> 解
    # n  --> 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
    # h  --> xの増加ステップ
    
    c = (0, 0.5, 0.5, 1)    #4次
    k = []
    f = []
    tmp = []
    
    for i in range(n):
        k.append([0, 0, 0, 0])
        f.append(0)
        tmp.append(0)

    for j in range(4):
        for i in range(n):
            tmp[i] = y[i] + k[i][j] * c[j]
        tx = x + c[j] * h
        _, _, f = grk(tx, tmp, f)
        for i in range(n):
            k[i][j] = h * f[i]

    x = x + h
    for i in range(n):
        y[i] = y[i] + (k[i][0] + k[i][3]) / 6 + (k[i][1] + k[i][2]) / 3

    return x, y


# MAIN ------------------------------------------------------------

h = 1.0e-3      #刻み
Nmax = 20000    #計算回数
step = 200      #表示間隔
n = 2
y = [0, 0]

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
        print(x, y[0], y[1])
        gx.append(x)
        gy.append(y[0])

#グラフ表示
plt.plot(gx, gy)
plt.show()


#  stop
#end program main

