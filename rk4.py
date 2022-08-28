# 4次ルンゲ・クッタ法
#   1階微分方程式
#
# sikinote(https://slpr.sakura.ne/jp/qp) 2020/04/22によります
#

import numpy as np
import math
import matplotlib.pyplot as plt

def grk(x, y):
    #計算する関数
    f = [0, 0]
    f[0] = -y[0] * math.sin(x)
    return f

def rk4(x, y, n, h):
    #4次ルンゲクッタ法によって解を求める
    # x  --> 変数
    # y  --> 解
    # n  --> 方程式数（1階微分方程式:n=1,連立1階微分方程式:n=2)
    # h  --> xの増加ステップ
    
    c = (0, 0.5, 0.5, 1)    #4次
    k = np.zeros((n, 4))
    f = np.zeros(n)
    tmp = np.zeros(n)

    for j in range(4):
        for i in range(n):
            tmp[i] = y[i] + k[i][j] * c[j]
        tx = x + c[j] * h
        f = grk(tx, tmp)
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

