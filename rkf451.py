# ルンゲ・クッタ=フェールベルグ法
#   刻み幅自動制御
#
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
    ganma = 0.15

    f[0] = y[1]                         #dy/dx = v
    f[1] = -2 * ganma * y[1] - y[0]     #dv/dx = -2 * ganma * v - y
    return f


def dkf45(x, y, N, tol, h):
    #帰り値
    #info = -2  異常　計算範囲をオーバー
    #     = -1  異常　刻み値が小さくなりすぎ　tolの再検討が必要
    #     =  0  正常進行中
    #     =  1  x終了値に到達し完了

    hmin = 1e-14
    hmax = 1
    info = 0
    key = 0

    #ルンゲクッタ法のブッチャーテーブル準備
    #global a, b, Rc

    #5次6段
    S = 6    #段数
    """
    a = np.zeros((S, S))
    a[0] = (0, 0, 0, 0, 0, 0)
    a[1] = (1/4, 0, 0, 0, 0, 0)
    a[2] = (3/32, 9/32, 0, 0, 0, 0)
    a[3] = (1932/2197, -7200/2197, 7296/2197, 0, 0, 0)
    a[4] = (439/216, -8, 3680/513, -845/4104, 0, 0)
    a[5] = (-8/27, 2, -3544/2565, 1859/4104, -11/40, 0)
    #cの計算
    for i in range(6):
        print("c", i, "=", sum(a[i]))
    """
    b = np.zeros((2, S))
    b[0] = (25/216, 0, 1408/2565, 2197/4104, -1/5, 0)           #4次
    b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)    #5次
    c = (0, 1/4, 3/8, 12/13, 1, 1/2)
    #Rc = b[1] - b[0]
    Rc = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)         # 検算用

    tx = np.zeros(N)
    ty = np.zeros(N)
    tf = np.zeros(N)
    K = np.zeros((N, S))
    R = np.zeros(N)
    Rcnt = 0

    if(abs(h) >= hmax):
        if(h <= 0):
            h = -hmax
        else:
            h = hmax
  
    FLAG = 1
    if(abs(x - xbound) <= hmin):
        info = 1
        FLAG = 0
    else:
        if(abs(h) <= 1e-15):
            print("刻み値が小さくなりすぎ。tolの変更が必要")
            print("x =", x)
            print("h =", h)
            for n in range(N):
                print("y(",n,") =", y(n))
            print("================")
            info = -1
            FLAG = 0
            raise Exception("stop")

    while(FLAG == 1):
        #係数Kの計算
        for s in range(S):
            tx = x + c[s] * h
            for n in range(N):
                if(s == 0):
                    ty[n] = y[n]
                else:
                    ty[n] = y[n] + c[s] * K[n][s - 1]
            tf = grk(tx, ty, N)
            for n in range(N):
                K[n][s] = h * tf[n]
                
        #step 4
        #4次での解と5次での解の差を求める
        # R = |x(4次)-x(5次)|/h
        R = 0
        r = np.zeros(N)
        for n in range(N):
            for s in range(S):
                r[n] += Rc[s] * K[n][s]
            R += r[n] ** 2
        R = abs(math.sqrt(R) / h / N)
        Rcnt += 1
        #print(n, "R =", R)

        Sy = 0
        for n in range(N):
            Sy += y[n] ** 2
        Sy = math.sqrt(Sy)
        if(Sy >= 1):
            err = tol * Sy
        else:
            err = tol
        
        #step 5
        #4次での解を計算
        y5 = np.zeros(N)
        y5 = y
        if((R <= err) or (key == 1)):
            #print('R計算回数', Rcnt)
            x = x + h  
            for n in range(N):
                for s in range(S):
                    y[n] += b[0][s] * K[n][s] 
                    #y5[n] += b[1][s] * K[n][s] #5次

                #print(y[n],y5[n], '差',n,y[n]-y5[n])

            FLAG = 0
        
        #step 6
        #  Avoid zero deviding.
        if(R >= 1e-20):
            delta = (err / (2 * R)) ** (1 / 4)
        else:
            delta = 4

        #step 7
        if(delta <= 0.1):
            #現在の刻み幅hは大きすぎる
            h = 0.1 * h
        elif(delta >= 4):
            #現在の刻み幅は小さすぎる
            h = 4 * h
        else:
            #刻み幅を適正値に修正する
            h = delta * h

        #step 8
        if(abs(h) >= hmax):
            #刻み幅がmaxを超えないように
            if(h <= 0):
                h = -hmax
            else:
                h = hmax

        #step 9
        if(abs(xbound - x) <= abs(h)):
            #hをx終端値に合わせる
            h = xbound - x
            if(abs(h) <= hmin):
                #正常終了
                info = 1
                FLAG = 0

        #計算が範囲外になっている場合
        if((h <= 0) and ((xbound - x) >= 0)):
            info = -2
            FLAG = 0
        elif((h > 0) and ((xbound - x) <= 0)):
            info = -2
            FLAG = 0
        
    if(key == 1):
        print("Strange point between ", x - h, " and ", x)
        info = -9

    return x, y, info, h



# MAIN ------------------------------------------------------------

x = 0
xbound = 20     #x終了値
#初期値
N = 2
y = np.zeros(N)
y[0] = 1
y[1] = -0.15

gx = []
gy = []

tol = 1e-8  #計算精度
h = xbound - x  #刻み
i = 0
info = 0
step = 200

while(info <= 0):
    x, y, info, h = dkf45(x, y, N, tol, h)
    i += 1
    if((i % step == 0) or (info == 1)):     #最終データも表示
        print('{:8d} x:{:9.5f} y:{:13.9f} v:{:13.9f} h:{:9.6f}'.format(i, x, y[0], y[1], h))
        gx.append(x)
        gy.append(y[0])
print('計算回数 {:8d}回'.format(i))    


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

