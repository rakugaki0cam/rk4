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

def grk(x, y):
    #計算する関数
    f = np.zeros(1)
    f[0] = -y[0] * math.sin(x)
    ################### f(0) = y(0) * math.cos(x)
    return f


##TEST ---------------------------------
#a, b, c, Rc = rk_preparation('rk4')
#print(a, b, c, Rc)


def dkf45(N, x, h, y, xbound, info, tol, ik):
   #------------
   #return
   #info = -2  (Failed. calclation range has already over.)
   #     = -1  (Failed.
   #             h becomes too small. change tol or check condition of func.)
   #     =  0  (Success. running now)
   #     =  1  (Success. x reach xbound normally)
   #------------

      #implicit none
       #interface
      #   def func(iN,ix,iy,is,iik)
      #     implicit none
      #     integer,intent(in)::iN,is,iik
      #      double precision,intent(in)::ix,iy(1:iN)
      #      double precision::func
      #    end def func
      # end interface
    #integer,intent(in)::N,ik
    #double precision,intent(in)::xbound
    #double precision,intent(in)::tol
    #double precision,intent(inout)::x,h,y(1:N)
    #integer,intent(inout)::info
    #double precision,parameter::
   hmin = 1e-14
   hmax = 1
    #integer::i,j,FLAG
    #double precision::R,delta,tx,tmp(1:N),K(1:s,1:N),Sy,err
   K = np.zeros((s, N))

   #ルンゲクッタ法のブッチャーテーブル準備
   global a, b, c, Rc, s

    #5次6段
    s = 6    #次数
    a = np.zeros((s, s))
    b = np.zeros((2, s))

    c = (0, 0.25, 3/8, 12/13, 1, 0.5)
    a[0] = (0, 0, 0, 0, 0, 0)
    a[1] = (0.25, 0, 0, 0, 0, 0)
    a[2] = (0.09375, 0.28125, 0, 0, 0, 0)
    a[3] = (1932/2197, -7200/2197, 7296/2197, 0, 0, 0)
    a[4] = (439/216, -8, 3680/513, -845/4104, 0, 0)
    a[5] = (-8/27, 2, -3544/2565, 1859/4104, -11/40, 0)

    b[0] = (25/216, 0, 1408/2565, 2197/4104, -0.2, 0)  
    b[1] = (16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)
    Rc = (1/360, 0, -128/4275, -2197/75240, 1/50, 2/55)
  
      
    if(abs(h) >= hmax):
        if(h <= 0):
            h = -hmax
        else:
            h = hmax
  
    FLAG=1
    if(abs(x - xbound) <= 1e-15):
        info = 1
        FLAG = 0
    else:
        if(abs(h) <= 1e-15):
            print("maybe overflow or underflow, please change tol.")
            print("====Err info====")
            print("x    --> ", x)
            print("h    --> ", h)

            for i in range(N):
                print("y(",i,") --> ",y(i))
            
            print("================")
            info = -1
            FLAG = 0
            raise Exception("stop")

    if((h <= 0) and (xbound - x >= 0)):
        info = -2
        FLAG = 0
    elif((h > 0) and (xbound - x <= 0)):
        info = -2
        FLAG = 0


    while(FLAG == 1):
        tx = x
        for j in range(s):
            tx = x + c[j] * h #############
            tmp = y
            for i in range(j - 1):
                #tmp(1:N) = tmp(1:N) + K(i, 1:N) * a(j, i) ##############################
                temp
            for i in range(N):
                K[j, i] = h * func(N, tx, tmp, i, ik) ###############################FUNC()
                

        #step 4
        R = 0
        for i in range(N):
            R += (Rc[1] * K[1][i] + Rc[3] * K[3][i] + Rc[4] * K[4][i] + Rc[5] * K[5][i] + Rc [6] * K[6][i]) ** 2
        
        R = abs(dsqrt(R) / h)

        Sy = 0
        for i in range(N):
            Sy = Sy + (y[i] * y[i])


        Sy = dsqrt(Sy)
        if(Sy >= 1):
            err = tol * Sy
        else:
            err = tol
        
        #step 5
        if(R <= err):
            x = x + h  
            for i in range(s):
                y(1:N) = y(1:N) + b[0][i] * K(i, 1:N) ##########################
            FLAG = 0
        
        #step 6
        #  Avoid zero deviding.
        if(R >= 1e-20):
            delta = (err/(2 * R)) ** 0.25
        else:
            delta = 4

        #step 7
        if(delta <= 0.1):
            #def changes dramatically.
            h = 0.1 * h
        elif(delta >= 4):
            #def changes loosely.
            h = 4 * h
        else:
            #def changes moderately.
            h = delta * h

        #step 8
        if(abs(h) <= hmax):
            if(h <= 0):
                h = -hmax
            else:
                h = hmax

        #step 9
        if(abs(xbound-x) <= abs(h)):
            h = xbound - x
            if(abs(h) <= hmin):
                info = 1
                FLAG = 0
      
    return info #################################



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
gx2 = []
gy2 = []
for i in range(0, 200):
    x = float(i / 10)
    y = math.exp(-0.15 * x) * math.cos(math.sqrt(1 - 0.15 ** 2) * x)
    gx2.append(x)
    gy2.append(y)

plt.plot(gx, gy, linestyle = 'None', marker = 'x', markeredgecolor = 'r')
plt.plot(gx2, gy2, color = 'g')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


#  stop
#end program main

