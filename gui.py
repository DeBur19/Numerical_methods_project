from tkinter import *
from tkinter import filedialog as fd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import sys


def integral_trap(x, y, N):
    ans = 0
    h = (x[N] - x[0]) / (2 * N)
    for i in range(N):
        ans += (y[i] + y[i + 1])
    return ans * h


def cubic_spline(xgrid, fgrid):
    x = np.array(xgrid)
    f = np.array(fgrid)
    N = x.shape[0] - 1
    a = np.zeros(N + 1)
    a[1:] = f[:N]
    b = np.zeros(N + 1)
    c = np.zeros(N + 1)
    d = np.zeros(N + 1)
    h = np.zeros(N + 1)
    h[1:] = x[1:] - x[:N]
    p = np.zeros(N)
    q = np.zeros(N)
    F = np.zeros(N)
    F[1:] = 3 * (((f[2:] - f[1:N]) / h[2:]) - ((f[1:N] - f[:N - 1]) / h[1:N]))
    p[2] = -h[2] / (2 * (h[1] + h[2]))
    q[2] = F[1] / (2 * (h[1] + h[2]))

    # Метод прогонки
    # Вычисление коэффицентов прогонки:
    for i in range(2, N - 1):
        p[i + 1] = -h[i + 1] / (h[i] * p[i] + 2 * (h[i] + h[i + 1]))
        q[i + 1] = (F[i] - h[i] * q[i]) / (h[i] * p[i] + 2 * (h[i] + h[i + 1]))

    # Вычисление решения системы:
    c[N] = (F[N - 1] - h[N - 1] * q[N - 1]) / (2 * (h[N - 1] + h[N]) + h[N - 1] * p[N - 1])
    for i in range(N - 1, 1, -1):
        c[i] = p[i] * c[i + 1] + q[i]

    # Вычисление коэффицентов b, d:
    b[N] = ((f[N] - f[N - 1]) / h[N]) - 2 * h[N] * c[N] / 3
    d[N] = -c[N] / (3 * h[N])
    for i in range(1, N):
        b[i] = ((f[i] - f[i - 1]) / h[i]) - (h[i] * (c[i + 1] + 2 * c[i])) / 3
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])
    return a, b, c, d


def interpolated_func(arg, xgrid, a, b, c, d):
    x = np.array(xgrid)
    for i in range(1, x.shape[0]):
        if arg <= x[i]:
            return a[i] + b[i] * (arg - x[i - 1]) + c[i] * (arg - x[i - 1]) ** 2 + d[i] * (arg - x[i - 1]) ** 3


def Koshi_Adams(t_, f_0, f_1, u_0):
    t = np.array(t_)
    n = t.shape[0]
    tao = t[1] - t[0]
    u = np.empty((2, n))
    u[:, 0] = np.reshape(u_0, (2))
    f0 = np.array([f_0(t[0], u_0), f_1(t[0], u_0)])
    f1 = np.empty(2)
    # вычислим u_1 с помощью явной схемы Эйлера
    u[:, 1] = u[:, 0] + tao * f0
    # вычислим остальные значения u_k с помощью явной схемы Адамса
    for k in range(1, n - 1):
        f1 = np.array([f_0(t[k], u[:, k]), f_1(t[k], u[:, k])])
        u[:, k + 1] = u[:, k] + (tao / 2) * (3 * f1 - f0)
        f0, f1 = f1, f0
    return u


def z_par():
    global a, b, T
    try:
        (a, b, T) = map(float, ent1.get().split())
        isok.set('Параметры z(t) введены успешно')
    except BaseException:
        isok.set('Ошибка ввода параметров z(t)')


def s_par():
    global c, d
    try:
        (c, d) = map(float, ent2.get().split())
        isok.set('Параметры S(t) введены успешно')
    except BaseException:
        isok.set('Ошибка ввода параметров S(t)')


def ro_par():
    global e, f
    try:
        (e, f) = map(float, ent3.get().split())
        isok.set('Параметры ro(w) введены успешно')
    except BaseException:
        isok.set('Ошибка ввода параметров ro(w)')


def u_par():
    global numy
    try:
        numy = int(ent4.get())
        isok.set('Параметр U(y) введён успешно')
    except BaseException:
        isok.set('Ошибка ввода параметра U(y)')


def k_par():
    global x0, y0, beta
    try:
        (x0, y0, beta) = map(float, ent5.get().split())
        isok.set('Параметры задачи Коши введены успешно')
    except BaseException:
        isok.set('Ошибка ввода параметров задачи Коши')


def comp_z():
    global z, a, b, T
    if a is None:
        isok.set('Ошибка вычисления сетки z(t): отсутствует a')
        return
    if b is None:
        isok.set('Ошибка вычисления сетки z(t): отсутствует b')
        return
    file_name = fd.asksaveasfilename(title='z(t)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    t = np.linspace(0, T, 1000)
    z_t = a * t + b * np.cos(t)
    z = pd.DataFrame({'t': t, 'z(t)': z_t})
    z.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для z(t) вычислена успешно')


def comp_s():
    global S, c, d, T
    if c is None:
        isok.set('Ошибка вычисления сетки S(t): отсутствует c')
        return
    if d is None:
        isok.set('Ошибка вычисления сетки S(t): отсутствует d')
        return
    file_name = fd.asksaveasfilename(title='S(t)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    t = np.linspace(0, T, 1000)
    S_t = c * t + d * np.sin(t)
    S = pd.DataFrame({'t': t, 'S(t)': S_t})
    S.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для S(t) вычислена успешно')


def comp_ro():
    global ro, e, f
    if e is None:
        isok.set('Ошибка вычисления сетки ro(w): отсутствует e')
        return
    if f is None:
        isok.set('Ошибка вычисления сетки ro(w): отсутствует f')
        return
    file_name = fd.asksaveasfilename(title='ro(w)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    w = np.linspace(0, 1, 1000)
    ro_w = e * w * (f - w)
    ro = pd.DataFrame({'w': w, 'ro(w)': ro_w})
    ro.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для ro(w) вычислена успешно')


def U_grid():
    global U, ro, numy, e, f, aU, bU, cU, dU, y, inter_coeffs
    if numy is None:
        isok.set('Ошибка вычисления сетки U(t): отсутствует numy')
        return
    if ro is None:
        isok.set('Ошибка расчёта сетки U(y): отсутствует сетка для ro(w)')
        return
    file_name = fd.asksaveasfilename(title='Коэффиценты интерполяции U(y)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    y = np.linspace(0, 1, numy + 2)
    U_y = np.empty(y.shape[0])
    for i in range(y.shape[0] - 1):
        w = np.linspace(y[i], 1, 1000)
        ro_w = e * w * (f - w)
        U_y[i] = integral_trap(w, ro_w, w.shape[0] - 1)
    U_y[y.shape[0] - 1] = 0
    U = pd.DataFrame({'y': y, 'U(y)': U_y})
    aU, bU, cU, dU = cubic_spline(y, U_y)
    inter_coeffs = pd.DataFrame({'a' : aU, 'b' : bU, 'c' : cU, 'd' : dU})
    inter_coeffs.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка и интерполяционные коэффиценты для U(y) вычислены успешно')


def func1(t, u):
    tmpu = u[1]
    if tmpu < 0:
        tmpu = 0
    elif tmpu > 1:
        tmpu = 1
    tmpz = a - b * np.sin(t)
    return tmpz * interpolated_func(tmpu, y, aU, bU, cU, dU)


def func2(t, u):
    global beta, a, b
    return beta * (u[0] - (a * t + b * np.cos(t)))


def Koshi():
    global U, inter_coeffs, S, z, ro, koshi, x0, y0, beta, T
    if ro is None:
        isok.set('Ошибка решения задачи Коши: отсутствует сетка для ro(w)')
        return
    if S is None:
        isok.set('Ошибка решения задачи Коши: отсутствует сетка для S(t)')
        return
    if z is None:
        isok.set('Ошибка решения задачи Коши: отсутствует сетка для z(t)')
        return
    if U is None:
        isok.set('Ошибка решения задачи Коши: отсутствует сетка для U(y)')
        return
    if inter_coeffs is None:
        isok.set('Ошибка решения задачи Коши: отсутствуют коэффиценты интерполяции для U(y)')
        return
    if x0 is None:
        isok.set('Ошибка решения задачи Коши: отсутствует x0')
        return
    if y0 is None:
        isok.set('Ошибка решения задачи Коши: отсутствует y0')
        return
    if beta is None:
        isok.set('Ошибка решения задачи Коши: отсутствует beta')
        return
    if T is None:
        isok.set('Ошибка решения задачи Коши: отсутствует T')
        return
    file_name = fd.asksaveasfilename(title='Koshi', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    koshi_ans = Koshi_Adams(np.linspace(0, T, 1000), func1, func2, np.array([x0, y0]))
    koshi = pd.DataFrame({'t' : np.linspace(0, T, 1000), 'x(t)' : koshi_ans[0], 'y(t)' : koshi_ans[1]})
    koshi.to_csv(file_name + '.csv', index=False)
    isok.set('Задача Коши решена успешно')


def plots():
    if ro is None:
        isok.set('Ошибка построения графиков: отсутствует сетка для ro(w)')
        return
    if S is None:
        isok.set('Ошибка построения графиков: отсутствует сетка для S(t)')
        return
    if koshi is None:
        isok.set('Ошибка построения графиков: отсутствует решение задачи Коши')
        return
    plt.clf()
    plt.subplot(141)
    plt.plot(ro.iloc[:,0], ro.iloc[:,1])
    plt.title('График ro(w)')
    plt.xlabel('w')
    plt.ylabel('ro(w)')
    plt.subplot(142)
    plt.plot(koshi['t'], koshi['x(t)'], label='x(t)')
    plt.plot(S['t'], S['S(t)'], label='S(t)')
    diffSx = np.array(koshi['x(t)']) - np.array(S['S(t)'])
    plt.plot(S['t'], diffSx, label='x(t)-S(t)')
    plt.title('Графики x(t), S(t), x(t)-S(t)')
    plt.xlabel('t')
    plt.legend()
    plt.subplot(143)
    plt.plot(z['t'], z['z(t)'], label='z(t)')
    plt.title('График z(t)')
    plt.xlabel('t')
    plt.ylabel('z(t)')
    plt.subplot(144)
    plt.plot(koshi['t'], koshi['y(t)'], label='y(t)')
    plt.title('График y(t)')
    plt.xlabel('t')
    plt.ylabel('y(t)')
    # if betas is not None:
    #    plt.subplot()
    fig.canvas.draw()


def crit_plots():
    None


def main():
    for arg in sys.argv[1:]:
        if arg == 'test1':
            ent1.insert(0, '1 2')
            ent2.insert(0, '3 4')
            ent3.insert(0, '5')
            z_par()
            s_par()
            ro_par()
            comp_z()
            comp_s()
            comp_ro()
            print('test1 - ok')
        elif arg == 'test2':
            ent3.insert(0, '5')
            ent4.insert(0, '6')
            ro_par()
            comp_ro()
            u_par()
            U_grid()
            print('test2 - ok')
        elif arg == 'test3':
            ent1.insert(0, '1 2')
            ent2.insert(0, '3 4')
            ent3.insert(0, '5')
            ent4.insert(0, '6')
            z_par()
            s_par()
            ro_par()
            u_par()
            comp_z()
            comp_s()
            comp_ro()
            U_grid()
            ent5.insert(0, '7 8 9 10')
            k_par()
            Koshi()
            print('test3 - ok')
        else:
            return 1
    ent1.delete(0, END)
    ent2.delete(0, END)
    ent3.delete(0, END)
    ent4.delete(0, END)
    ent5.delete(0, END)
    myApp()


def myApp():
    root.mainloop()

warnings.filterwarnings('ignore')

matplotlib.use('TkAgg')

root = Tk()

a = None
b = None
c = None
d = None
e = None
f = None
z = None
S = None
ro = None
U = None
y = None
inter_coeffs = None
koshi = None
numy = None
x0 = None
y0 = None
beta = None
betas = None
T = None
aU = None
bU = None
cU = None
dU = None

fig = plt.figure(1, figsize=(20, 6))
canvas = FigureCanvasTkAgg(fig, master=root)
plot_widget = canvas.get_tk_widget()

isok = StringVar()
isok.set('')
lab1 = Label(textvariable=isok, bg="white", fg="red")
lab1.place(relx=0.4, rely=0.05)

lab1 = Label(text='Введите параметры a, b функции z(t)=at+bcos(t),\nа также T (через пробел):', bg="white", fg="black")
lab1.place(relx=0.1, rely=0.1)

ent1 = Entry()
ent1.place(relx=0.4, rely=0.1)

b11 = Button(text="Ввести параметры", command=z_par)
b11.place(relx=0.6, rely=0.1)

b1 = Button(text="Вычислить сетку для z(t)", command=comp_z)
b1.place(relx=0.7, rely=0.1)

lab2 = Label(text='Введите параметры c, d функции S(t)=ct+dsin(t)\n(через пробел):', bg="white", fg="black")
lab2.place(relx=0.1, rely=0.15)

ent2 = Entry()
ent2.place(relx=0.4, rely=0.15)

b22 = Button(text="Ввести параметры", command=s_par)
b22.place(relx=0.6, rely=0.15)

b2 = Button(text="Вычислить сетку для S(t)", command=comp_s)
b2.place(relx=0.7, rely=0.15)

lab3 = Label(text='Введите параметры e, f функции ro(t)=ew(f-w):', bg="white", fg="black")
lab3.place(relx=0.1, rely=0.2)

ent3 = Entry()
ent3.place(relx=0.4, rely=0.2)

b33 = Button(text="Ввести параметры", command=ro_par)
b33.place(relx=0.6, rely=0.2)

b3 = Button(text="Вычислить сетку для ro(w)", command=comp_ro)
b3.place(relx=0.7, rely=0.2)

lab4 = Label(text='Введите количество точек разбиения numy:', bg="white", fg="black")
lab4.place(relx=0.1, rely=0.25)

ent4 = Entry()
ent4.place(relx=0.4, rely=0.25)

b44 = Button(text="Ввести параметры", command=u_par)
b44.place(relx=0.6, rely=0.25)

b4 = Button(text="Вычислить сетку для U(y)", command=U_grid)
b4.place(relx=0.7, rely=0.25)

lab5 = Label(text='Введите параметры x0, y0, beta для задачи Коши\n(через пробел):', bg="white", fg="black")
lab5.place(relx=0.1, rely=0.3)

ent5 = Entry()
ent5.place(relx=0.4, rely=0.3)

b55 = Button(text="Ввести параметры", command=k_par)
b55.place(relx=0.6, rely=0.3)

b5 = Button(text="Решить задачу Коши для вычисленных сеток", command=Koshi)
b5.place(relx=0.7, rely=0.3)

b6 = Button(text="Показать графики", command=plots)
b6.place(relx=0.3, rely=0.35)

b7 = Button(text="Показать графики критериев", command=crit_plots)
b7.place(relx=0.5, rely=0.35)

plot_widget.place(relx=0, rely=0.4)

if __name__ == "__main__":
    main()
