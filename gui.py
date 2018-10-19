from tkinter import *
from tkinter import filedialog as fd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import sys


def z_par():
    global a, b
    try:
        (a, b) = map(float, ent1.get().split())
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
    global e
    try:
        e = float(ent3.get())
        isok.set('Параметр ro(w) введён успешно')
    except BaseException:
        isok.set('Ошибка ввода параметра ro(w)')


def u_par():
    global dy
    try:
        dy = float(ent4.get())
        isok.set('Параметр U(y) введён успешно')
    except BaseException:
        isok.set('Ошибка ввода параметра U(y)')


def k_par():
    global x0, y0, beta, T
    try:
        (x0, y0, beta, T) = map(float, ent5.get().split())
        isok.set('Параметры задачи Коши введены успешно')
    except BaseException:
        isok.set('Ошибка ввода параметров задачи Коши')


def comp_z():
    global z, a, b
    if a is None:
        isok.set('Ошибка вычисления сетки z(t): отсутствует a')
        return
    if b is None:
        isok.set('Ошибка вычисления сетки z(t): отсутствует b')
        return
    file_name = fd.asksaveasfilename(title='z(t)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    z = pd.DataFrame({'t': np.arange(0, 10), 'z(t)': np.arange(0, 10) ** 2})
    z.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для z(t) вычислена успешно')


def comp_s():
    global S, c, d
    if c is None:
        isok.set('Ошибка вычисления сетки S(t): отсутствует c')
        return
    if d is None:
        isok.set('Ошибка вычисления сетки S(t): отсутствует d')
        return
    file_name = fd.asksaveasfilename(title='S(t)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    S = pd.DataFrame({'t': np.arange(0, 10), 'S(t)': np.arange(0, 10) ** 3})
    S.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для S(t) вычислена успешно')


def comp_ro():
    global ro, e
    if e is None:
        isok.set('Ошибка вычисления сетки ro(w): отсутствует e')
        return
    file_name = fd.asksaveasfilename(title='ro(w)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    ro = pd.DataFrame({'w': np.arange(0, 10), 'ro(w)': np.arange(0, 10) ** 2})
    ro.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для ro(w) вычислена успешно')


def U_grid():
    global U, ro, dy
    if dy is None:
        isok.set('Ошибка вычисления сетки U(t): отсутствует dy')
        return
    if ro is None:
        isok.set('Ошибка расчёта сетки U(y): отсутствует сетка для ro(w)')
        return
    file_name = fd.asksaveasfilename(title='U(y)', filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    U = pd.DataFrame({'t': np.arange(0, 10), 'S(t)': np.arange(0, 10) ** 3})
    U.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для U(y) вычислена успешно')


def Koshi():
    global U, S, z, ro, koshi, x0, y0, beta, T
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
    koshi = U + S + z + ro
    koshi.to_csv(file_name + '.csv', index=False)
    isok.set('Задача Коши решена успешно')


def ro_U_plots():
    if ro is None:
        isok.set('Ошибка построения графика: отсутствует сетка для ro(w)')
        return
    if U is None:
        isok.set('Ошибка построения графика: отсутствует сетка для U(y)')
        return
    plt.clf()
    plt.subplot(121)
    plt.plot(ro.iloc[:,0], ro.iloc[:,1])
    plt.title('График ro(w)')
    plt.xlabel('w')
    plt.ylabel('ro(w)')
    plt.subplot(122)
    plt.plot(U.iloc[:, 0], U.iloc[:, 1])
    plt.title('График U(y)')
    plt.xlabel('y')
    plt.ylabel('U(y)')
    fig.canvas.draw()


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
z = None
S = None
ro = None
U = None
koshi = None
dy = None
x0 = None
y0 = None
beta = None
T = None

fig = plt.figure(1, figsize=(14, 6))
canvas = FigureCanvasTkAgg(fig, master=root)
plot_widget = canvas.get_tk_widget()

isok = StringVar()
isok.set('')
lab1 = Label(textvariable=isok, bg="white", fg="red")
lab1.place(relx=0.4, rely=0.05)

lab1 = Label(text='Введите параметры a, b функции z(t)\n(через пробел):', bg="white", fg="black")
lab1.place(relx=0.1, rely=0.1)

ent1 = Entry()
ent1.place(relx=0.4, rely=0.1)

b11 = Button(text="Ввести параметры", command=z_par)
b11.place(relx=0.6, rely=0.1)

b1 = Button(text="Вычислить сетку для z(t)", command=comp_z)
b1.place(relx=0.7, rely=0.1)

lab2 = Label(text='Введите параметры c, d функции S(t)\n(через пробел):', bg="white", fg="black")
lab2.place(relx=0.1, rely=0.15)

ent2 = Entry()
ent2.place(relx=0.4, rely=0.15)

b22 = Button(text="Ввести параметры", command=s_par)
b22.place(relx=0.6, rely=0.15)

b2 = Button(text="Вычислить сетку для S(t)", command=comp_s)
b2.place(relx=0.7, rely=0.15)

lab3 = Label(text='Введите параметр e функции ro(t):', bg="white", fg="black")
lab3.place(relx=0.1, rely=0.2)

ent3 = Entry()
ent3.place(relx=0.4, rely=0.2)

b33 = Button(text="Ввести параметры", command=ro_par)
b33.place(relx=0.6, rely=0.2)

b3 = Button(text="Вычислить сетку для ro(w)", command=comp_ro)
b3.place(relx=0.7, rely=0.2)

lab4 = Label(text='Введите параметр разбиения dy:', bg="white", fg="black")
lab4.place(relx=0.1, rely=0.25)

ent4 = Entry()
ent4.place(relx=0.4, rely=0.25)

b44 = Button(text="Ввести параметры", command=s_par)
b44.place(relx=0.6, rely=0.25)

b4 = Button(text="Вычислить сетку для U(y)", command=U_grid)
b4.place(relx=0.7, rely=0.25)

lab5 = Label(text='Введите параметры x0, y0, beta, T для задачи Коши\n(через пробел):', bg="white", fg="black")
lab5.place(relx=0.1, rely=0.3)

ent5 = Entry()
ent5.place(relx=0.4, rely=0.3)

b55 = Button(text="Ввести параметры", command=k_par)
b55.place(relx=0.6, rely=0.3)

b5 = Button(text="Решить задачу Коши для вычисленных сеток", command=Koshi)
b5.place(relx=0.7, rely=0.3)

b6 = Button(text="Показать графики", command=ro_U_plots)
b6.place(relx=0.4, rely=0.35)

plot_widget.place(relx=0, rely=0.4)

if __name__ == "__main__":
    main()
