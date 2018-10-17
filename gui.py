from tkinter import *
from tkinter import filedialog as fd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def read_z():
    global z
    file_name = fd.askopenfilename(filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    z = pd.read_csv(file_name, engine='python')


def read_s():
    global S
    file_name = fd.askopenfilename(filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    S = pd.read_csv(file_name, engine='python')

def read_ro():
    global ro
    file_name = fd.askopenfilename(filetypes=[("csv files", "*.csv")])
    if file_name == '':
        return
    ro = pd.read_csv(file_name, engine='python')


def U_grid():
    global U, ro
    if ro is None:
        isok.set('Ошибка расчёта сетки U(y): отсутствует сетка для ro(w)')
        return
    file_name = fd.asksaveasfilename(filetypes=[("csv files", "*.csv")])
    U = ro + ro
    U.to_csv(file_name + '.csv', index=False)
    isok.set('Сетка для U(y) вычислена успешно')


def Koshi():
    global U, S, z, ro, koshi
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
    file_name = fd.asksaveasfilename(filetypes=[("csv files", "*.csv")])
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


matplotlib.use('TkAgg')

root = Tk()

z = None
S = None
ro = None
U = None
koshi = None

fig = plt.figure(1, figsize=(14, 6))
canvas = FigureCanvasTkAgg(fig, master=root)
plot_widget = canvas.get_tk_widget()

isok = StringVar()
isok.set('')
lab1 = Label(textvariable=isok, bg="white", fg="red")
lab1.pack()

b1 = Button(text="Загрузить сетку для Z(t)", command=read_z)
b1.pack()

b2 = Button(text="Загрузить сетку для S(t)", command=read_s)
b2.pack()

b3 = Button(text="Загрузить сетку для ro(w)", command=read_ro)
b3.pack()

b4 = Button(text="Вычислить сетку для U(y)", command=U_grid)
b4.pack()

b5 = Button(text="Решить задачу Коши для вычисленных сеток", command=Koshi)
b5.pack()

b6 = Button(text="Показать графики", command=ro_U_plots)
b6.pack()

plot_widget.pack()

root.mainloop()