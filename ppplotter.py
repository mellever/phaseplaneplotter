from tkinter import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from matplotlib.backend_bases import MouseEvent
from matplotlib.figure import Figure
from scipy.integrate import solve_ivp

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk


counter = 0
sol_p = 0 
sol_n = 0
iv = [0,0]
t_span = 40
t_span_n = [t_span/2, 0] #time interval over which the solutions are integrated backwards in time
t_span_p = [0, t_span/2] #time interval over which the solutions are integrated forwards in time
t_p = np.linspace(t_span_p[0], t_span_p[1], 5000)
t_n = np.linspace(t_span_n[0], t_span_n[1], 100)

def f(t, z):
    u, v = z
    return [eval(u_in.get()), eval(v_in.get())]

def on_click(event):
    if event.inaxes is not None:
        print(event.xdata, event.ydata)
    else:
        print('Clicked ouside axes bounds but inside plot window')

def callback(event):
    if event.inaxes is not None:
        iv = [event.xdata, event.ydata]
        sol_p = solve_ivp(f, t_span_p, iv, t_eval = t_p, dense_output=True)
        sol_n = solve_ivp(f, t_span_n, iv, t_eval = t_n, dense_output=True)
        plot(sol_p, sol_n, iv)
    else: print('Clicked ouside axes bounds but inside plot window')

def normalize(a,b):
    a = a / np.sqrt(a**2+b**2)
    b = b / np.sqrt(a**2+b**2)
    return a,b

def plot(sol_p, sol_n, iv):
    norm = True
    global counter 
    counter = 1
    plotter(norm, sol_p, sol_n, iv)


def switch():
    global counter
    counter += 1
    if counter%2 == 1:
        norm = True
        plotter(norm, sol_p, sol_n, iv)
    else:
        norm = False
        plotter(norm, sol_p, sol_n, iv)

def ppsolver():
    u_min = float(umin.get()) #minimal value of u to plot
    u_max = float(umax.get())  #maximal value of u to plot
    
    v_min = float(vmin.get())  #minimal value of v to plot
    v_max = float(vmax.get())  #maximal value of v to plot

    x = np.linspace(u_min, u_max, int(np.sqrt(eval(str(arrows.get())))))
    y = np.linspace(v_min, v_max, int(np.sqrt(eval(str(arrows.get())))))

    X, Y = np.meshgrid(x, y) #make matrix from vectors

    t = 0 #starting at t=0

    a, b = np.zeros(X.shape), np.zeros(Y.shape)

    NI, NJ = X.shape

    for i in range(NI):
        for j in range(NJ):
            x = X[i, j]
            y = Y[i, j]
            u = x
            v = y
            yprime = [eval(str(u_in.get())), eval(str(v_in.get()))]
            a[i,j] = yprime[0]
            b[i,j] = yprime[1]

    return X,Y,a,b

def plotter(norm, sol_p, sol_n, iv): 
    fig = Figure(figsize=(5,5), dpi=100)
    X,Y,u,v = ppsolver()
    if norm: u,v = normalize(u,v)
    p = fig.add_subplot(111)
    #p.title.set_text(r'$input$')
    p.set_xlabel('u')
    p.set_ylabel('v')
    p.set_xlim(eval(umin.get()),eval(umax.get()))
    p.set_ylim(eval(vmin.get()),eval(vmax.get()))
    p.quiver(X,Y,u,v, color='black', scale = eval('50/'+str(scale.get())))
    if (sol_p != 0):
        print('1')
        z_p = sol_p.sol(t_p)
        z_n = sol_n.sol(t_n)
        p.scatter(iv[0], iv[1], marker='o', c='blue', label = 'iv = ('+str(round(iv[0], 2)) +', ' +str(round(iv[1],2))+')')
        p.plot(z_p[0], z_p[1], color = 'green', label = r'$t_+$')
        p.plot(z_n[0], z_n[1], color = 'red', label = r'$t_-$')
        p.legend()
    p.plot()

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=3, column=0,rowspan=5, columnspan=5)
    canvas.mpl_connect('button_press_event',callback)
    toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=True)
    toolbar.update()
    canvas._tkcanvas.grid(row=3, column=0,rowspan=5, columnspan=5)
    fig, ax = plt.subplots()
    fig.canvas.callbacks.connect('button_press_event', on_click)
    plt.show()


#create screen
root = Tk()
root.title='Phase Plane Plotter'

#Create label widget
uLabel = Label(root, text='du/dt = ')
vLabel = Label(root, text='dv/dt = ')

#Putting label onto screen
uLabel.grid(row=0, column=0)
vLabel.grid(row=1, column=0)

u_in = Entry(root, borderwidth=1)
u_in.grid(row=0, column=1)
u_in.insert(0, 'v - 0.5*u')

v_in = Entry(root, borderwidth=1)
v_in.grid(row=1, column=1)
v_in.insert(0, 'np.sin(u)')

plotBotton = Button(root, text='plot', padx=10, pady=10, command= lambda: plot(sol_p, sol_n, iv))
plotBotton.grid(row=0, column=2)

umin = Entry(root)
umin.grid(row=2, column=0)
umin.insert(0, -5)

umax = Entry(root)
umax.grid(row=2, column=1)
umax.insert(0, 5)

vmin = Entry(root)
vmin.grid(row=2, column=2)
vmin.insert(0, -5)

vmax = Entry(root)
vmax.grid(row=2, column=3)
vmax.insert(0, 5)

normBotton = Button(root, text='Variable', command=switch)
normBotton.grid(row=0, column=3)

arrows = Entry(root)
arrows.grid(row=1, column=2)
arrows.insert(0, '1000')

scale = Entry(root)
scale.grid(row=1, column=3)
scale.insert(0, '1')


root.mainloop()
