#Import packages
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import *

import numpy as np

from sympy import sympify, lambdify, solve, solveset, Interval, Eq, latex, sin, cos, Matrix, re
from sympy.abc import u, v

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure 


# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        #Set window title
        self.setWindowTitle("Phase Plane Plotter")

        #Create labels
        LabelF1 = QLabel("du/dt")
        LabelF2 = QLabel("dv/dt")
        LabelDensity = QLabel("Density")
        LabelXmin = QLabel("Xmin")
        LabelXmax = QLabel("Xmax")
        LabelYmin = QLabel("Ymin")
        LabelYmax = QLabel("Ymax")
        LabelNullclines = QLabel("Nullclines")
        LabelEquilibria = QLabel("Equilibria")

        #Create checkbox for nullclines
        self.NullclinesCheckbox = QCheckBox()
        self.EquilibriaCheckbox = QCheckBox()


        #Create button for plotting
        PlotButton = QPushButton("Plot")
        PlotButton.setCheckable(True)
        PlotButton.clicked.connect(self.plot_button)
        
        #Create input for F1
        self.InputF1 = QLineEdit()
        self.InputF1.setText("v-0.5*u")

        #Create input for F2
        self.InputF2 = QLineEdit()
        self.InputF2.setText("sin(u)")

        #Create input for Xmin, Xmax, Ymin, Ymax
        self.InputXmin = QLineEdit()
        self.InputXmin.setText("-5")
        self.InputXmax = QLineEdit()
        self.InputXmax.setText("5")
        self.InputYmin = QLineEdit()
        self.InputYmin.setText("-5")
        self.InputYmax = QLineEdit()
        self.InputYmax.setText("5")

        #Create canvas to display plot
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel("u")
        self.ax.set_ylabel("v")
        plt.tight_layout()
        self.canvas = FigureCanvasQTAgg(self.fig)
       
        #Create matplotlib toolbar
        Toolbar = NavigationToolbar2QT(self.canvas, self)

        #Create slider for density lines
        self.Slider = QSlider(Qt.Orientation.Horizontal)
        self.Slider.setMinimum(1)
        self.Slider.setMaximum(5)
        self.Slider.setValue(2)

        #Set layout
        layout = QGridLayout()

        #Add widgets
        layout.addWidget(LabelF1, 0, 0)
        layout.addWidget(self.InputF1, 0, 1)

        layout.addWidget(LabelF2, 1, 0)
        layout.addWidget(self.InputF2, 1, 1)
        
        layout.addWidget(LabelXmin, 2, 0)
        layout.addWidget(self.InputXmin, 2, 1)

        layout.addWidget(LabelXmax, 3, 0)
        layout.addWidget(self.InputXmax, 3, 1)

        layout.addWidget(LabelYmin, 4, 0)
        layout.addWidget(self.InputYmin, 4, 1)

        layout.addWidget(LabelYmax, 5, 0)
        layout.addWidget(self.InputYmax, 5, 1)

        layout.addWidget(LabelDensity, 6, 0)
        layout.addWidget(self.Slider, 6, 1)

        layout.addWidget(LabelEquilibria, 7, 0)
        layout.addWidget(self.EquilibriaCheckbox, 7, 1)
        
        layout.addWidget(LabelNullclines, 8, 0)
        layout.addWidget(self.NullclinesCheckbox, 8, 1)
             
        layout.addWidget(PlotButton, 10, 0, 1, 2)
        layout.addWidget(self.canvas, 0, 2, 10, 1)
        layout.addWidget(Toolbar, 10, 2)

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

    def plot_button(self):
        Xmin = float(self.InputXmin.text())
        Xmax = float(self.InputXmax.text())
        Ymin = float(self.InputYmin.text())
        Ymax = float(self.InputYmax.text())

        F1 = self.InputF1.text()
        F2 = self.InputF2.text()

        StepSize = 0.2

        x = np.arange(Xmin, Xmax, StepSize)
        y = np.arange(Ymin, Ymax, StepSize)

        X, Y = np.meshgrid(x, y)

        F1_sympy = sympify(F1)
        F2_sympy = sympify(F2)


        F1 = lambdify([u,v], F1_sympy)
        F2 = lambdify([u,v], F2_sympy)


        U = F1(X, Y)
        V = F2(X, Y)
        
        self.ax.clear()
        self.ax.set_xlabel("u")
        self.ax.set_ylabel("v")
        self.ax.set_xlim(Xmin, Xmax)
        self.ax.set_ylim(Ymin, Ymax)
        plt.tight_layout()
        self.ax.streamplot(X, Y, U, V, density=self.Slider.value(), color="tab:gray")


        plotColors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"] #Add more colors to be sure


        if self.EquilibriaCheckbox.isChecked():
            SolutionSet = []

            #Compute Jacobian matrix
            F = Matrix([F1_sympy, F2_sympy])
            J = F.jacobian(Matrix([u,v]))

            if F1_sympy.has(cos) or F1_sympy.has(sin):
                if F2_sympy.has(cos) or F2_sympy.has(sin):
                    print("double trig functions")
                else:
                    F1_symbols = F1_sympy.free_symbols
                    if len(F1_symbols) == 2:
                        sol = solveset(F1_sympy, (u,v))
                        print("How to handle this?")
                    elif len(F1_symbols) == 1:
                        if u in F1_symbols:
                            sol1 = solveset(F1_sympy, u, domain=Interval(Xmin, Xmax))
                            for s in sol1:
                                sol = solve(F2_sympy.subs(u, s))
                                for s2 in sol:
                                    SolutionSet.append([s, s2])
                        else:
                            sol1 = solveset(F1_sympy, v, domain=Interval(Ymin, Ymax))
                            for s in sol1:
                                sol = solve(F2_sympy.subs(v, s))
                                for s2 in sol:
                                    SolutionSet.append([s2, s])
                    else: print("Not interesting")
            elif F2_sympy.has(cos) or F2_sympy.has(sin):
                if F1_sympy.has(cos) or F1_sympy.has(sin):
                    print("double trig functions, but now handled by the above")
                else:
                    F2_symbols = F2_sympy.free_symbols
                    if len(F2_symbols) == 2:
                        sol = solve(F2_sympy, [u,v])
                        print(sol)
                        print("How to handle this?")
                    elif len(F2_symbols) == 1:
                        if u in F2_symbols:
                            sol1 = solveset(F2_sympy, u, domain=Interval(Xmin, Xmax))
                            for s in sol1:
                                sol = solve(F1_sympy.subs(u, s))
                                for s2 in sol:
                                    SolutionSet.append([s, s2])
                        else:
                            sol1 = solveset(F2_sympy, v, domain=Interval(Ymin, Ymax))
                            for s in sol1:
                                sol = solve(F1_sympy.subs(v, s))
                                for s2 in sol:
                                    SolutionSet.append([s2, s])
                    else: print("Not interesting")
            else:
                sol = solve([F1_sympy, F2_sympy], [u,v])
            
            for eq in SolutionSet:
                J_eq = J.subs([(u, eq[0]), (v, eq[1])])
                ev_J = list(J_eq.eigenvals().keys())
                re_ev_J = np.array([re(i) for i in ev_J])
                if len(re_ev_J[re_ev_J>0])>0:
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=50, color="black", zorder=10)
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=20, color="white", zorder=10)
                else:
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=50, color="black", zorder=10)


        if self.NullclinesCheckbox.isChecked():
            solCounter = 0
            for F in [F1_sympy, F2_sympy]:
                F_symbols = F.free_symbols
                if F.has(cos) or F.has(sin):
                    if len(F_symbols) == 2:
                        F_sol = solveset(F, v)
                        for sol in F_sol:
                            F_sol_num = lambdify(u, sol, "numpy")
                            self.ax.plot(x, F_sol_num(x), c=plotColors[solCounter], label=r"$v = " + latex(sol) + r"$")
                            solCounter+=1
                    elif len(F_symbols) == 1:
                        if u in F_symbols:
                            F_sol = solveset(F, u, domain=Interval(Xmin, Xmax))
                            for sol in F_sol:
                                self.ax.vlines(sol, Ymin, Ymax, color=plotColors[solCounter], label=r"$u = " + latex(sol) + r"$")
                                solCounter+=1
                        else:
                            F_sol = solveset(F, v, domain=Interval(Ymin, Ymax))
                            for sol in F_sol:
                                self.ax.hlines(sol, Xmin, Xmax, color=plotColors[solCounter], label=r"$v = " + latex(sol) + r"$")
                                solCounter+=1
                    else: print("Not interesting")
                else:
                    if len(F_symbols) == 2:
                        F_sol = solve(F, v)
                        for sol in F_sol:
                            F_sol_num = lambdify(u, sol, "numpy")
                            self.ax.plot(x, F_sol_num(x), c=plotColors[solCounter], label=r"$v = " + latex(sol) + r"$")
                            solCounter+=1
                    elif len(F_symbols) == 1:
                        if u in F_symbols:
                            F_sol = solve(F, u, domain=Interval(Xmin, Xmax))
                            for sol in F_sol:
                                self.ax.vlines(sol, Ymin, Ymax, color=plotColors[solCounter], label=r"$u = " + latex(sol) + r"$")
                                solCounter+=1
                        else:
                            F_sol = solve(F, v, domain=Interval(Ymin, Ymax))
                            for sol in F_sol:
                                self.ax.hlines(sol, Xmin, Xmax, color=plotColors[solCounter], label=r"$v = " + latex(sol) + r"$")
                                solCounter+=1
                    else: print("Not interesting")
            self.ax.legend()
        
        self.canvas.draw()


app = QApplication([])

window = MainWindow()
window.show()

# Start the event loop.
app.exec()