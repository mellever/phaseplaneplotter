#Import packages
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import *
from PyQt6.QtGui import QPixmap

import numpy as np
import io

from sympy import sympify, lambdify, solve, solveset, Interval, Eq, latex, sin, cos, Matrix, re, pi, nsimplify
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
        LabelXmin = QLabel("Umin")
        LabelXmax = QLabel("Umax")
        LabelYmin = QLabel("Vmin")
        LabelYmax = QLabel("Vmax")
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
        #self.InputF1.setText("-u+3*v")
        #self.InputF1.setText("u*(-2*u-v+180)")
        #self.InputF1.setText("cos(u)")
        self.InputF1.setText("v-0.5*u")
        #self.InputF1.setText("sin(u+v)")


        #Create input for F2
        self.InputF2 = QLineEdit()
        self.InputF2.setText("-3*v")
        self.InputF2.setText("v*(-u-2*v+120)")
        #self.InputF2.setText("cos(v)")
        self.InputF2.setText("sin(u)")
        #self.InputF2.setText("sin(u+v)")
        #self.InputF2.setText("v-0.5*u")

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
        self.fig, self.ax = plt.subplots(figsize=(5,5))
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

        #Create place to render latex equations
        self.latexPanel = QLabel(self)
        self.latex_code = []
        self.render_latex()

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
        layout.addWidget(self.canvas, 0, 2, 10, 5)
        layout.addWidget(Toolbar, 10, 2)

        layout.addWidget(self.latexPanel, 0, 7, 10, 3)

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

        self.plotColors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:cyan"] #Add more colors to be sure

    def render_latex(self):
        # Use matplotlib to create an image from the LaTeX string
        fig, ax = plt.subplots(figsize=(3,5))
        yline = 1
        i = 0
        for line in self.latex_code:
            if line == "Equilibria:" or line == "Nullclines:":
                yline -= 0.05
                ax.text(0, yline, line, fontsize=12, va="center")
            else:
                if "Unstable" in line or "Stable" in line:
                    line_split = line.split("_")
                    ax.text(0, yline, f"${line_split[0]}$", fontsize=12, va="center")
                    ax.text(0.5, yline, line_split[1], fontsize=12, va="center")

                else: 
                    ax.text(0.1, yline, f"${line}$", fontsize=12, va="center")
                    ax.text(0, yline-0.02, "\u2022", c=self.plotColors[i], fontsize=16)
                    i+=1
            yline -= 0.05
        ax.axis('off')
        plt.tight_layout()

        # Save the figure to a buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        plt.close(fig)

        # Load the image into QPixmap and display it in the QLabel
        pixmap = QPixmap()
        pixmap.loadFromData(buf.getvalue())
        self.latexPanel.setPixmap(pixmap)

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

        if self.EquilibriaCheckbox.isChecked():
            SolutionSet = []

            #Compute Jacobian matrix
            F = Matrix([F1_sympy, F2_sympy])
            J = F.jacobian(Matrix([u,v]))

            if F1_sympy.has(cos) or F1_sympy.has(sin):
                if F2_sympy.has(cos) or F2_sympy.has(sin):
                    print("Warning, does not find all equilibria!")
                    sol = solve([F1_sympy, F2_sympy])
                    if type(sol) is list:
                        for s in sol:
                            SolutionSet.append([s[u], s[v]])
                    else:
                        SolutionSet.append([sol[u], sol[v]])
                else:
                    F1_symbols = F1_sympy.free_symbols
                    if len(F1_symbols) == 2:
                        sol = solve(F1_sympy)
                        for s in sol:
                            if u in s:
                                sol1 = solve(F2_sympy.subs(u, s[u]), v)[0]
                                SolutionSet.append([s[u].subs(v, sol1), sol1])
                            else:
                                sol1 = solve(F2_sympy.subs(v, s[v]))[0]
                                SolutionSet.append([sol1, s[v].subs(u, sol1)])
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
                F2_symbols = F2_sympy.free_symbols
                if len(F2_symbols) == 2:
                    sol = solve(F2_sympy)
                    for s in sol:
                        if u in s:
                            sol1 = solve(F1_sympy.subs(u, s[u]), v)[0]
                            SolutionSet.append([s[u].subs(v, sol1), sol1])
                        else:
                            sol1 = solve(F1_sympy.subs(v, s[v]))[0]
                            SolutionSet.append([sol1, s[v].subs(u, sol1)])

                elif len(F2_symbols) == 1:
                    if u in F2_symbols:
                        sol1 = solveset(F2_sympy, u, domain=Interval(Xmin, Xmax))
                        for s in sol1:
                            sol = solve(F1_sympy.subs(u, s), v)
                            for s2 in sol:
                                SolutionSet.append([s, nsimplify(s2, [pi])])
                    else:
                        sol1 = solveset(F2_sympy, v, domain=Interval(Ymin, Ymax))
                        for s in sol1:
                            sol = solve(F1_sympy.subs(v, s))
                            for s2 in sol:
                                SolutionSet.append([nsimplify(s2, [pi]), s])
                else: print("Not interesting")
            else:
                sol = solve([F1_sympy, F2_sympy])
                if type(sol) is list:
                    for s in sol:
                        SolutionSet.append([s[u], s[v]])
                else:
                    SolutionSet.append([sol[u], sol[v]])

            if "Equilibria:" in self.latex_code: lableGen = False
            else: lableGen = True
            if lableGen: self.latex_code.append("Equilibria:")


            for eq in SolutionSet:
                J_eq = J.subs([(u, eq[0]), (v, eq[1])])
                ev_J = list(J_eq.eigenvals().keys())
                re_ev_J = np.array([re(i) for i in ev_J])
                if len(re_ev_J[re_ev_J>0])>0:
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=50, color="black", zorder=10)
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=20, color="white", zorder=10)
                    if lableGen: self.latex_code.append("("+latex(eq[0])+","+latex(eq[1])+")_Unstable")
                else:
                    self.ax.scatter(eq[0], eq[1], marker = 'o', s=50, color="black", zorder=10)
                    if lableGen: self.latex_code.append("("+latex(eq[0])+","+latex(eq[1])+")_Stable")

        if self.NullclinesCheckbox.isChecked():
            if "Nullclines:" in self.latex_code: lableGen = False
            else: lableGen = True
            if lableGen: self.latex_code.append("Nullclines:")
            solCounter = 0
            for F in [F1_sympy, F2_sympy]:
                F_symbols = F.free_symbols
                if F.has(cos) or F.has(sin):
                    if len(F_symbols) == 2:
                        print("Warning, does not find all nullclines")
                        F_sol = solve(F, v)
                        for sol in F_sol:
                            F_sol_num = lambdify(u, sol, "numpy")
                            self.ax.plot(x, F_sol_num(x), c=self.plotColors[solCounter])
                            if lableGen: self.latex_code.append("v = " + latex(sol))
                            solCounter+=1
                    elif len(F_symbols) == 1:
                        if u in F_symbols:
                            F_sol = solveset(F, u, domain=Interval(Xmin, Xmax))
                            for sol in F_sol:
                                self.ax.vlines(sol, Ymin, Ymax, color=self.plotColors[solCounter])
                                if lableGen: self.latex_code.append("u = " + latex(sol))
                                solCounter+=1
                        else:
                            F_sol = solveset(F, v, domain=Interval(Ymin, Ymax))
                            for sol in F_sol:
                                self.ax.hlines(sol, Xmin, Xmax, color=self.plotColors[solCounter])
                                if lableGen: self.latex_code.append("v = " + latex(sol))
                                solCounter+=1
                    else: print("Not interesting")
                else:
                    if len(F_symbols) == 2:
                        F_sol = solve(F)
                        for sol in F_sol:
                            if u in sol:
                                if sol[u].has(v):
                                    F_sol_num = lambdify(v, sol[u], "numpy")
                                    self.ax.plot(F_sol_num(y), y, c=self.plotColors[solCounter])
                                    if lableGen: self.latex_code.append("u = " + latex(sol[u]))
                                else:
                                    self.ax.vlines(sol[u], Ymin, Ymax, color=self.plotColors[solCounter])
                                    if lableGen: self.latex_code.append("u = " + latex(sol[u]))
                            else:
                                if sol[v].has(u):
                                    F_sol_num = lambdify(u, sol[v], "numpy")
                                    self.ax.plot(x, F_sol_num(x), c=self.plotColors[solCounter])
                                    if lableGen: self.latex_code.append("v = " + latex(sol[v]))
                                else:
                                    self.ax.hlines(sol[v], Xmin, Xmax, color=self.plotColors[solCounter])
                                    if lableGen: self.latex_code.append("v = " + latex(sol[v]))
                            solCounter+=1

                    elif len(F_symbols) == 1:
                        if u in F_symbols:
                            F_sol = solve(F, u, domain=Interval(Xmin, Xmax))
                            for sol in F_sol:
                                self.ax.vlines(sol, Ymin, Ymax, color=self.plotColors[solCounter])
                                if lableGen: self.latex_code.append("u = " + latex(sol))
                                solCounter+=1
                        else:
                            F_sol = solve(F, v, domain=Interval(Ymin, Ymax))
                            for sol in F_sol:
                                self.ax.hlines(sol, Xmin, Xmax, color=self.plotColors[solCounter])
                                if lableGen: self.latex_code.append("v = " + latex(sol))
                                solCounter+=1
                    else: print("Not interesting")
        
        self.render_latex()
        self.latex_code = []
        self.canvas.draw()


app = QApplication([])

window = MainWindow()
window.show()

# Start the event loop.
app.exec()