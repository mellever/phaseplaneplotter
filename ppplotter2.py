#Import packages
from PyQt6.QtCore import QSize, Qt
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import *

import numpy as np

from sympy import sympify, lambdify
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


        #Create button for plotting
        PlotButton = QPushButton("Plot")
        PlotButton.setCheckable(True)
        PlotButton.clicked.connect(self.plot_button)
        
        #Create input for F1
        self.InputF1 = QLineEdit()
        self.InputF1.setText("cos(u)")

        #Create input for F2
        self.InputF2 = QLineEdit()
        self.InputF2.setText("sin(v)")

        #Create input for Xmin, Xmax, Ymin, Ymax
        self.InputXmin = QLineEdit()
        self.InputXmin.setText("0")
        self.InputXmax = QLineEdit()
        self.InputXmax.setText("5")
        self.InputYmin = QLineEdit()
        self.InputYmin.setText("0")
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

        layout.addWidget(LabelDensity, 6, 0)
        layout.addWidget(self.Slider, 6, 1)
        
        layout.addWidget(LabelXmin, 2, 0)
        layout.addWidget(self.InputXmin, 2, 1)

        layout.addWidget(LabelXmax, 3, 0)
        layout.addWidget(self.InputXmax, 3, 1)

        layout.addWidget(LabelYmin, 4, 0)
        layout.addWidget(self.InputYmin, 4, 1)

        layout.addWidget(LabelYmax, 5, 0)
        layout.addWidget(self.InputYmax, 5, 1)
             
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

        F1 = lambdify([u,v], sympify(F1))
        F2 = lambdify([u,v], sympify(F2))

        U = F1(X, Y)
        V = F2(X, Y)
        
        self.ax.clear()
        self.ax.set_xlabel("u")
        self.ax.set_ylabel("v")
        plt.tight_layout()
        self.ax.streamplot(X, Y, U, V, density=self.Slider.value())
        self.canvas.draw()


app = QApplication([])

window = MainWindow()
window.show()

# Start the event loop.
app.exec()