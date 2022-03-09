#!/usr/bin/env python

import os
import sys
import pathlib
import appdirs
import configparser

from PyQt5.QtWidgets import (
    QApplication, QFormLayout, QWidget, QLineEdit,
    QMainWindow, QMenu, QAction, QVBoxLayout, QDialog,
    QDialogButtonBox, QTextEdit, QMessageBox, QFileDialog, QPushButton,
)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QIcon

from simulator import Simulator


class ClickableLineEdit(QLineEdit):
    clicked = pyqtSignal()

    def mousePressEvent(self, e):
        if e.button() == Qt.LeftButton:
            self.clicked.emit()
        else:
            super().mousePressEvent(e)


class Window(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mosquito lifecycle")

        # icon
        if getattr(sys, 'frozen', False):
            fname = pathlib.Path(sys._MEIPASS) / 'andytest.ico'
        else:
            me = pathlib.Path(sys.argv[0]).resolve()
            fname = me.parent.parent / 'etc' / 'andytest.ico'
        # on linux we need a ppm/pgm file, not a windows ico, just ignore
        try:
            self.setWindowIcon(QIcon(str(fname)))
        except Exception:
            pass

        # vlayout contains the form layout with Go button at the bottom
        vlayout = QVBoxLayout()
        flayout = QFormLayout()

        self.in_params = in_params = ClickableLineEdit()
        in_params.setToolTip('CSV file quantifying the parameters')
        flayout.addRow("Parameters file:", in_params)
        in_params.clicked.connect(self.set_in_params)

        self.out_fn = out_fn = ClickableLineEdit()
        out_fn.setToolTip('Output will be written to this file')
        flayout.addRow("Output file:", out_fn)
        out_fn.clicked.connect(self.set_out_fn)

        vlayout.addLayout(flayout)

        self.but = QPushButton('GO', self)
        self.but.setToolTip('Run mosquito lifecycle')
        self.but.setStyleSheet('color: green')
        self.but.clicked.connect(self.go)
        vlayout.addWidget(self.but)

        widget = QWidget()
        widget.setLayout(vlayout)

        # config file
        adirs = appdirs.AppDirs('Mosquito lifecycle', 'CSIRO')
        self.cdir = cdir = pathlib.Path(adirs.user_data_dir)
        os.makedirs(cdir, mode=0o755, exist_ok=True)
        self.cfile = pathlib.Path(adirs.user_data_dir) / 'config.ini'
        print(self.cfile)
        self.__load_config()

        self.__create_menus()

        self.setCentralWidget(widget)

    def __create_menus(self):
        bar = self.menuBar()
        self.setMenuBar(bar)

        menu = QMenu("&File", self)
        bar.addMenu(menu)
        action = QAction("&Quit", self)
        action.triggered.connect(self.close)
        menu.addAction(action)

        menu = QMenu("&Help", self)
        bar.addMenu(menu)
        action = QAction("&About", self)
        action.triggered.connect(lambda: self.msg('About', 'version.txt'))
        menu.addAction(action)
        action = QAction("&Help", self)
        action.triggered.connect(lambda: self.msg('Help', 'help.html'))
        menu.addAction(action)

    def __load_config(self):
        if not self.cfile.exists():
            cp = configparser.ConfigParser()
            cp['DEFAULT'] = {
                'in_dir': pathlib.Path.home() / 'Downloads',
                'out_dir': pathlib.Path.home() / 'Downloads',
            }
            with open(self.cfile, 'w') as cfh:
                cp.write(cfh)

        self.cp = configparser.ConfigParser()
        self.cp.read(self.cfile)

    def save_config(self):
        with open(self.cfile, 'w') as cfh:
            self.cp.write(cfh)

    def msg(self, title, fname):

        class ScrollableDialog(QDialog):

            def __init__(self, title, txt, parent=None, html=False):
                super().__init__(parent)
                self.setWindowTitle(title)
                btns = QDialogButtonBox(QDialogButtonBox.Ok)
                btns.accepted.connect(self.accept)
                layout = QVBoxLayout()
                te = QTextEdit()
                if html:
                    te.setHtml(txt)
                else:
                    te.setText(txt)
                te.setMinimumWidth(600)
                te.setMinimumHeight(400)
                te.setReadOnly(True)
                layout.addWidget(te)
                layout.addWidget(btns)
                self.setLayout(layout)

        if getattr(sys, 'frozen', False):
            fname = pathlib.Path(sys._MEIPASS) / fname
        else:
            me = pathlib.Path(sys.argv[0]).resolve()
            fname = me.parent / fname

        with open(fname) as fh:
            txt = fh.read()
            # if txt is short a messagebox will do
            if txt.count('\n') < 5:
                QMessageBox.information(self, 'About', txt)
            else:
                dl = ScrollableDialog(
                    title, txt, parent=self,
                    html=pathlib.Path(fname).suffix == '.html'
                )
                dl.exec()

    def set_in_params(self):
        fname, check = QFileDialog.getOpenFileName(
            self, "CSV input parameters", self.cp['DEFAULT']['in_dir'],
            "CSV File (*.csv)",
        )
        if not check:
            return

        # make sure extension is .csv
        if pathlib.Path(fname).suffix != '.csv':
            fname += '.csv'

        # so we use the same directory next time
        self.cp['DEFAULT']['in_dir'] = str(pathlib.Path(fname).parent)
        self.save_config()

        self.in_params.setText(fname)

    def set_out_fn(self):
        fname, check = QFileDialog.getSaveFileName(
            self, "Output file", self.cp['DEFAULT']['out_dir'],
            "CSV File (*.csv)",
        )
        if not check:
            return

        # make sure extension is .csv
        if pathlib.Path(fname).suffix != '.csv':
            fname += '.csv'

        # so we use the same directory next time
        self.cp['DEFAULT']['out_dir'] = str(pathlib.Path(fname).parent)
        self.save_config()

        self.out_fn.setText(fname)

    def go(self):
        in_params = self.in_params.text()
        out_fn = self.out_fn.text()

        #fibo = routines.fib(10)
        #QMessageBox.information(self, 'Fibo', ' '.join(map(str, fibo)))
        #return 

        if not in_params:
            QMessageBox.critical(
                self, "Error", "Provide a CSV file quantifying the parameters"
            )
            return False
        if not out_fn:
            QMessageBox.critical(
                self, "Error", "Must provide an output filename"
            )
            return False

        simulator = Simulator(in_params)
        simulator.go()
        simulator.output(out_fn)
        




if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec_())
