#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import pylab
import Tree as ts
import dendropy

from dendropy.calculate import treecompare
from cStringIO import StringIO
from Bio.Nexus import Trees
from PyQt4 import QtGui
from Bio import Phylo

# zewnetrzne pliki
import ConsensusWindow as cs
import InfoWindow as iw

"""
>TODO
>> all looks completed, nothing to add. Tomasz
"""


# Okno wyboru kolejnych opcji. Wybierz co chcesz zrobic dalej
"""
███    ██ ███████ ██   ██ ████████ ██     ██ ██ ███    ██ ██████   ██████  ██     ██
████   ██ ██       ██ ██     ██    ██     ██ ██ ████   ██ ██   ██ ██    ██ ██     ██
██ ██  ██ █████     ███      ██    ██  █  ██ ██ ██ ██  ██ ██   ██ ██    ██ ██  █  ██
██  ██ ██ ██       ██ ██     ██    ██ ███ ██ ██ ██  ██ ██ ██   ██ ██    ██ ██ ███ ██
██   ████ ███████ ██   ██    ██     ███ ███  ██ ██   ████ ██████   ██████   ███ ███
"""
class NextWindow(QtGui.QWidget):
    chosenFileName = ''
    tree = 0
    path1 = ''

    def __init__(self):
        super(NextWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.cpath = ''
        self.cpath2 = ''
        self.initUI()

    def initUI(self):
        self.showOpenFileWindow()
        #Trees.Tree.display(self.tree)

        # WYSZUIKWANIE DROGI PO KRAWEDZIACH
        # labels
        self.label1 = QtGui.QLabel('Pierwszy wezel')
        self.label2 = QtGui.QLabel('Drugi wezel')

        # combo boxes
        self.combobox1 = QtGui.QComboBox(self)
        self.combobox2 = QtGui.QComboBox(self)

        # DODAJE PUSTY WEZEL
        self.combobox1.addItem('')
        self.combobox2.addItem('')

        for clade in self.terminals:
            self.combobox1.addItem(clade.name)
            self.combobox2.addItem(clade.name)

        self.combobox1.activated[str].connect(self.combo1)
        self.combobox2.activated[str].connect(self.combo2)

        self.combobox1.activated[str].connect(self.combo1)
        self.combobox2.activated[str].connect(self.combo2)

        self.pathBtn = QtGui.QPushButton('wyznacz')
        self.pathBtn.clicked.connect(self.showPathWindow)

        # buttony
        self.info = QtGui.QLabel('Kalkulator drzew filogenicznych. Mozliwosc wczytania zarowno drzew ukorzenionych\n'
                                 'jak i nieukorzenionych. Wybierz jedna z opcji dostepnych ponizej.')

        self.lay0 = QtGui.QLabel('>> OPERACJA NA JEDNYM DRZEWIE')
        self.lay1 = QtGui.QLabel('Pozwala uzyskac podstawowe informacje o wybranych drzewie')
        self.btn1 = QtGui.QPushButton('Uzyskaj informacje o drzewie')
        self.btn1.clicked.connect(self.OpenInfoWindow)

        self.lay2 = QtGui.QLabel('Pozwala na prosta wizualizacje drzewa')
        self.btn2 = QtGui.QPushButton('Rysuj drzewo')
        self.btn2.clicked.connect(self.DrawSimple)

        self.btn22 = QtGui.QPushButton('Bardziej zaawansowane z wartosciami')
        self.btn22.clicked.connect(self.DrawAdvan)

        self.lay3 = QtGui.QLabel('Pozwala na wyznaczenie najkrotszej trasy pomiedzy dwoma liscmi')
        self.btn3 = QtGui.QPushButton('Wyznacz trase')

        self.lay4 = QtGui.QLabel('>> OPERACJE NA DWOCH DRZEWACH')
        self.lay5 = QtGui.QLabel('Pozwala na wyrysowanie drzewa konsensusu')
        self.btn4 = QtGui.QPushButton('Wybierz drzewa')
        self.btn4.clicked.connect(self.showConsensusWindow)


        self.anotherTree = QtGui.QPushButton('Wczytaj inne drzewo')
        self.anotherTree.clicked.connect(self.showOpenFileWindow)

        # obraz
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)

        self.grid.addWidget(self.info)
        self.grid.addWidget(self.anotherTree, 2, 0)
        self.grid.addWidget(self.lay0, 3, 0)
        self.grid.addWidget(self.lay1, 4, 0)
        self.grid.addWidget(self.btn1, 5, 0)

        self.grid.addWidget(self.lay2, 6, 0)
        self.grid.addWidget(self.btn2, 7, 0)
        self.grid.addWidget(self.btn22, 8, 0)

        self.grid.addWidget(self.lay3, 9, 0)
        self.grid.addWidget(self.label1, 9, 1)
        self.grid.addWidget(self.label2, 9, 2)
        self.grid.addWidget(self.pathBtn, 10, 0)
        self.grid.addWidget(self.combobox1, 10, 1)
        self.grid.addWidget(self.combobox2, 10, 2)

        self.grid.addWidget(self.lay4, 11, 0)
        self.grid.addWidget(self.lay5, 12, 0)
        self.grid.addWidget(self.btn4, 13, 0)

        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 50)
        self.setWindowTitle('Menu')
        self.show()

    def showConsensusWindow(self):
        self.consWin = cs.ConsensusWindow()
        self.consWin.show()

    def showOpenFileWindow(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik zawierajacy drzewo', './Trees')
        if fname != '':
            self.chosenFileName = str(fname)
            f = open(str(fname), 'r')
            self.tree1 = Trees.Tree(f.read())
            self.moredata = Trees.Tree.display(self.tree1)
            print self.moredata
            self.tmpf = open('more.txt', 'w')
            self.tmpf.write(str(self.moredata))
            self.tmpf.close()
            # with f:
            # path1 = f.read()
            # self.pathInfo.setText(path1)

            self.fileExtension = (os.path.splitext(str(fname))[1])[1:]
            self.tree = Phylo.read(str(fname), self.fileExtension)
            self.tree1s = self.tree
            self.terminals = self.tree1s.get_terminals()

    def OpenInfoWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileWindow()

        if self.tree != 0:
            self.infoWin = iw.InfoWindow(self.tree)
            self.infoWin.show()

    # RYSUJE PODSTAWOWY GRAF
    def DrawSimple(self):
        if self.chosenFileName == '':
            self.showOpenFileWindow()

        if self.tree != 0:
            self.tree.root.color = '#808080'
            Phylo.draw_graphviz(self.tree, node_size=2500)
            pylab.show()

    # RYSUJE BARDZIEJ ZAAWANSOWANY WYKRES
    def DrawAdvan(self):
        if self.chosenFileName == '':
            self.showOpenFileWindow()

        if self.tree != 0:
            self.tree.root.color = '#808080'
            Phylo.draw(self.tree, branch_labels=lambda c: c.branch_length)

    def combo1(self, text):
        self.cpath = str(text)

    def combo2(self, text):
        self.cpath2 = str(text)

    def showPathWindow(self):
        if self.cpath != '' and self.cpath2 != '':
            self.start = self.tree.find_clades(self.cpath).next()
            self.end = self.tree.find_clades(name=self.cpath2).next()

            for clade in self.tree.trace(self.start, self.end):
                clade.color = 'red'

            for clade in self.tree.find_clades():
                if not (clade in self.tree.trace(self.start, self.end)):
                    clade.color = 'grey'

            self.start.color = 'blue'

            # RYSOWANIE
            Phylo.draw_graphviz(self.tree, node_size=2500)
            pylab.plot(0, 0)

            # DAJE CZYSTY OBRAZ BEZ OSI ** PEWNIE MOZNA PROSCIEJ
            frame1 = pylab.gca()
            for xlabel_i in frame1.axes.get_xticklabels():
                xlabel_i.set_visible(False)
                xlabel_i.set_fontsize(0.0)
            for xlabel_i in frame1.axes.get_yticklabels():
                xlabel_i.set_fontsize(0.0)
                xlabel_i.set_visible(False)
            for tick in frame1.axes.get_xticklines():
                tick.set_visible(False)
            for tick in frame1.axes.get_yticklines():
                tick.set_visible(False)
            pylab.show()

        else:
            print "Nie wybrano punktow"

            # WYSWIETLA INFORMACJE
            img = pylab.imread('img/wally.png', 'rb')
            pylab.imshow(img)
            pylab.plot(0, 0)

            # DAJE CZYSTY OBRAZ BEZ OSI ** PEWNIE MOZNA PROSCIEJ
            frame1 = pylab.gca()
            for xlabel_i in frame1.axes.get_xticklabels():
                xlabel_i.set_visible(False)
                xlabel_i.set_fontsize(0.0)
            for xlabel_i in frame1.axes.get_yticklabels():
                xlabel_i.set_fontsize(0.0)
                xlabel_i.set_visible(False)
            for tick in frame1.axes.get_xticklines():
                tick.set_visible(False)
            for tick in frame1.axes.get_yticklines():
                tick.set_visible(False)

            # SHOWTIME
            pylab.show()

"""
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
"""
def main():
    app = QtGui.QApplication(sys.argv)
    # mw = MainWindow()
    wd = NextWindow()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

