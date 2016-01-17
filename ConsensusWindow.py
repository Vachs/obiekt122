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

import MainClass as mc
import InfoWindow as iw

class ConsensusWindow(QtGui.QWidget):
    def __init__(self):
        super(ConsensusWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.initUI()

    def initUI(self):
        # sciezki
        self.label1 = QtGui.QLabel(self.path1)
        self.label2 = QtGui.QLabel(self.path2)

        # buttony
        self.path1Btn = QtGui.QPushButton('Wybierz pierwsze drzewo')
        self.path1Btn.clicked.connect(self.path2File1)

        self.path2Btn = QtGui.QPushButton('Wybierz drugie drzewo')
        self.path2Btn.clicked.connect(self.path2File2)

        self.consTreeBtn = QtGui.QPushButton('Wygeneruj drzewo konsensusu')
        self.consTreeBtn.clicked.connect(self.drawConsensusTreeBio)

        self.consInfoBtn = QtGui.QPushButton('Pokaz informacje o tym drzewie')
        self.consInfoBtn.clicked.connect(self.showConsensusTreeBio)

        self.otherRF = QtGui.QPushButton('Oblicz odleglosci')
        self.otherRF.clicked.connect(self.showOtherRF)

        self.lay1 = QtGui.QLabel('Metryka Robinsona-Fouldsa: ')
        self.lay2 = QtGui.QLabel('Odleglosc euklidesowa: ')
        #self.lay2 = QtGui.QLabel('Roznica symetryczna')
        #self.lay3 = QtGui.QLabel('Falszywe pozytywy i negatywy')

        self.res1 = QtGui.QLabel('')
        self.res2 = QtGui.QLabel('')
        #self.res3 = QtGui.QPushButton('')
        #self.res4 = QtGui.QPushButton('')

        # obraz
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(12)

        self.grid.addWidget(self.path1Btn, 1, 0)
        self.grid.addWidget(self.path2Btn, 2, 0)
        self.grid.addWidget(self.label1, 3, 0)
        self.grid.addWidget(self.label2, 4, 0)

        self.grid.addWidget(self.consTreeBtn, 5, 0)
        self.grid.addWidget(self.consInfoBtn, 6, 0)
        self.grid.addWidget(self.otherRF, 7, 0)

        self.grid.addWidget(self.lay1, 1, 1)
        self.grid.addWidget(self.lay2, 2, 1)
        self.grid.addWidget(self.res1, 1, 2)
        self.grid.addWidget(self.res2, 2, 2)

        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 50)
        self.setWindowTitle('Drzewa konsensusu')
        self.show()

    def showOpenFileDialog(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik z drzewem')
        return str(fname)

    def path2File1(self):
        self.path1 = self.showOpenFileDialog()
        self.label1.setText(self.path1)

    def path2File2(self):
        self.path2 = self.showOpenFileDialog()
        self.label2.setText(self.path2)

    def OpenInfoWindow(self):
        if self.tree != 0:
            self.infoWin = iw.InfoWindow(self.tree)
            self.infoWin.show()

    def showOtherRF(self):
        if self.path1 != '' and self.path2 != '':
            self.calcDistance()
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

    def showConsensusTreeBio(self):
        if self.path1 != '' and self.path2 != '':
            # odczytaj rozszerzenie
            self.fileEx1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileEx2 = (os.path.splitext(self.path2)[1])[1:]

            self.trees = []

            # pierwsze
            self.f = open(self.path1, 'r')
            self.tree1 = Trees.Tree(self.f.read())
            self.trees.append(self.tree1)
            self.f.close()

            # drugie
            self.f = open(self.path2, 'r')
            self.tree2 = Trees.Tree(self.f.read())
            self.trees.append(self.tree2)
            self.f.close()

            self.consensus_tree = Trees.consensus(self.trees)

            # SHOWTIME
            # rysuj
            self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            #self.tree.root.color = '#808080'
            self.OpenInfoWindow()
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

    def drawConsensusTreeBio(self):
        if self.path1 != '' and self.path2 != '':
            self.fileEx1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileEx2 = (os.path.splitext(self.path2)[1])[1:]

            self.trees = []

            self.f = open(self.path1, 'r')
            self.tree1 = Trees.Tree(self.f.read())
            self.trees.append(self.tree1)
            self.f.close()

            self.f = open(self.path2, 'r')
            self.tree2 = Trees.Tree(self.f.read())
            self.trees.append(self.tree2)
            self.f.close()

            self.consensus_tree = Trees.consensus(self.trees)

            # SHOWTIME
            self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            self.tree.root.color = '#808080'
            #self.OpenInfoWindow()
            Phylo.draw(self.tree)
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

    def calcDistance(self):
        if self.path1 != '' and self.path2 != '':
            self.fileEx1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileEx2 = (os.path.splitext(self.path2)[1])[1:]

            tns = dendropy.TaxonNamespace()
            self.tree1 = dendropy.Tree.get_from_path(self.path1, self.fileEx1, taxon_namespace=tns)
            self.tree2 = dendropy.Tree.get_from_path(self.path2, self.fileEx2, taxon_namespace=tns)

            self.tree1.encode_bipartitions()
            self.tree2.encode_bipartitions()

            print(treecompare.false_positives_and_negatives(self.tree1, self.tree2))

            # self.tree1 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')
            # self.tree2 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')

            # self.tree1.encode_bipartitions()
            # self.tree2.encode_bipartitions()

            # oblicz dystans
            # self.symDist = self.tree1.symmetric_difference(self.tree2)
            self.symDist = treecompare.symmetric_difference(self.tree1, self.tree2)
            self.fpnDist = treecompare.false_positives_and_negatives(self.tree1, self.tree2)
            self.eucDist = treecompare.euclidean_distance(self.tree1, self.tree2)
            self.rfDist = treecompare.robinson_foulds_distance(self.tree1, self.tree2)

            # pokaz wyniki
            self.res1.setText(str(self.eucDist)) #eucDist
            self.res2.setText(str(self.rfDist))  #rfDist
            #self.lay3.setText(str(self.symDist)) #symDist
            #self.lay4.setText(str(self.fpnDist)) #fpnDist

    # ROOT-uje drzewo
    # Usuwa/Dodaje korzen
    def makeRootUnroot(self, mod):
        if self.path1 != '' and self.path2 == '':
            # get files extensions
            self.fileEx1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileEx2 = (os.path.splitext(self.path2)[1])[1:]

            # open tree files
            self.trees = []
            self.drzewo = []

            # first tree
            self.f = open(self.path1, 'r')
            self.miss = self.f.read()
            self.tree1 = Trees.Tree(self.miss)
            self.dre = ts.Tree(self.miss)
            print "# Before modification"
            print self.tree1

            if mod == 0:
                print "# After modification -- Rooting (at midpoint):"
                self.dre.root_midpoint()
            elif mod == 1:
                print "# After modification -- UnRooting:"
                self.dre.unroot()
            elif mod == 2:
                print "# After modification -- Rooting (balanced):"
                self.dre.root_balanced()
            print self.dre
            print "\nDetails about tree:"
            self.dre.display()
            Phylo.draw_ascii(self.tree1)
            self.show()
            self.f.close()
