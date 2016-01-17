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
##from PyQt4 import QtCore
from PyQt4 import QtGui
from Bio import Phylo

# zewnetrzne pliki
import ConvertClass
import ConsensusWindow as cs
import MainClass as ms

# Klasa wyswietlajaca informacje o drzewie
# Napisana w innym terminie -> w j. angielskim
"""
██ ███    ██ ███████  ██████  ██     ██ ██ ███    ██ ██████   ██████  ██     ██
██ ████   ██ ██      ██    ██ ██     ██ ██ ████   ██ ██   ██ ██    ██ ██     ██
██ ██ ██  ██ █████   ██    ██ ██  █  ██ ██ ██ ██  ██ ██   ██ ██    ██ ██  █  ██
██ ██  ██ ██ ██      ██    ██ ██ ███ ██ ██ ██  ██ ██ ██   ██ ██    ██ ██ ███ ██
██ ██   ████ ██       ██████   ███ ███  ██ ██   ████ ██████   ██████   ███ ███
"""
class InfoWindow(QtGui.QMainWindow):
    def __init__(self, tree):
        super(InfoWindow, self).__init__()
        self.tree = tree
        self.textt = 'Information about tree' + "\n" + str(tree)
        #self.tree.display()
        self.initUI()

    def showAscii(self):
        # print tree
        self.tmpf = open('/tmp/ascii.txt', 'w')
        Phylo.draw_ascii(self.tree, self.tmpf)

        self.tmpf = open('/tmp/ascii.txt', 'r')
        with self.tmpf:
            self.textt += "\n" + self.tmpf.read()

    def initUI(self):
        self.getMinMaxDepths()
        self.getMinMaxLengths()
        self.numberOfTermianls = self.tree.count_terminals()
        self.numberOfNonterminals = len(self.tree.get_nonterminals())



        self.isBifurcating = ''
        if self.tree.is_bifurcating():
            self.isBifurcating = 'Yes'
        else:
            self.isBifurcating = 'No :<'

        self.textEdit = QtGui.QTextEdit()
        self.textEdit.setReadOnly(True)
        self.textEdit.setFontFamily('Courier')
        self.textEdit.setWordWrapMode(True)

        # self.makeStringInfo()
        self.textEdit.setText(self.textt)
        self.setCentralWidget(self.textEdit)

        # labels
        self.showAscii()
        self.textt += '\nNumber of terminals and non-terminals: ' + str(
            self.numberOfTermianls + self.numberOfNonterminals)
        self.textt += '\nNumber of terminals: ' + str(self.numberOfTermianls)
        self.textt += '\nNumber of non-terminals: ' + str(self.numberOfNonterminals)
        self.textt += '\nMin depth: ' + str(self.minDepth)
        self.textt += '\nMax depth: ' + str(self.maxDepth)
        self.textt += '\nMin length: ' + str(self.minLength)
        self.textt += '\nMax length: ' + str(self.maxLength)
        self.textt += '\nTotal branch length: ' + str(self.tree.total_branch_length())
        self.textt += '\nIs tree is bifurcanting: ' + self.isBifurcating

        self.textEdit.setText(self.textt)
        self.setGeometry(200, 200, 800, 800)
        self.setWindowTitle('Information')
        self.show()

    def getMinMaxDepths(self):
        self.minDepth = -1
        self.maxDepth = -1
        for clade in self.tree.depths(unit_branch_lengths=True):
            if not (clade.name is None):
                self.currentDepth = self.tree.depths(unit_branch_lengths=True)[clade]
                if self.currentDepth < self.minDepth or self.minDepth == -1:
                    self.minDepth = int(self.currentDepth)

                if self.currentDepth > self.maxDepth or self.maxDepth == -1:
                    self.maxDepth = int(self.currentDepth)


    def getMinMaxLengths(self):
        self.minLength = -1.0
        self.maxLength = -1.0
        for clade in self.tree.depths():
            if not (clade.name is None):
                self.currentLength = self.tree.depths()[clade]
                if self.currentLength < self.minLength or self.minLength == -1.0:
                    self.minLength = self.currentLength

                if self.currentLength > self.maxLength or self.maxLength == -1.0:
                    self.maxLength = self.currentLength
