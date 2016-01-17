#!/usr/bin/env python
import sys
import os
import pylab
import dendropy

from dendropy.calculate import treecompare
from cStringIO import StringIO
from Bio.Nexus import Trees
##from PyQt4 import QtCore
from PyQt4 import QtGui
from Bio import Phylo


# Klasa odpowiadajaca za konwersje
class ConvertWindow(QtGui.QWidget):
    chosenFileName = ''
    chosenInputFormat = ''
    chosenOutputFormat = ''

    def __init__(self):
        super(ConvertWindow, self).__init__()
        self.initUI()

    def initUI(self):
        # labels
        self.inputLabel = QtGui.QLabel('Wejscie')
        self.outputLabel = QtGui.QLabel('Wyjscie')

        # fields for files content
        self.textEditIn = QtGui.QTextEdit()
        self.textEditIn.setReadOnly(True)
        self.textEditOut = QtGui.QTextEdit()
        self.textEditOut.setReadOnly(True)

        # format groups
        self.inputFormatGroup = self.createInputFormatExclusiveGroup()
        self.outpuFormatGroup = self.createOutputFormatExclusiveGroup()

        # button for opening file to convert
        self.openBtn = QtGui.QPushButton('Wybierz plik')
        self.openBtn.clicked.connect(lambda: self.showOpenFileDialog(self.textEditIn))
        self.convBtn = QtGui.QPushButton('Konwertuj')
        self.convBtn.clicked.connect(lambda: self.convertTreeFile(self.textEditIn, self.textEditOut))

        # window layout
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)

        self.grid.addWidget(self.openBtn, 1, 0)
        self.grid.addWidget(self.convBtn, 1, 1)

        self.grid.addWidget(self.inputLabel, 2, 0)
        self.grid.addWidget(self.inputFormatGroup, 2, 1)
        self.grid.addWidget(self.outputLabel, 2, 6)
        self.grid.addWidget(self.outpuFormatGroup, 2, 7)

        self.grid.addWidget(self.textEditIn, 3, 0, 5, 6)
        self.grid.addWidget(self.textEditOut, 3, 6, 5, 10)

        self.setLayout(self.grid)
        self.setGeometry(200, 200, 600, 500)
        self.setWindowTitle('Konwerter')
        self.show()

    def createInputFormatExclusiveGroup(self):
        self.groupBox = QtGui.QGroupBox("Format pliku wejsciowego")

        self.radio1 = QtGui.QRadioButton("Newick")
        self.radio2 = QtGui.QRadioButton("Nexus")
        self.radio3 = QtGui.QRadioButton("PhyloXML")

        self.radio1.clicked.connect(lambda: self.rememberInputFormat(self.radio1.text().toLower()))
        self.radio2.clicked.connect(lambda: self.rememberInputFormat(self.radio2.text().toLower()))
        self.radio3.clicked.connect(lambda: self.rememberInputFormat(self.radio3.text().toLower()))

        self.radio1.setChecked(True)
        self.rememberInputFormat(self.radio1.text().toLower())

        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addWidget(self.radio1)
        self.vbox.addWidget(self.radio2)
        self.vbox.addWidget(self.radio3)
        self.groupBox.setLayout(self.vbox)

        return self.groupBox

    def createOutputFormatExclusiveGroup(self):
        self.groupBox = QtGui.QGroupBox("Format pliku wyjsciowego")

        self.radio1 = QtGui.QRadioButton("Newick")
        self.radio2 = QtGui.QRadioButton("Nexus")
        self.radio3 = QtGui.QRadioButton("PhyloXML")

        self.radio1.clicked.connect(lambda: self.rememberOutputFormat(self.radio1.text().toLower()))
        self.radio2.clicked.connect(lambda: self.rememberOutputFormat(self.radio2.text().toLower()))
        self.radio3.clicked.connect(lambda: self.rememberOutputFormat(self.radio3.text().toLower()))

        self.radio2.setChecked(True)
        self.rememberOutputFormat(self.radio2.text().toLower())

        self.vbox = QtGui.QVBoxLayout()
        self.vbox.addWidget(self.radio1)
        self.vbox.addWidget(self.radio2)
        self.vbox.addWidget(self.radio3)
        self.groupBox.setLayout(self.vbox)

        return self.groupBox

    def rememberInputFormat(self, formatName):
        self.chosenInputFormat = formatName

    def rememberOutputFormat(self, formatName):
        self.chosenOutputFormat = formatName

    def showOpenFileDialog(self, textEditIn=None):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik', './Trees')
        if fname != '':
            self.chosenFileName = str(fname)
            f = open(fname, 'r')

            with f:
                data = f.read()
                if not (textEditIn is None):
                    textEditIn.setText(data)

    def convertTreeFile(self, inputTextEdit, outputTextEdit):
        if self.chosenFileName is '':
            self.showOpenFileDialog(inputTextEdit)

        # zamiana
        if self.chosenInputFormat != self.chosenOutputFormat:
            if self.chosenFileName != '' and self.chosenInputFormat != '':
                self.convertedFileName = str(self.chosenFileName).replace(
                        '.' + str(self.chosenInputFormat), '.' + str(self.chosenOutputFormat))

                Phylo.convert(str(self.chosenFileName), str(self.chosenInputFormat),
                              self.convertedFileName, str(self.chosenOutputFormat))

                f = open(self.convertedFileName, 'r')

                with f:
                    data = f.read()
                    outputTextEdit.setText(data)
