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

class ConsensusWindow(QtGui.QWidget):    
    def __init__(self):
        super(ConsensusWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.initUI()
        
    def initUI(self):
        #labels
        self.path1Label = QtGui.QLabel(self.path1)
        self.path2Label = QtGui.QLabel(self.path2)        
        
        #buttons
        self.path1Btn = QtGui.QPushButton('Wybierz pierwsze drzewo')
        self.path1Btn.clicked.connect(self.chooseFile1)
        
        self.path2Btn = QtGui.QPushButton('Wybierz drugie drzewo')
        self.path2Btn.clicked.connect(self.chooseFile2)
        
        self.consDenBtn = QtGui.QPushButton('Wygeneruj drzewo konsensusu (dendropy)')
        self.consDenBtn.clicked.connect(self.drawConsensusTreeDendropy)
        
        self.consNexBtn = QtGui.QPushButton('Wygeneruj drzewo konsensusu (Bio.Nexus)')
        self.consNexBtn.clicked.connect(self.drawConsensusTreeBioNexus)
        
        #window layout
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)
        
        self.grid.addWidget(self.path1Label, 1, 0)
        self.grid.addWidget(self.path1Btn, 1, 1)
        
        self.grid.addWidget(self.path2Label, 2, 0)
        self.grid.addWidget(self.path2Btn, 2, 1)

        self.grid.addWidget(self.consDenBtn, 3, 0)
        self.grid.addWidget(self.consNexBtn, 3, 1)
        
        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 50)
        self.setWindowTitle('Drzewa konsensusu')
        self.show()
        
    def showOpenFileDialog(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik', './Trees')
        return str(fname)
    
    
    def chooseFile1(self):
        self.path1 = self.showOpenFileDialog()
        self.path1Label.setText(self.path1)
    
    
    def chooseFile2(self):
        self.path2 = self.showOpenFileDialog()
        self.path2Label.setText(self.path2)
        
    
    def drawConsensusTreeDendropy(self):
        if self.path1 != '' and self.path2 != '':
            #get files extensions
            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]
            
            #open tree files
            self.tree1 = dendropy.Tree.get_from_path(self.path1, self.fileExtension1)
            self.tree2 = dendropy.Tree.get_from_path(self.path2, self.fileExtension2)
            
            #prepare tree list
            self.trees = dendropy.TreeList()
            self.trees.append(self.tree1)
            self.trees.append(self.tree2)
            
            #generate consensus tree
            self.consensus_tree = self.trees.consensus(min_freq=0.2)
            
            #draw tree
            self.handle = StringIO(self.consensus_tree._as_newick_string())
            # POPRAWIONY BLAD Z KONWERSJA DO BUFORA

            #self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            self.tree.root.color = '#808080'
            Phylo.draw(self.tree)
            
            
    def drawConsensusTreeBioNexus(self):
        if self.path1 != '' and self.path2 != '':
            #get files extensions
            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]
            
            #open tree files            
            self.trees = []
            
            #first tree
            self.f = open(self.path1, 'r')
            self.tree1 = Trees.Tree(self.f.read())
            self.trees.append(self.tree1)
            self.f.close()
            
            #second tree
            self.f = open(self.path2, 'r')
            self.tree2 = Trees.Tree(self.f.read())
            self.trees.append(self.tree2)
            self.f.close()


            #generate consensus tree
            self.consensus_tree = Trees.consensus(self.trees)
            
            #draw tree
            self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            self.tree.root.color = '#808080'
            Phylo.draw(self.tree)
            
            

class DistancesWindow(QtGui.QWidget):    
    def __init__(self):
        super(DistancesWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.initUI()
        
    def initUI(self):
        #labels
        self.path1Label = QtGui.QLabel(self.path1)
        self.path2Label = QtGui.QLabel(self.path2)
        
        self.dist1Label = QtGui.QLabel('Odleglosc euklidesowa')
        self.dist1Value = QtGui.QLabel('')
        
        self.dist2Label = QtGui.QLabel('Metryka Robinsona-Fouldsa')
        self.dist2Value = QtGui.QLabel('')
        
        self.dist3Label = QtGui.QLabel('Roznica symetryczna')
        self.dist3Value = QtGui.QLabel('')
        
        self.dist4Label = QtGui.QLabel('Falszywe pozytywy i negatywy')
        self.dist4Value = QtGui.QLabel('')
        
        
        #buttons
        self.path1Btn = QtGui.QPushButton('Wybierz pierwsze drzewo')
        self.path1Btn.clicked.connect(self.chooseFile1)
        
        self.path2Btn = QtGui.QPushButton('Wybierz drugie drzewo')
        self.path2Btn.clicked.connect(self.chooseFile2)
        
        self.calcBtn = QtGui.QPushButton('Oblicz odleglosci')
        self.calcBtn.clicked.connect(self.calculateDistance)
        
        #window layout
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)
        
        self.grid.addWidget(self.path1Label, 1, 0)
        self.grid.addWidget(self.path1Btn, 1, 1)
        
        self.grid.addWidget(self.path2Label, 2, 0)
        self.grid.addWidget(self.path2Btn, 2, 1)

        self.grid.addWidget(self.calcBtn, 3, 1)
        
        self.grid.addWidget(self.dist1Label, 4, 0)
        self.grid.addWidget(self.dist1Value, 4, 1)
        
        self.grid.addWidget(self.dist2Label, 5, 0)
        self.grid.addWidget(self.dist2Value, 5, 1)
        
        self.grid.addWidget(self.dist3Label, 6, 0)
        self.grid.addWidget(self.dist3Value, 6, 1)
        
        self.grid.addWidget(self.dist4Label, 7, 0)
        self.grid.addWidget(self.dist4Value, 7, 1)
        
        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 50)
        self.setWindowTitle('Obliczanie odleglosci')
        self.show()
        
    def showOpenFileDialog(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik', '.newick')
        return str(fname)
    
    
    def chooseFile1(self):
        self.path1 = self.showOpenFileDialog()
        self.path1Label.setText(self.path1)
    
    
    def chooseFile2(self):
        self.path2 = self.showOpenFileDialog()
        self.path2Label.setText(self.path2)
        
    
    def calculateDistance(self):
        if self.path1 != '' and self.path2 != '':
            #get files extensions

            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

            #open tree files
            tns = dendropy.TaxonNamespace()
            self.tree1 = dendropy.Tree.get_from_path(self.path1, self.fileExtension1, taxon_namespace=tns)
            self.tree2 = dendropy.Tree.get_from_path(self.path2, self.fileExtension2, taxon_namespace=tns)

            self.tree1.encode_bipartitions()
            self.tree2.encode_bipartitions()

            print(treecompare.false_positives_and_negatives(self.tree1, self.tree2))

            # self.tree1 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')
            # self.tree2 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')

            # self.tree1.encode_bipartitions()
            #self.tree2.encode_bipartitions()


            #calculate distances
            #self.symDist = self.tree1.symmetric_difference(self.tree2)
            self.symDist = treecompare.symmetric_difference(self.tree1, self.tree2)
            self.fpnDist = treecompare.false_positives_and_negatives(self.tree1, self.tree2)
            self.eucDist = treecompare.euclidean_distance(self.tree1, self.tree2)
            self.rfDist  = treecompare.robinson_foulds_distance(self.tree1, self.tree2)
            
            #show distances
            self.dist1Value.setText(str(self.eucDist))
            self.dist2Value.setText(str(self.rfDist))
            self.dist3Value.setText(str(self.symDist))
            self.dist4Value.setText(str(self.fpnDist))


class NodesWindow(QtGui.QWidget):    
    def __init__(self, tree):
        super(NodesWindow, self).__init__()
        self.tree = tree
        self.terminals = tree.get_terminals()
        self.cb1 = ''
        self.cb2 = ''
        self.initUI()
        
    def initUI(self):
        self.items = { 'nothing':'' }
        #labels
        self.node1Label = QtGui.QLabel('Pierwszy wezel')
        self.node2Label = QtGui.QLabel('Drugi wezel')
        
        #combo boxes
        self.node1ComboBox = QtGui.QComboBox(self)
        self.node2ComboBox = QtGui.QComboBox(self)

        # DODAJE PUSTY WEZEL
        self.node1ComboBox.addItem('')
        self.node2ComboBox.addItem('')


        for clade in self.terminals:
            self.node1ComboBox.addItem(clade.name)
            self.node2ComboBox.addItem(clade.name)


        self.node1ComboBox.activated[str].connect(self.onActivatedCB1)
        self.node2ComboBox.activated[str].connect(self.onActivatedCB2)

        #if self.node1ComboBox.activated[str].connect(self.onActivatedCB1) == '':
        #    print '2'

        self.pathBtn = QtGui.QPushButton('Pokaz sciezke')
        self.pathBtn.clicked.connect(self.showPathWindow)
        
        
        #window layout
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)
        
        self.grid.addWidget(self.node1Label, 1, 0)
        self.grid.addWidget(self.node2Label, 1, 1)
        
        self.grid.addWidget(self.node1ComboBox, 2, 0)
        self.grid.addWidget(self.node2ComboBox, 2, 1)

        self.grid.addWidget(self.pathBtn, 3, 1)        
        
        self.setLayout(self.grid)
        self.setGeometry(200, 200, 100, 50)
        self.setWindowTitle('Wybor wezlow')
        self.show()
        #self.showPathWindow()
    
    
    def onActivatedCB1(self, text):
        self.cb1 = str(text)
        
        
    def onActivatedCB2(self, text):
        self.cb2 = str(text)
        
        
    def showPathWindow(self):
        if self.cb1 != '' and self.cb2 != '':            
            self.start = self.tree.find_clades(self.cb1).next()
            self.end = self.tree.find_clades(name = self.cb2).next()
            
            for clade in self.tree.trace(self.start, self.end):
                clade.color = 'red'
                
            for clade in self.tree.find_clades():
                if not(clade in self.tree.trace(self.start, self.end)):
                    clade.color = 'grey'
            
            self.start.color = 'blue'

            # RYSOWANIE
            Phylo.draw_graphviz(self.tree, node_size = 2500)
            pylab.plot(0,0)

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
            img = pylab.imread('wally.png', 'rb')
            pylab.imshow(img)
            pylab.plot(0,0)

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


class InfoWindow(QtGui.QWidget):    
    def __init__(self, tree):
        super(InfoWindow, self).__init__()
        self.tree = tree
        self.initUI()
        
    def initUI(self):
        self.getMinMaxDepths()
        self.getMinMaxLengths()
        self.numberOfTermianls = self.tree.count_terminals()
        self.numberOfNonterminals = len(self.tree.get_nonterminals())
        
        self.isBifurcating = ''
        if self.tree.is_bifurcating():
            self.isBifurcating = 'tak'
        else:
            self.isBifurcating = 'nie'
                
        
        #labels
        self.nodesLabel = QtGui.QLabel('Liczba wezlow w drzewie')
        self.nodesValue = QtGui.QLabel(str(self.numberOfTermianls + self.numberOfNonterminals))
        
        self.leafsLabel = QtGui.QLabel('Liczba lisci w drzewie')
        self.leafsValue = QtGui.QLabel(str(self.numberOfTermianls))
        
        self.nonleafsLabel = QtGui.QLabel('Liczba wezlow niebedacych lisciami')
        self.nonleafsValue = QtGui.QLabel(str(self.numberOfNonterminals))
        
        self.minDepthLabel = QtGui.QLabel('Najmniejsza glebokosc')
        self.minDepthValue = QtGui.QLabel(str(self.minDepth))
        
        self.maxDepthLabel = QtGui.QLabel('Najwieksza glebokosc')
        self.maxDepthValue = QtGui.QLabel(str(self.maxDepth))
        
        self.minLengthLabel = QtGui.QLabel('Najmniejsza odleglosc od korzenia do liscia')
        self.minLengthValue = QtGui.QLabel(str(self.minLength))
        
        self.maxLengthLabel = QtGui.QLabel('Najwieksza odleglosc od korzenia do liscia')
        self.maxLengthValue = QtGui.QLabel(str(self.maxLength))
        
        self.totalLengthLabel = QtGui.QLabel('Laczna dlugosc galezi w drzewie')
        self.totalLengthValue = QtGui.QLabel(str(self.tree.total_branch_length()))
        
        self.bifurcatingLabel = QtGui.QLabel('Drzewo scisle dwudzielne')
        self.bifurcatingValue = QtGui.QLabel(self.isBifurcating)
        
        #window layout
        self.grid = QtGui.QGridLayout(self)
        self.grid.setSpacing(10)
        
        self.grid.addWidget(self.nodesLabel, 1, 0)
        self.grid.addWidget(self.nodesValue, 1, 1)

        self.grid.addWidget(self.leafsLabel, 2, 0)
        self.grid.addWidget(self.leafsValue, 2, 1)

        self.grid.addWidget(self.nonleafsLabel, 3, 0)
        self.grid.addWidget(self.nonleafsValue, 3, 1)
        
        self.grid.addWidget(self.maxDepthLabel, 4, 0)
        self.grid.addWidget(self.maxDepthValue, 4, 1)
        
        self.grid.addWidget(self.minDepthLabel, 5, 0)
        self.grid.addWidget(self.minDepthValue, 5, 1)
        
        self.grid.addWidget(self.maxLengthLabel, 6, 0)
        self.grid.addWidget(self.maxLengthValue, 6, 1)
        
        self.grid.addWidget(self.minLengthLabel, 7, 0)
        self.grid.addWidget(self.minLengthValue, 7, 1)
        
        self.grid.addWidget(self.totalLengthLabel, 8, 0)
        self.grid.addWidget(self.totalLengthValue, 8, 1)
        
        self.grid.addWidget(self.bifurcatingLabel, 9, 0)
        self.grid.addWidget(self.bifurcatingValue, 9, 1)
        
        
        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 200)
        self.setWindowTitle('Informacje')
        self.show()
        
    def getMinMaxDepths(self):
        self.minDepth = -1
        self.maxDepth = -1
        for clade in self.tree.depths(unit_branch_lengths=True):
            if not(clade.name is None):
                self.currentDepth = self.tree.depths(unit_branch_lengths=True)[clade]
                if self.currentDepth < self.minDepth or self.minDepth == -1:
                    self.minDepth = int(self.currentDepth)
                
                if self.currentDepth > self.maxDepth or self.maxDepth == -1:
                    self.maxDepth = int(self.currentDepth)
                    
                    
    def getMinMaxLengths(self):
        self.minLength = -1.0
        self.maxLength = -1.0
        for clade in self.tree.depths():
            if not(clade.name is None):
                self.currentLength = self.tree.depths()[clade]
                if self.currentLength < self.minLength or self.minLength == -1.0:
                    self.minLength = self.currentLength
                
                if self.currentLength > self.maxLength or self.maxLength == -1.0:
                    self.maxLength = self.currentLength
        

class AsciiTreeWindow(QtGui.QWidget):
    tree = 0
    
    def __init__(self, tree):
        super(AsciiTreeWindow, self).__init__()
        self.tree = tree
        self.initUI()
        
    def initUI(self):        
        #field for drawing Ascii tree
        self.textEdit = QtGui.QTextEdit()
        self.textEdit.setReadOnly(True)
        self.textEdit.setFontFamily('Courier')
        self.textEdit.setWordWrapMode(True)
        #self.textEdit.setStyleSheet('')
        
        # layout
        self.layout = QtGui.QVBoxLayout(self)
        self.layout.addWidget(self.textEdit)
        self.setLayout(self.layout)
        
        #print tree
        self.tmpf = open('/tmp/ascii.txt', 'w')
        Phylo.draw_ascii(self.tree, self.tmpf)
        
        self.tmpf = open('/tmp/ascii.txt', 'r')
        with self.tmpf:        
                self.data = self.tmpf.read()
                self.textEdit.setText(self.data)
        
        self.setGeometry(200, 200, 700, 400)
        self.setWindowTitle('Tekstowe wyswietlanie')
        self.show()
        

class ConvertWindow(QtGui.QWidget):
    chosenFileName = ''
    chosenInputFormat = ''
    chosenOutputFormat = ''
    
    def __init__(self):
        super(ConvertWindow, self).__init__()
        self.initUI()
        
    def initUI(self):
        #labels
        self.inputLabel = QtGui.QLabel('Wejscie')
        self.outputLabel = QtGui.QLabel('Wyjscie')
        
        #fields for files content
        self.textEditIn = QtGui.QTextEdit()
        self.textEditIn.setReadOnly(True)
        self.textEditOut = QtGui.QTextEdit()
        self.textEditOut.setReadOnly(True)
        
        #format groups
        self.inputFormatGroup = self.createInputFormatExclusiveGroup()
        self.outpuFormatGroup = self.createOutputFormatExclusiveGroup()
        
        #button for opening file to convert
        self.openBtn = QtGui.QPushButton('Wybierz plik')
        self.openBtn.clicked.connect(lambda : self.showOpenFileDialog(self.textEditIn))
        self.convBtn = QtGui.QPushButton('Konwertuj')
        self.convBtn.clicked.connect(lambda : self.convertTreeFile(self.textEditIn, self.textEditOut))
        
        #window layout
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
        
        self.radio1.clicked.connect(lambda : self.rememberInputFormat(self.radio1.text().toLower()))
        self.radio2.clicked.connect(lambda : self.rememberInputFormat(self.radio2.text().toLower()))
        self.radio3.clicked.connect(lambda : self.rememberInputFormat(self.radio3.text().toLower()))

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

        self.radio1.clicked.connect(lambda : self.rememberOutputFormat(self.radio1.text().toLower()))
        self.radio2.clicked.connect(lambda : self.rememberOutputFormat(self.radio2.text().toLower()))
        self.radio3.clicked.connect(lambda : self.rememberOutputFormat(self.radio3.text().toLower()))

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
        
        # convert
        if self.chosenInputFormat != self.chosenOutputFormat:
            if self.chosenFileName != '' and self.chosenInputFormat != '':
                
                self.convertedFileName = str(self.chosenFileName).replace(
                    '.' + str(self.chosenInputFormat),'.' + str(self.chosenOutputFormat))
                
                Phylo.convert(str(self.chosenFileName), str(self.chosenInputFormat),
                              self.convertedFileName, str(self.chosenOutputFormat))
                
                f = open(self.convertedFileName, 'r')
        
                with f:        
                    data = f.read()
                    outputTextEdit.setText(data)


class MainWindow(QtGui.QMainWindow):
    chosenFileName = ''
    tree = 0
    
    def __init__(self):
        super(MainWindow, self).__init__()
        
        self.initUI()
        
    def initUI(self):      
        self.textEdit = QtGui.QTextEdit()
        self.setCentralWidget(self.textEdit)
        self.statusBar()

        #menu bar
        self.openFile = QtGui.QAction(QtGui.QIcon(''), 'Otworz', self)
        self.openFile.setShortcut('Ctrl+O')
        self.openFile.setStatusTip('Otworz plik')
        self.openFile.triggered.connect(self.showOpenFileDialog)
        
        self.convertFile = QtGui.QAction(QtGui.QIcon(''), 'Konwertuj', self)
        self.convertFile.setShortcut('Ctrl+W')
        self.convertFile.setStatusTip('Przekonwertuj plik na inny format')
        self.convertFile.triggered.connect(self.showConvertWindow)
        
        self.drawAsciiTree = QtGui.QAction(QtGui.QIcon(''), 'Tekstowo - ASCII', self)
        self.drawAsciiTree.setShortcut('Ctrl+T')
        self.drawAsciiTree.setStatusTip('Rysuj drzewo tekstowo - ASCII')
        self.drawAsciiTree.triggered.connect(self.showAsciiTreeWindow)
        
        self.matplotlibDrawTree = QtGui.QAction(QtGui.QIcon(''), 'Graficznie - matplotlib', self)
        self.matplotlibDrawTree.setShortcut('Ctrl+M')
        self.matplotlibDrawTree.setStatusTip('Rysuj ukorzenione drzewo - matplotlib')
        self.matplotlibDrawTree.triggered.connect(self.showMatplotlibTreeWindow)
        
        self.graphvizDrawTree = QtGui.QAction(QtGui.QIcon(''), 'Graficznie - graphviz', self)
        self.graphvizDrawTree.setShortcut('Ctrl+G')
        self.graphvizDrawTree.setStatusTip('Rysuj nieukorzenione drzewo - graphviz')
        self.graphvizDrawTree.triggered.connect(self.showGraphvizRootedTreeWindow)
        
        self.pathNodesTree = QtGui.QAction(QtGui.QIcon(''), 'Znajdz sciezke', self)
        self.pathNodesTree.setShortcut('Ctrl+P')
        self.pathNodesTree.setStatusTip('Znajduje sciezke miedzy dwoma wezlami w drzewie')
        self.pathNodesTree.triggered.connect(self.showChooseNodesWindow)
        
        self.distanceTree = QtGui.QAction(QtGui.QIcon(''), 'Wyznacz odleglosci', self)
        self.distanceTree.setShortcut('Ctrl+D')
        self.distanceTree.setStatusTip('Oblicza rozne metryki odleglosci miedzy dwoma drzewami')
        self.distanceTree.triggered.connect(self.showDistancesWindow)
        
        self.consensusTree = QtGui.QAction(QtGui.QIcon(''), 'Wyznacz drzewo konsensusu', self)
        self.consensusTree.setShortcut('Ctrl+S')
        self.consensusTree.setStatusTip('Wyznacza drzewa konsensusu za pomoca Dendropy i Bio.Nexus')
        self.consensusTree.triggered.connect(self.showConsensusWindow)
        
        self.treeInfo = QtGui.QAction(QtGui.QIcon(''), 'Informacje', self)
        self.treeInfo.setShortcut('Ctrl+I')
        self.treeInfo.setStatusTip('Pokaz podstawowe informacje o wczytanym drzewie')
        self.treeInfo.triggered.connect(self.showTreeInfoWindow)

        self.menubar = self.menuBar()
        self.fileMenu = self.menubar.addMenu('&Plik')
        self.fileMenu.addAction(self.openFile)
        self.fileMenu.addAction(self.convertFile)
        self.drawMenu = self.menubar.addMenu('&Rysuj')
        self.drawMenu.addAction(self.drawAsciiTree)
        self.drawMenu.addAction(self.matplotlibDrawTree)
        self.drawMenu.addAction(self.graphvizDrawTree)
        self.toolMenu = self.menubar.addMenu('&Narzedzia')
        self.toolMenu.addAction(self.pathNodesTree)
        self.toolMenu.addAction(self.distanceTree)
        self.toolMenu.addAction(self.consensusTree)
        self.otherMenu = self.menubar.addMenu('&Inne')
        self.otherMenu.addAction(self.treeInfo)
        
        
        self.setGeometry(200, 200, 550, 400)
        self.setWindowTitle('Drzewa filogenetyczne')
        self.show()
        
        
    def showConvertWindow(self):
        self.convWin = ConvertWindow()
        self.convWin.show()
        
        
    def showAsciiTreeWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()
        
        if self.tree != 0:
            self.asciiWin = AsciiTreeWindow(self.tree)
            self.asciiWin.show()
            
            
    def showMatplotlibTreeWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()
        
        if self.tree != 0:
            self.tree.root.color = '#808080'
            Phylo.draw(self.tree, branch_labels = lambda c: c.branch_length)
            
    ## showGraphvizUnRootedTreeWindow
    def showGraphvizRootedTreeWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()
        
        if self.tree != 0:
            #self.tree.rooted = True
            self.tree.root.color = '#808080'
            Phylo.draw_graphviz(self.tree, node_size = 2500)
            pylab.show()
        
        
    def showTreeInfoWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()
        
        if self.tree != 0:
            self.infoWin = InfoWindow(self.tree)
            self.infoWin.show()
            
    
    def showChooseNodesWindow(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()
            
        if self.tree != 0:
            self.nodesWin = NodesWindow(self.tree)
            self.nodesWin.show()
            
            
    def showDistancesWindow(self):
        self.distancesWin = DistancesWindow()
        self.distancesWin.show()
    
    
    def showConsensusWindow(self):
        self.consWin = ConsensusWindow()
        self.consWin.show()
        

    def showOpenFileDialog(self):
          
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Otworz plik', './Trees')
        
        if fname != '':
            self.chosenFileName = str(fname)
            
            f = open(str(fname), 'r')
        
            with f:        
                data = f.read()
                self.textEdit.setText(data)
                
            self.fileExtension = (os.path.splitext(str(fname))[1])[1:]
            
            self.tree = Phylo.read(str(fname), self.fileExtension)
                                
         
def main():
    app = QtGui.QApplication(sys.argv)
    mw = MainWindow()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()