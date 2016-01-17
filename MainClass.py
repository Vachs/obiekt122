#!/usr/bin/env python
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


class ConsensusWindow(QtGui.QWidget):
    def __init__(self):
        super(ConsensusWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.initUI()

    def initUI(self):
        # sciezki
        self.path1Label = QtGui.QLabel(self.path1)
        self.path2Label = QtGui.QLabel(self.path2)

        # buttony
        self.path1Btn = QtGui.QPushButton('Wybierz pierwsze drzewo')
        self.path1Btn.clicked.connect(self.chooseFile1)

        self.path2Btn = QtGui.QPushButton('Wybierz drugie drzewo')
        self.path2Btn.clicked.connect(self.chooseFile2)

        self.consDenBtn = QtGui.QPushButton('Wygeneruj drzewo konsensusu')
        self.consDenBtn.clicked.connect(self.drawConsensusTreeBioNexus)

        self.consNexBtn = QtGui.QPushButton('Pokaz informacje o tym drzewie')
        self.consNexBtn.clicked.connect(self.showConsensusTreeBioNexus)

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
        self.grid.addWidget(self.path1Label, 3, 0)
        self.grid.addWidget(self.path2Label, 4, 0)


        self.grid.addWidget(self.consDenBtn, 5, 0)
        self.grid.addWidget(self.consNexBtn, 6, 0)
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

    def chooseFile1(self):
        self.path1 = self.showOpenFileDialog()
        self.path1Label.setText(self.path1)

    def chooseFile2(self):
        self.path2 = self.showOpenFileDialog()
        self.path2Label.setText(self.path2)

    def OpenInfoWindow(self):
        if self.tree != 0:
            self.infoWin = InfoWindow(self.tree)
            self.infoWin.show()

    def showOtherRF(self):
        if self.path1 != '' and self.path2 != '':
            self.calculateDistance()



    def showConsensusTreeBioNexus(self):
        if self.path1 != '' and self.path2 != '':
            # get files extensions
            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

            # open tree files
            self.trees = []

            # first tree
            self.f = open(self.path1, 'r')
            self.tree1 = Trees.Tree(self.f.read())
            self.trees.append(self.tree1)
            self.f.close()

            # second tree
            self.f = open(self.path2, 'r')
            self.tree2 = Trees.Tree(self.f.read())
            self.trees.append(self.tree2)
            self.f.close()

            # generate consensus tree
            self.consensus_tree = Trees.consensus(self.trees)

            # draw tree
            self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            #self.tree.root.color = '#808080'
            self.OpenInfoWindow()
            #Phylo.draw(self.tree)

    def drawConsensusTreeBioNexus(self):
        if self.path1 != '' and self.path2 != '':
            #self.calculateDistance()
            # get files extensions
            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

            # open tree files
            self.trees = []

            # first tree
            self.f = open(self.path1, 'r')
            self.tree1 = Trees.Tree(self.f.read())
            self.trees.append(self.tree1)
            self.f.close()

            # second tree
            self.f = open(self.path2, 'r')
            self.tree2 = Trees.Tree(self.f.read())
            self.trees.append(self.tree2)
            self.f.close()

            # generate consensus tree
            self.consensus_tree = Trees.consensus(self.trees)

            # draw tree
            self.handle = StringIO(self.consensus_tree.to_string(plain_newick=True))
            self.tree = Phylo.read(self.handle, 'newick')
            self.tree.root.color = '#808080'
            #self.OpenInfoWindow()
            Phylo.draw(self.tree)

    def calculateDistance(self):
        if self.path1 != '' and self.path2 != '':
            # get files extensions

            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

            # open tree files
            tns = dendropy.TaxonNamespace()
            self.tree1 = dendropy.Tree.get_from_path(self.path1, self.fileExtension1, taxon_namespace=tns)
            self.tree2 = dendropy.Tree.get_from_path(self.path2, self.fileExtension2, taxon_namespace=tns)

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

    # ROOTuje drzewo
    def makeRootUnroot(self, mod):
        if self.path1 != '' and self.path2 == '':
            # get files extensions
            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

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


# Okno wyboru kolejnych opcji. Wybierz co chcesz zrobic dalej
class NextWindow(QtGui.QWidget):
    chosenFileName = ''
    tree = 0
    path1 = ''

    def __init__(self):
        super(NextWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        # self.tree = tree
        # self.terminals = tree.get_terminals()
        self.cb1 = ''
        self.cb2 = ''
        self.initUI()

    def initUI(self):
        self.showOpenFileDialog()

        #Trees.Tree.display(self.tree)
        # WYSZUIKWANIE DROGI PO KRAWEDZIACH
        # labels
        self.node1Label = QtGui.QLabel('Pierwszy wezel')
        self.node2Label = QtGui.QLabel('Drugi wezel')

        # combo boxes
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

        self.node1ComboBox.activated[str].connect(self.onActivatedCB1)
        self.node2ComboBox.activated[str].connect(self.onActivatedCB2)

        self.pathBtn = QtGui.QPushButton('wyznacz')
        self.pathBtn.clicked.connect(self.showPathWindow)

        # buttony
        self.info = QtGui.QLabel('Kalkulator drzew filogenicznych. Mozliwosc wczytania zarowno drzew ukorzenionych\n'
                                 'jak i nieukorzenionych. Wybierz jedna z opcji dostepnych ponizej.')

        self.lay0 = QtGui.QLabel('>> OPERACJA NA JEDNYM DRZEWIE')
        self.lay1 = QtGui.QLabel('Pozwala uzyskac podstawowe informacje o wybranych drzewie')
        self.btn1 = QtGui.QPushButton('Uzyskaj informacje o drzewie')
        # self.btn1.whatsThis()
        # self.btn1.setWhatsThis('sadsa')
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
        self.anotherTree.clicked.connect(self.showOpenFileDialog)

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
        self.grid.addWidget(self.node1Label, 9, 1)
        self.grid.addWidget(self.node2Label, 9, 2)
        self.grid.addWidget(self.pathBtn, 10, 0)
        self.grid.addWidget(self.node1ComboBox, 10, 1)
        self.grid.addWidget(self.node2ComboBox, 10, 2)

        self.grid.addWidget(self.lay4, 11, 0)
        self.grid.addWidget(self.lay5, 12, 0)
        self.grid.addWidget(self.btn4, 13, 0)

        # self.grid.addWidget(self.consDenBtn, 3, 0)
        # self.grid.addWidget(self.consNexBtn, 3, 1)
        # self.grid.addWidget(self.consNexBtn2, 3, 2)

        self.setLayout(self.grid)
        self.setGeometry(200, 200, 200, 50)
        self.setWindowTitle('Menu')
        self.show()

    def showConsensusWindow(self):
        self.consWin = ConsensusWindow()
        self.consWin.show()


    def showOpenFileDialog(self):
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
            self.showOpenFileDialog()

        if self.tree != 0:
            self.infoWin = InfoWindow(self.tree)
            self.infoWin.show()

    # RYSUJE PODSTAWOWY GRAF
    def DrawSimple(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()

        if self.tree != 0:
            self.tree.root.color = '#808080'
            Phylo.draw_graphviz(self.tree, node_size=2500)
            pylab.show()

    # RYSUJE BARDZIEJ ZAAWANSOWANY WYKRES
    def DrawAdvan(self):
        if self.chosenFileName == '':
            self.showOpenFileDialog()

        if self.tree != 0:
            self.tree.root.color = '#808080'
            Phylo.draw(self.tree, branch_labels=lambda c: c.branch_length)

    def onActivatedCB1(self, text):
        self.cb1 = str(text)

    def onActivatedCB2(self, text):
        self.cb2 = str(text)

    def showPathWindow(self):
        if self.cb1 != '' and self.cb2 != '':
            self.start = self.tree.find_clades(self.cb1).next()
            self.end = self.tree.find_clades(name=self.cb2).next()

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
            img = pylab.imread('wally.png', 'rb')
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

# WCIELIC TO DO INFORMACJI

class DistancesWindow(QtGui.QWidget):
    def __init__(self):
        super(DistancesWindow, self).__init__()
        self.path1 = ''
        self.path2 = ''
        self.initUI()

    def initUI(self):
        # labels
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

        # buttons
        self.path1Btn = QtGui.QPushButton('Wybierz pierwsze drzewo')
        self.path1Btn.clicked.connect(self.chooseFile1)

        self.path2Btn = QtGui.QPushButton('Wybierz drugie drzewo')
        self.path2Btn.clicked.connect(self.chooseFile2)

        self.calcBtn = QtGui.QPushButton('Oblicz odleglosci')
        self.calcBtn.clicked.connect(self.calculateDistance)

        # window layout
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
            # get files extensions

            self.fileExtension1 = (os.path.splitext(self.path1)[1])[1:]
            self.fileExtension2 = (os.path.splitext(self.path2)[1])[1:]

            # open tree files
            tns = dendropy.TaxonNamespace()
            self.tree1 = dendropy.Tree.get_from_path(self.path1, self.fileExtension1, taxon_namespace=tns)
            self.tree2 = dendropy.Tree.get_from_path(self.path2, self.fileExtension2, taxon_namespace=tns)

            self.tree1.encode_bipartitions()
            self.tree2.encode_bipartitions()

            print(treecompare.false_positives_and_negatives(self.tree1, self.tree2))

            # self.tree1 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')
            # self.tree2 = dendropy.Tree.get_from_string('((A, B), (C, D))', 'newick')

            # self.tree1.encode_bipartitions()
            # self.tree2.encode_bipartitions()


            # calculate distances
            # self.symDist = self.tree1.symmetric_difference(self.tree2)
            self.symDist = treecompare.symmetric_difference(self.tree1, self.tree2)
            self.fpnDist = treecompare.false_positives_and_negatives(self.tree1, self.tree2)
            self.eucDist = treecompare.euclidean_distance(self.tree1, self.tree2)
            self.rfDist = treecompare.robinson_foulds_distance(self.tree1, self.tree2)

            # show distances
            self.dist1Value.setText(str(self.eucDist))
            self.dist2Value.setText(str(self.rfDist))
            self.dist3Value.setText(str(self.symDist))
            self.dist4Value.setText(str(self.fpnDist))


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


def main():
    app = QtGui.QApplication(sys.argv)
    # mw = MainWindow()
    wd = NextWindow()
    sys.exit(app.exec_())

## TODO:
## potrzeba dodac jeszcze usuwania danego wierzcholka
if __name__ == '__main__':
    main()
