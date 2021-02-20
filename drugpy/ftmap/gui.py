import traceback
from os.path import expanduser

import pymol.gui
import pymol.Qt
from pymol import cmd as pm

from .core import load_atlas, load_ftmap

QFileDialog = pymol.Qt.QtWidgets.QFileDialog
QFormLayout = pymol.Qt.QtWidgets.QFormLayout
QPushButton = pymol.Qt.QtWidgets.QPushButton
QLineEdit = pymol.Qt.QtWidgets.QLineEdit
QDialog = pymol.Qt.QtWidgets.QDialog
QRegExp = pymol.Qt.QtCore.QRegExp
QtCore = pymol.Qt.QtCore
QRegExpValidator = pymol.Qt.QtGui.QRegExpValidator


class LoadFTMapServerResultDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Setup window
        self.setModal(True)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)

        self.layout = QFormLayout(self)

        # Input id
        self.result_id_line = QLineEdit("", self)
        self.result_id_line.setValidator(QRegExpValidator(QRegExp("^[0-9]+$")))
        self.layout.addRow("FTMap result id:", self.result_id_line)

        # Ok / Cancel buttons
        self.fetch_button = QPushButton("Fetch")
        self.layout.addWidget(self.fetch_button)
        self.fetch_button.clicked.connect(self.start)

    def start(self):
        result_id = self.result_id_line.text()
        if not result_id:
            return
        self.fetch_button.setDisabled(True)
        try:
            ret = load_ftmap(result_id)
        except Exception as exc:
            self.done(QDialog.Rejected)
            raise Exception("The result cannot be loaded.")
        self.done(QDialog.Accepted)


def init_plugin_gui(menu):
    load_server_result_action = menu.addAction("Load FTMap server result")

    @load_server_result_action.triggered.connect
    def triggered():
        LoadFTMapServerResultDialog().exec_()

    load_atlas_pdb_action = menu.addAction("Load Atlas PDB")

    @load_atlas_pdb_action.triggered.connect
    def triggered():
        atlas_pdb = QFileDialog.getOpenFileNames(
            menu, "Atlas result file", expanduser("~"), "PDB file (*.pdb)"
        )[0]
        load_atlas(*atlas_pdb)

    load_ftmap_pdb_action = menu.addAction("Load FTMap PDB")

    @load_ftmap_pdb_action.triggered.connect
    def triggered():
        ftmap_pdb = QFileDialog.getOpenFileNames(
            menu, "FTMap result file", expanduser("~"), "PDB file (*.pdb)"
        )[0]
        load_ftmap(*ftmap_pdb)
