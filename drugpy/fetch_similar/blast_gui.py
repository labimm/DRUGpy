import pymol.gui
import pymol.Qt

from .blast import fetch_similar_blast, fetch_similar_blast_update

QWidget = pymol.Qt.QtWidgets.QWidget
QFileDialog = pymol.Qt.QtWidgets.QFileDialog
QFormLayout = pymol.Qt.QtWidgets.QFormLayout
QPushButton = pymol.Qt.QtWidgets.QPushButton
QSpinBox = pymol.Qt.QtWidgets.QSpinBox
QDoubleSpinBox = pymol.Qt.QtWidgets.QDoubleSpinBox
QDockWidget = pymol.Qt.QtWidgets.QDockWidget
QLineEdit = pymol.Qt.QtWidgets.QLineEdit
QCheckBox = pymol.Qt.QtWidgets.QCheckBox
QApplication = pymol.Qt.QtWidgets.QApplication
QMessageBox = pymol.Qt.QtWidgets.QMessageBox
QVBoxLayout = pymol.Qt.QtWidgets.QVBoxLayout
QTextEdit = pymol.Qt.QtWidgets.QTextEdit
QDialog = pymol.Qt.QtWidgets.QDialog
QDialogButtonBox = pymol.Qt.QtWidgets.QDialogButtonBox
QDesktopWidget = pymol.Qt.QtWidgets.QDesktopWidget
QProgressBar = pymol.Qt.QtWidgets.QProgressBar
QAction = pymol.Qt.QtWidgets.QAction
QComboBox = pymol.Qt.QtWidgets.QComboBox

LeftDockWidgetArea = pymol.Qt.QtCore.Qt.LeftDockWidgetArea
QRegExp = pymol.Qt.QtCore.QRegExp
QtCore = pymol.Qt.QtCore
QThread = pymol.Qt.QtCore.QThread
pyqtSignal = pymol.Qt.QtCore.Signal

QRegExpValidator = pymol.Qt.QtGui.QRegExpValidator
QPalette = pymol.Qt.QtGui.QPalette
QTextDocument = pymol.Qt.QtGui.QTextDocument
QIntValidator = pymol.Qt.QtGui.QIntValidator
QTextCursor = pymol.Qt.QtGui.QTextCursor
QIcon = pymol.Qt.QtGui.QIcon


class FetchBlastClusterDialog(QDialog):
    def __init__(self, *vina_args, parent=None):
        super().__init__(parent)

        # Setup window
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)

        self.layout = QFormLayout(self)

        # Form
        self.ref_pdb_line = QLineEdit("", self)
        self.ref_pdb_line.setValidator(
            QRegExpValidator(QRegExp("[a-zA-Z0-9]{4}"), self)
        )
        self.layout.addRow("Reference PDB:", self.ref_pdb_line)

        self.ref_chain_line = QLineEdit("", self)
        self.ref_chain_line.setValidator(QRegExpValidator(QRegExp("[A-Z]"), self))
        self.layout.addRow("Reference chain:", self.ref_chain_line)

        self.similarity_combo = QComboBox(self)
        self.similarity_combo.addItem("100", 100)
        self.similarity_combo.addItem("95", 95)
        self.similarity_combo.addItem("90", 90)
        self.similarity_combo.addItem("70", 70)
        self.similarity_combo.addItem("50", 50)
        self.similarity_combo.addItem("40", 40)
        self.similarity_combo.addItem("30", 30)
        self.layout.addRow("Similarity threshold:", self.similarity_combo)

        self.ligand_line = QLineEdit("", self)
        self.layout.addRow("Reference ligand:", self.ligand_line)

        self.dist_spin = QDoubleSpinBox(self)
        self.dist_spin.setRange(0.0, 10.0)
        self.dist_spin.setValue(5.0)
        self.layout.addRow("Distance from ligand:", self.dist_spin)

        self.compounds_line = QLineEdit("organic or inorganic", self)
        self.layout.addRow("Ligand compounds:", self.compounds_line)

        self.prosthetics_line = QLineEdit("HEM FAD NAP NDP ADP FMN", self)
        self.layout.addRow("Prosthetic groups:", self.prosthetics_line)

        self.max_resol_spin = QDoubleSpinBox(self)
        self.max_resol_spin.setRange(0.0, 10.0)
        self.max_resol_spin.setValue(0.0)
        self.layout.addRow("Max resolution:", self.max_resol_spin)

        self.max_structs_spin = QSpinBox(self)
        self.max_structs_spin.setRange(0, 1000)
        self.max_structs_spin.setValue(50)
        self.layout.addRow("Max structures:", self.max_structs_spin)

        # Fetch button
        self.fetch_button = QPushButton("Fetch")
        self.layout.addWidget(self.fetch_button)
        self.fetch_button.clicked.connect(self.start)

        # Update DB button
        self.update_button = QPushButton("Update DB")
        self.layout.addWidget(self.update_button)
        self.update_button.clicked.connect(fetch_similar_blast_update)

    def start(self):
        pdb = self.ref_pdb_line.text()
        chain = self.ref_chain_line.text()
        pdb_chain_id = f"{pdb}_{chain}"

        similarity = self.similarity_combo.currentData()
        ligand = self.ligand_line.text().strip() or None
        dist = self.dist_spin.value()
        compounds = self.compounds_line.text()
        prosthetic_groups = self.prosthetics_line.text()
        max_resolution = self.max_resol_spin.value()
        max_structures = self.max_structs_spin.value()

        fetch_similar_blast(
            pdb_chain_id,
            similarity,
            ligand,
            dist,
            compounds,
            prosthetic_groups,
            max_resolution,
            max_structures,
        )

        self.done(QDialog.Accepted)


def init_gui(menu):

    action = menu.addAction("Fetch similar (Blast)")

    @action.triggered.connect
    def triggered():
        dialog = FetchBlastClusterDialog()
        dialog.exec_()
