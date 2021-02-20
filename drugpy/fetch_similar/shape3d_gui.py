import pymol.gui
import pymol.Qt

from .shape3d import fetch_similar_shape3d

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

        self.similarity_spin = QDoubleSpinBox(self)
        self.similarity_spin.setRange(0.0, 1.0)
        self.similarity_spin.setValue(0.5)
        self.layout.addRow("Min similarity:", self.similarity_spin)

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
        self.max_resol_spin.setValue(2.0)
        self.layout.addRow("Max resolution:", self.max_resol_spin)

        # Fetch button
        self.fetch_button = QPushButton("Fetch")
        self.layout.addWidget(self.fetch_button)
        self.fetch_button.clicked.connect(self.start)

    def start(self):
        pdb_id = self.ref_pdb_line.text()

        min_similarity = self.similarity_spin.value()
        max_resolution = self.max_resol_spin.value()
        ligand = self.ligand_line.text().strip() or None
        dist = self.dist_spin.value()
        compounds = self.compounds_line.text()
        prosthetic_groups = self.prosthetics_line.text()

        fetch_similar_shape3d(
            pdb_id,
            min_similarity,
            max_resolution,
            ligand,
            dist,
            compounds,
            prosthetic_groups,
        )

        self.done(QDialog.Accepted)


def init_gui(menu):

    action = menu.addAction("Fetch similar (Shape 3D)")

    @action.triggered.connect
    def triggered():
        dialog = FetchBlastClusterDialog()
        dialog.exec_()
