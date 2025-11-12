"""
progress_dialog.py - Progress dialog with spinner
"""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QProgressBar
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from PyQt5.QtGui import QFont


class SpinnerLabel(QLabel):
    """Custom spinning animation label"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.angle = 0
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.rotate)
        self.setFixedSize(60, 60)
        self.setAlignment(Qt.AlignCenter)
        self.update_spinner()

    def start(self):
        """Start spinning animation"""
        self.timer.start(50)  # Update every 50ms

    def stop(self):
        """Stop spinning animation"""
        self.timer.stop()

    def rotate(self):
        """Rotate the spinner"""
        self.angle = (self.angle + 15) % 360
        self.update_spinner()

    def update_spinner(self):
        """Update spinner display"""
        # Create a simple rotating circle with dots
        dots = "●" * 8
        rotation_index = (self.angle // 45) % 8

        # Create visual effect
        display = ""
        for i in range(8):
            if i == rotation_index:
                display += "◉"
            elif abs(i - rotation_index) == 1 or abs(i - rotation_index) == 7:
                display += "◍"
            else:
                display += "○"

        self.setText(display[0])
        self.setStyleSheet(f"""
            QLabel {{
                font-size: 48px;
                color: #2196F3;
                transform: rotate({self.angle}deg);
            }}
        """)


class ProgressDialog(QDialog):
    """Progress dialog with spinner and messages"""

    cancelled = pyqtSignal()

    def __init__(self, parent=None, dark_mode=False):
        super().__init__(parent)
        self.dark_mode = dark_mode
        self.setWindowTitle("Computing UMAP")
        self.setModal(True)
        self.setFixedSize(400, 250)
        self.setWindowFlags(Qt.Dialog | Qt.CustomizeWindowHint | Qt.WindowTitleHint)

        self.setup_ui()
        self.apply_theme()

    def setup_ui(self):
        """Setup the UI"""
        layout = QVBoxLayout(self)
        layout.setSpacing(20)
        layout.setContentsMargins(30, 30, 30, 30)

        # Title
        self.title_label = QLabel("🔬 Computing UMAP")
        self.title_label.setFont(QFont("Arial", 16, QFont.Bold))
        self.title_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.title_label)

        # Spinner
        self.spinner = SpinnerLabel(self)
        layout.addWidget(self.spinner, alignment=Qt.AlignCenter)

        # Progress bar (indeterminate)
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(0)  # Indeterminate mode
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setFixedHeight(8)
        layout.addWidget(self.progress_bar)

        # Status message
        self.status_label = QLabel("Initializing...")
        self.status_label.setFont(QFont("Arial", 11))
        self.status_label.setAlignment(Qt.AlignCenter)
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Info label
        self.info_label = QLabel("This may take a few moments depending on dataset size")
        self.info_label.setFont(QFont("Arial", 9))
        self.info_label.setAlignment(Qt.AlignCenter)
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        layout.addStretch()

    def apply_theme(self):
        """Apply theme to dialog"""
        if self.dark_mode:
            self.setStyleSheet("""
                QDialog {
                    background-color: #2b2b2b;
                    color: #e0e0e0;
                }
                QLabel {
                    color: #e0e0e0;
                }
                QProgressBar {
                    border: 2px solid #424242;
                    border-radius: 4px;
                    background-color: #1e1e1e;
                }
                QProgressBar::chunk {
                    background-color: #2196F3;
                    border-radius: 2px;
                }
            """)
            self.info_label.setStyleSheet("color: #b0b0b0;")
        else:
            self.setStyleSheet("""
                QDialog {
                    background-color: white;
                    color: #333;
                }
                QLabel {
                    color: #333;
                }
                QProgressBar {
                    border: 2px solid #e0e0e0;
                    border-radius: 4px;
                    background-color: #f5f5f5;
                }
                QProgressBar::chunk {
                    background-color: #2196F3;
                    border-radius: 2px;
                }
            """)
            self.info_label.setStyleSheet("color: #666;")

    def update_status(self, message):
        """Update status message"""
        self.status_label.setText(message)

    def showEvent(self, event):
        """Start spinner when dialog is shown"""
        super().showEvent(event)
        self.spinner.start()

    def closeEvent(self, event):
        """Stop spinner when dialog is closed"""
        self.spinner.stop()
        super().closeEvent(event)
