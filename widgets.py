"""
widgets.py - Custom widgets for Spatial Viewer
"""

from PyQt5.QtWidgets import QPushButton


class ModernButton(QPushButton):
    """Custom styled button with theme support"""

    def __init__(self, text, primary=False, parent=None):
        super().__init__(text)
        self.primary = primary
        self.parent_widget = parent
        self.setMinimumHeight(36)
        self.apply_theme()

    def apply_theme(self, dark_mode=False):
        """Apply theme-specific styles"""
        if self.primary:
            if dark_mode:
                self.setStyleSheet("""
                    QPushButton {
                        background-color: #1976D2;
                        color: white;
                        border: none;
                        border-radius: 6px;
                        padding: 8px 16px;
                        font-weight: bold;
                        font-size: 13px;
                    }
                    QPushButton:hover {
                        background-color: #2196F3;
                    }
                    QPushButton:pressed {
                        background-color: #0D47A1;
                    }
                    QPushButton:disabled {
                        background-color: #424242;
                        color: #757575;
                    }
                """)
            else:
                self.setStyleSheet("""
                    QPushButton {
                        background-color: #2196F3;
                        color: white;
                        border: none;
                        border-radius: 6px;
                        padding: 8px 16px;
                        font-weight: bold;
                        font-size: 13px;
                    }
                    QPushButton:hover {
                        background-color: #1976D2;
                    }
                    QPushButton:pressed {
                        background-color: #0D47A1;
                    }
                    QPushButton:disabled {
                        background-color: #e0e0e0;
                        color: #9e9e9e;
                    }
                """)
        else:
            if dark_mode:
                self.setStyleSheet("""
                    QPushButton {
                        background-color: #424242;
                        color: #e0e0e0;
                        border: 1px solid #616161;
                        border-radius: 6px;
                        padding: 8px 16px;
                        font-size: 13px;
                    }
                    QPushButton:hover {
                        background-color: #525252;
                        border-color: #757575;
                    }
                    QPushButton:pressed {
                        background-color: #323232;
                    }
                """)
            else:
                self.setStyleSheet("""
                    QPushButton {
                        background-color: #f5f5f5;
                        color: #333;
                        border: 1px solid #ddd;
                        border-radius: 6px;
                        padding: 8px 16px;
                        font-size: 13px;
                    }
                    QPushButton:hover {
                        background-color: #e0e0e0;
                        border-color: #bbb;
                    }
                    QPushButton:pressed {
                        background-color: #d0d0d0;
                    }
                """)
