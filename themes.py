"""
themes.py - Theme management for Spatial Viewer
"""


class ThemeManager:
    """Manages light and dark themes"""

    @staticmethod
    def get_light_theme():
        """Returns light theme stylesheet"""
        return """
            QMainWindow {
                background-color: #fafafa;
            }
            QWidget {
                background-color: #fafafa;
                color: #333;
            }
            QLabel {
                color: #333;
            }
            QGroupBox {
                font-weight: bold;
                border: 2px solid #e0e0e0;
                border-radius: 8px;
                margin-top: 12px;
                padding-top: 12px;
                background-color: white;
                color: #333;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #1976D2;
            }
            QRadioButton {
                color: #333;
                spacing: 8px;
            }
            QRadioButton::indicator {
                width: 18px;
                height: 18px;
            }
            QCheckBox {
                color: #333;
                spacing: 8px;
            }
            QLineEdit {
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 5px 10px;
                background-color: white;
                color: #333;
                font-size: 12px;
            }
            QLineEdit:focus {
                border-color: #2196F3;
                border-width: 2px;
            }
            QComboBox {
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 5px;
                background-color: white;
                color: #333;
            }
            QComboBox:hover {
                border-color: #2196F3;
            }
            QComboBox::drop-down {
                border: none;
                width: 20px;
            }
            QComboBox::down-arrow {
                image: url(none);
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #666;
                margin-right: 5px;
            }
            QComboBox QAbstractItemView {
                background-color: white;
                color: #333;
                selection-background-color: #2196F3;
                border: 1px solid #ddd;
            }
            QSlider::groove:horizontal {
                height: 6px;
                background: #e0e0e0;
                border-radius: 3px;
            }
            QSlider::handle:horizontal {
                background: #2196F3;
                width: 16px;
                margin: -5px 0;
                border-radius: 8px;
            }
            QSlider::handle:horizontal:hover {
                background: #42A5F5;
            }
            QToolBar {
                background-color: #f5f5f5;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 4px;
            }
            QSplitter::handle {
                background-color: #e0e0e0;
            }
        """

    @staticmethod
    def get_dark_theme():
        """Returns dark theme stylesheet"""
        return """
            QMainWindow {
                background-color: #1e1e1e;
            }
            QWidget {
                background-color: #1e1e1e;
                color: #e0e0e0;
            }
            QLabel {
                color: #e0e0e0;
            }
            QGroupBox {
                font-weight: bold;
                border: 2px solid #424242;
                border-radius: 8px;
                margin-top: 12px;
                padding-top: 12px;
                background-color: #2b2b2b;
                color: #e0e0e0;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
                color: #2196F3;
            }
            QRadioButton {
                color: #e0e0e0;
                spacing: 8px;
            }
            QRadioButton::indicator {
                width: 18px;
                height: 18px;
            }
            QCheckBox {
                color: #e0e0e0;
                spacing: 8px;
            }
            QLineEdit {
                border: 1px solid #424242;
                border-radius: 4px;
                padding: 5px 10px;
                background-color: #2b2b2b;
                color: #e0e0e0;
                font-size: 12px;
            }
            QLineEdit:focus {
                border-color: #2196F3;
                border-width: 2px;
            }
            QComboBox {
                border: 1px solid #424242;
                border-radius: 4px;
                padding: 5px;
                background-color: #2b2b2b;
                color: #e0e0e0;
            }
            QComboBox:hover {
                border-color: #2196F3;
            }
            QComboBox::drop-down {
                border: none;
                width: 20px;
            }
            QComboBox::down-arrow {
                image: url(none);
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #e0e0e0;
                margin-right: 5px;
            }
            QComboBox QAbstractItemView {
                background-color: #2b2b2b;
                color: #e0e0e0;
                selection-background-color: #2196F3;
                border: 1px solid #424242;
            }
            QSlider::groove:horizontal {
                height: 6px;
                background: #424242;
                border-radius: 3px;
            }
            QSlider::handle:horizontal {
                background: #2196F3;
                width: 16px;
                margin: -5px 0;
                border-radius: 8px;
            }
            QSlider::handle:horizontal:hover {
                background: #42A5F5;
            }
            QToolBar {
                background-color: #2b2b2b;
                border: 1px solid #424242;
                border-radius: 4px;
                padding: 4px;
            }
            QSplitter::handle {
                background-color: #424242;
            }
        """

    @staticmethod
    def get_panel_style(dark_mode):
        """Returns panel-specific styles"""
        if dark_mode:
            return "background-color: #2b2b2b; border-right: 1px solid #424242;"
        return "background-color: white; border-right: 1px solid #e0e0e0;"

    @staticmethod
    def get_info_label_style(dark_mode):
        """Returns info label styles"""
        if dark_mode:
            return "color: #b0b0b0; padding: 8px; font-size: 12px;"
        return "color: #666; padding: 8px; font-size: 12px;"

    @staticmethod
    def get_plot_colors(dark_mode):
        """Returns colors for plot elements"""
        if dark_mode:
            return {"bg_color": "#2b2b2b", "text_color": "#e0e0e0", "spine_color": "#424242"}
        return {"bg_color": "white", "text_color": "#333", "spine_color": "#ddd"}
