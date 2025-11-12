"""
main.py - Main application for Spatial RNA-seq Viewer
"""

import sys
import traceback
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFileDialog,
    QLabel,
    QMessageBox,
    QComboBox,
    QGroupBox,
    QSplitter,
    QRadioButton,
    QButtonGroup,
    QSlider,
    QCheckBox,
    QLineEdit,
    QCompleter,
    QSpinBox,
    QDoubleSpinBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import anndata as ad

# Import custom modules
from widgets import ModernButton
from themes import ThemeManager
from plotting import PlotManager
from compute import UMAPWorker, ComputeManager
from progress_dialog import ProgressDialog


class SpatialViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.adata = None
        self.current_view = "spatial"
        self.current_color = None
        self.dark_mode = False
        self.buttons = []
        self.theme_manager = ThemeManager()
        self.plot_manager = PlotManager()
        self.compute_manager = ComputeManager()
        self.umap_worker = None
        self.progress_dialog = None
        self.initUI()

    def initUI(self):
        self.setWindowTitle("🔬 Spatial RNA-seq Viewer Pro")
        self.setGeometry(50, 50, 1400, 900)

        # Main widget and splitter
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        self.splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(self.splitter)

        # Create panels
        self.left_panel = self.create_control_panel()
        self.splitter.addWidget(self.left_panel)

        self.right_panel = self.create_visualization_panel()
        self.splitter.addWidget(self.right_panel)

        # Set splitter proportions
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 3)
        self.splitter.setSizes([350, 1050])

        # Apply initial theme
        self.apply_theme()

    def create_control_panel(self):
        """Create left control panel"""
        panel = QWidget()
        panel.setMaximumWidth(400)

        layout = QVBoxLayout(panel)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(15)

        # Title with theme toggle
        title_layout = QHBoxLayout()

        title = QLabel("🔬 Spatial Viewer")
        title.setFont(QFont("Arial", 18, QFont.Bold))
        title_layout.addWidget(title)

        title_layout.addStretch()

        # Theme toggle button
        self.theme_button = ModernButton("🌙", parent=self)
        self.theme_button.setMaximumWidth(50)
        self.theme_button.setToolTip("Toggle dark mode")
        self.theme_button.clicked.connect(self.toggle_theme)
        self.buttons.append(self.theme_button)
        title_layout.addWidget(self.theme_button)

        layout.addLayout(title_layout)

        # Separator line
        separator = QLabel()
        separator.setFixedHeight(2)
        layout.addWidget(separator)
        layout.addSpacing(5)

        # File loading group
        file_group = QGroupBox("Data Loading")
        file_group.setFont(QFont("Arial", 14, QFont.Bold))
        file_layout = QVBoxLayout()

        self.load_button = ModernButton("📁 Load AnnData File (.h5ad)", primary=True, parent=self)
        self.load_button.clicked.connect(self.load_file)
        self.buttons.append(self.load_button)
        file_layout.addWidget(self.load_button)

        self.info_label = QLabel("No data loaded")
        self.info_label.setWordWrap(True)
        self.info_label.setFont(QFont("Arial", 12))
        file_layout.addWidget(self.info_label)

        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        # View selection group
        view_group = QGroupBox("View Mode")
        view_group.setFont(QFont("Arial", 14, QFont.Bold))
        view_layout = QVBoxLayout()

        self.view_button_group = QButtonGroup()
        self.spatial_radio = QRadioButton("📍 Spatial Distribution")
        self.umap_2d_radio = QRadioButton("📊 UMAP 2D Projection")
        self.umap_3d_radio = QRadioButton("🎲 UMAP 3D Projection")

        for radio in [self.spatial_radio, self.umap_2d_radio, self.umap_3d_radio]:
            radio.setFont(QFont("Arial", 12))

        self.spatial_radio.setChecked(True)
        self.view_button_group.addButton(self.spatial_radio)
        self.view_button_group.addButton(self.umap_2d_radio)
        self.view_button_group.addButton(self.umap_3d_radio)

        self.spatial_radio.toggled.connect(lambda: self.change_view("spatial"))
        self.umap_2d_radio.toggled.connect(lambda: self.change_view("umap_2d"))
        self.umap_3d_radio.toggled.connect(lambda: self.change_view("umap_3d"))

        for radio in [self.spatial_radio, self.umap_2d_radio, self.umap_3d_radio]:
            radio.setEnabled(False)
            view_layout.addWidget(radio)

        # Compute UMAP button
        self.compute_umap_button = ModernButton("⚡ Compute UMAP", parent=self)
        self.compute_umap_button.clicked.connect(self.show_umap_dialog)
        self.compute_umap_button.setEnabled(False)
        self.buttons.append(self.compute_umap_button)
        view_layout.addWidget(self.compute_umap_button)

        view_group.setLayout(view_layout)
        layout.addWidget(view_group)

        # Color by group
        color_group = QGroupBox("Color By")
        color_group.setFont(QFont("Arial", 14, QFont.Bold))
        color_layout = QVBoxLayout()

        # Gene search
        gene_label = QLabel("🧬 Gene Expression")
        gene_label.setFont(QFont("Arial", 12, QFont.Bold))
        color_layout.addWidget(gene_label)

        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("Type gene name...")
        self.gene_search.setMinimumHeight(36)
        self.gene_search.setFont(QFont("Arial", 12))
        self.gene_search.returnPressed.connect(self.apply_gene_coloring)
        color_layout.addWidget(self.gene_search)

        # Gene completer
        self.gene_completer = QCompleter()
        self.gene_completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.gene_completer.setFilterMode(Qt.MatchContains)
        self.gene_search.setCompleter(self.gene_completer)

        apply_gene_btn = ModernButton("✓ Apply Gene", parent=self)
        apply_gene_btn.clicked.connect(self.apply_gene_coloring)
        self.buttons.append(apply_gene_btn)
        color_layout.addWidget(apply_gene_btn)

        # Separator
        separator2 = QLabel("─" * 35)
        separator2.setAlignment(Qt.AlignCenter)
        color_layout.addWidget(separator2)

        # Obs metadata
        obs_label = QLabel("📊 Cell Metadata")
        obs_label.setFont(QFont("Arial", 12, QFont.Bold))
        color_layout.addWidget(obs_label)

        self.obs_combo = QComboBox()
        self.obs_combo.setMinimumHeight(36)
        self.obs_combo.setFont(QFont("Arial", 12))
        self.obs_combo.currentTextChanged.connect(self.apply_obs_coloring)
        color_layout.addWidget(self.obs_combo)

        # Clear button
        clear_btn = ModernButton("✖ Clear Coloring", parent=self)
        clear_btn.clicked.connect(self.clear_coloring)
        self.buttons.append(clear_btn)
        color_layout.addWidget(clear_btn)

        color_group.setLayout(color_layout)
        layout.addWidget(color_group)

        # Colormap group
        cmap_group = QGroupBox("Color Scheme")
        cmap_group.setFont(QFont("Arial", 14, QFont.Bold))
        cmap_layout = QVBoxLayout()

        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(
            ["viridis", "plasma", "inferno", "magma", "cividis", "RdYlBu_r", "RdBu_r", "Spectral", "coolwarm", "seismic"]
        )
        self.cmap_combo.setMinimumHeight(36)
        self.cmap_combo.setFont(QFont("Arial", 12))
        self.cmap_combo.currentTextChanged.connect(self.update_plot)
        cmap_layout.addWidget(self.cmap_combo)

        cmap_group.setLayout(cmap_layout)
        layout.addWidget(cmap_group)

        # Visualization options
        options_group = QGroupBox("Display Options")
        options_group.setFont(QFont("Arial", 14, QFont.Bold))
        options_layout = QVBoxLayout()

        self.show_colorbar = QCheckBox("Show color scale")
        self.show_colorbar.setFont(QFont("Arial", 12))
        self.show_colorbar.setChecked(True)
        self.show_colorbar.stateChanged.connect(self.update_plot)
        options_layout.addWidget(self.show_colorbar)

        # Point size slider
        size_label = QLabel("Point Size")
        size_label.setFont(QFont("Arial", 12, QFont.Bold))
        options_layout.addWidget(size_label)

        self.size_slider = QSlider(Qt.Horizontal)
        self.size_slider.setMinimum(1)
        self.size_slider.setMaximum(100)
        self.size_slider.setValue(10)
        self.size_slider.valueChanged.connect(self.update_plot)
        options_layout.addWidget(self.size_slider)

        self.size_value_label = QLabel("Current: 10")
        self.size_value_label.setFont(QFont("Arial", 9))
        self.size_slider.valueChanged.connect(lambda v: self.size_value_label.setText(f"Current: {v}"))
        options_layout.addWidget(self.size_value_label)

        options_group.setLayout(options_layout)
        layout.addWidget(options_group)

        # Export button
        self.export_button = ModernButton("💾 Export Plot", parent=self)
        self.export_button.clicked.connect(self.export_plot)
        self.export_button.setEnabled(False)
        self.buttons.append(self.export_button)
        layout.addWidget(self.export_button)

        layout.addStretch()
        return panel

    def create_visualization_panel(self):
        """Create right visualization panel"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(15, 15, 15, 15)

        # Matplotlib figure
        self.figure = Figure(figsize=(12, 10))
        self.canvas = FigureCanvas(self.figure)

        # Navigation toolbar
        self.toolbar = NavigationToolbar(self.canvas, panel)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        return panel

    def toggle_theme(self):
        """Toggle between light and dark theme"""
        self.dark_mode = not self.dark_mode
        self.theme_button.setText("☀️" if self.dark_mode else "🌙")
        self.apply_theme()

    def apply_theme(self):
        """Apply the current theme to all UI elements"""
        if self.dark_mode:
            self.setStyleSheet(self.theme_manager.get_dark_theme())
            self.figure.patch.set_facecolor("#2b2b2b")
            self.left_panel.setStyleSheet(self.theme_manager.get_panel_style(True))
            self.info_label.setStyleSheet(self.theme_manager.get_info_label_style(True))
        else:
            self.setStyleSheet(self.theme_manager.get_light_theme())
            self.figure.patch.set_facecolor("white")
            self.left_panel.setStyleSheet(self.theme_manager.get_panel_style(False))
            self.info_label.setStyleSheet(self.theme_manager.get_info_label_style(False))

        # Update all buttons
        for button in self.buttons:
            button.apply_theme(self.dark_mode)

        # Redraw plot with new background
        if self.adata is not None:
            self.update_plot()
        else:
            self.canvas.draw()

    def load_file(self):
        """Open file dialog and load AnnData file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select AnnData File", "", "AnnData Files (*.h5ad);;All Files (*)"
        )

        if file_path:
            try:
                self.info_label.setText("⏳ Loading data...")
                QApplication.processEvents()

                self.adata = ad.read_h5ad(file_path)

                # Update info
                n_obs = self.adata.n_obs
                n_vars = self.adata.n_vars
                file_name = file_path.split("/")[-1]
                self.info_label.setText(f"✅ {file_name}\n📊 {n_obs:,} cells × {n_vars:,} genes")

                # Enable controls
                self.spatial_radio.setEnabled(True)
                self.umap_2d_radio.setEnabled(True)
                self.umap_3d_radio.setEnabled(True)
                self.export_button.setEnabled(True)
                self.compute_umap_button.setEnabled(True)

                # Populate color options
                self.populate_gene_list()
                self.populate_obs_list()

                # Check UMAP availability
                self.check_umap_availability()

                # Plot
                self.update_plot()

            except Exception as e:
                error_msg = f"Error loading file:\n{str(e)}\n\n{traceback.format_exc()}"
                QMessageBox.critical(self, "Error", error_msg)
                self.info_label.setText("❌ Failed to load file")

    def populate_gene_list(self):
        """Populate gene search with autocomplete"""
        if self.adata.var_names is not None:
            gene_list = list(self.adata.var_names)
            from PyQt5.QtCore import QStringListModel

            model = QStringListModel(gene_list)
            self.gene_completer.setModel(model)

    def populate_obs_list(self):
        """Populate obs dropdown with available metadata columns"""
        self.obs_combo.clear()
        self.obs_combo.addItem("None")

        if self.adata.obs is not None:
            for col in self.adata.obs.columns:
                self.obs_combo.addItem(col)

    def apply_gene_coloring(self):
        """Apply coloring by gene expression"""
        gene_name = self.gene_search.text().strip()

        if not gene_name:
            QMessageBox.warning(self, "Warning", "Please enter a gene name")
            return

        if gene_name not in self.adata.var_names:
            QMessageBox.warning(
                self,
                "Gene Not Found",
                f"Gene '{gene_name}' not found in dataset.\nTotal genes: {len(self.adata.var_names)}",
            )
            return

        self.current_color = ("gene", gene_name)
        self.obs_combo.setCurrentIndex(0)
        self.update_plot()

    def apply_obs_coloring(self, text):
        """Apply coloring by obs metadata"""
        if text == "None":
            if self.current_color and self.current_color[0] != "gene":
                self.current_color = None
                self.update_plot()
        else:
            self.current_color = ("obs", text)
            self.gene_search.clear()
            self.update_plot()

    def clear_coloring(self):
        """Clear all coloring"""
        self.current_color = None
        self.gene_search.clear()
        self.obs_combo.setCurrentIndex(0)
        self.update_plot()

    def change_view(self, view):
        """Change visualization mode"""
        # Check if UMAP is available for UMAP views
        if view in ["umap_2d", "umap_3d"]:
            if not self.compute_manager.has_umap(self.adata):
                reply = QMessageBox.question(
                    self,
                    "UMAP Not Found",
                    "UMAP coordinates not found in the data.\n\nWould you like to compute UMAP now?",
                    QMessageBox.Yes | QMessageBox.No,
                )
                if reply == QMessageBox.Yes:
                    self.show_umap_dialog()
                    return
                else:
                    # Revert to spatial view
                    self.spatial_radio.setChecked(True)
                    return

            # Check dimensions for 3D
            if view == "umap_3d":
                dims = self.compute_manager.get_umap_dimensions(self.adata)
                if dims < 3:
                    reply = QMessageBox.question(
                        self,
                        "UMAP 3D Not Available",
                        f"Current UMAP has only {dims} dimensions.\n\nWould you like to compute UMAP 3D?",
                        QMessageBox.Yes | QMessageBox.No,
                    )
                    if reply == QMessageBox.Yes:
                        self.show_umap_dialog(n_components=3)
                        return
                    else:
                        # Revert to previous view
                        if self.current_view == "spatial":
                            self.spatial_radio.setChecked(True)
                        else:
                            self.umap_2d_radio.setChecked(True)
                        return

        self.current_view = view
        self.update_plot()

    def check_umap_availability(self):
        """Check and update UMAP availability status"""
        if not self.compute_manager.has_umap(self.adata):
            self.compute_umap_button.setText("⚡ Compute UMAP")
        else:
            dims = self.compute_manager.get_umap_dimensions(self.adata)
            self.compute_umap_button.setText(f"⚡ Recompute UMAP ({dims}D)")

    def show_umap_dialog(self, n_components=None):
        """Show dialog to configure and run UMAP"""
        dialog = UMAPConfigDialog(self, self.dark_mode, n_components)
        if dialog.exec_() == QDialog.Accepted:
            params = dialog.get_parameters()
            self.run_umap(**params)

    def run_umap(self, n_components=2, n_neighbors=15, min_dist=0.5):
        """Run UMAP computation in background thread"""
        if self.adata is None:
            return

        # Create progress dialog
        self.progress_dialog = ProgressDialog(self, self.dark_mode)

        # Create worker thread
        self.umap_worker = UMAPWorker(self.adata, n_components, n_neighbors, min_dist)

        # Connect signals
        self.umap_worker.progress.connect(self.progress_dialog.update_status)
        self.umap_worker.finished.connect(self.on_umap_finished)
        self.umap_worker.error.connect(self.on_umap_error)

        # Start computation
        self.umap_worker.start()
        self.progress_dialog.show()

    def on_umap_finished(self, adata):
        """Handle UMAP computation completion"""
        # Update adata with new UMAP coordinates
        self.adata.obsm.update(adata.obsm)
        self.adata.uns.update(adata.uns)

        # Close progress dialog
        if self.progress_dialog:
            self.progress_dialog.close()
            self.progress_dialog = None

        # Update UI
        self.check_umap_availability()

        # Show success message
        dims = self.compute_manager.get_umap_dimensions(self.adata)
        QMessageBox.information(
            self,
            "Success",
            f"UMAP computation completed!\n\nGenerated {dims}D UMAP coordinates for {self.adata.n_obs:,} cells.",
        )

        # Update plot if in UMAP view
        if self.current_view in ["umap_2d", "umap_3d"]:
            self.update_plot()

    def on_umap_error(self, error_msg):
        """Handle UMAP computation error"""
        # Close progress dialog
        if self.progress_dialog:
            self.progress_dialog.close()
            self.progress_dialog = None

        # Show error
        QMessageBox.critical(self, "UMAP Error", error_msg)

    def get_color_values(self):
        """Get values for coloring"""
        if self.current_color is None:
            return None, None

        color_type, value = self.current_color

        if color_type == "gene":
            if value in self.adata.var_names:
                idx = list(self.adata.var_names).index(value)
                vals = (
                    self.adata.X[:, idx].toarray().flatten()
                    if hasattr(self.adata.X, "toarray")
                    else self.adata.X[:, idx].flatten()
                )
                return vals, f"{value}"

        elif color_type == "obs":
            if value in self.adata.obs.columns:
                vals = self.adata.obs[value]
                if vals.dtype == "object" or vals.dtype.name == "category":
                    return vals.astype("category").cat.codes.values, f"{value}"
                return vals.values, f"{value}"

        return None, None

    def update_plot(self):
        """Update the current plot"""
        if self.adata is None:
            return

        try:
            color_vals, color_label = self.get_color_values()
            size = self.size_slider.value()
            cmap = self.cmap_combo.currentText()
            show_colorbar = self.show_colorbar.isChecked()
            plot_colors = self.theme_manager.get_plot_colors(self.dark_mode)

            if self.current_view == "spatial":
                self.plot_manager.plot_spatial(
                    self.figure, self.adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors
                )
            elif self.current_view == "umap_2d":
                self.plot_manager.plot_umap_2d(
                    self.figure, self.adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors
                )
            elif self.current_view == "umap_3d":
                self.plot_manager.plot_umap_3d(
                    self.figure, self.adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors
                )

            self.canvas.draw()

        except Exception as e:
            QMessageBox.warning(self, "Plot Error", f"Error updating plot:\n{str(e)}")

    def export_plot(self):
        """Export current plot to file"""
        if self.adata is None:
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg);;All Files (*)"
        )

        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches="tight", facecolor=self.figure.get_facecolor())
                QMessageBox.information(self, "Success", f"Plot saved successfully:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error saving plot:\n{str(e)}")


class UMAPConfigDialog(QDialog):
    """Dialog for configuring UMAP parameters"""

    def __init__(self, parent=None, dark_mode=False, default_components=None):
        super().__init__(parent)
        self.dark_mode = dark_mode
        self.setWindowTitle("⚙️ UMAP Configuration")
        self.setModal(True)
        self.setFixedWidth(400)

        self.setup_ui(default_components)
        self.apply_theme()

    def setup_ui(self, default_components):
        """Setup the UI"""
        layout = QVBoxLayout(self)
        layout.setSpacing(15)
        layout.setContentsMargins(20, 20, 20, 20)

        # Title
        title = QLabel("Configure UMAP Parameters")
        title.setFont(QFont("Arial", 14, QFont.Bold))
        layout.addWidget(title)

        # Description
        desc = QLabel("Adjust parameters for UMAP computation:")
        desc.setFont(QFont("Arial", 12))
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Form layout
        form_layout = QFormLayout()
        form_layout.setSpacing(10)

        # Number of components
        self.n_components_spin = QSpinBox()
        self.n_components_spin.setRange(2, 3)
        self.n_components_spin.setValue(default_components if default_components else 2)
        self.n_components_spin.setFont(QFont("Arial", 12))
        form_layout.addRow("Dimensions:", self.n_components_spin)

        # Number of neighbors
        self.n_neighbors_spin = QSpinBox()
        self.n_neighbors_spin.setRange(5, 100)
        self.n_neighbors_spin.setValue(15)
        self.n_neighbors_spin.setFont(QFont("Arial", 12))
        form_layout.addRow("Neighbors:", self.n_neighbors_spin)

        # Minimum distance
        self.min_dist_spin = QDoubleSpinBox()
        self.min_dist_spin.setRange(0.0, 1.0)
        self.min_dist_spin.setSingleStep(0.1)
        self.min_dist_spin.setValue(0.5)
        self.min_dist_spin.setFont(QFont("Arial", 12))
        form_layout.addRow("Min Distance:", self.min_dist_spin)

        layout.addLayout(form_layout)

        # Info text
        info = QLabel(
            "• Dimensions: 2D for visualization, 3D for advanced analysis\n"
            "• Neighbors: Higher values preserve global structure\n"
            "• Min Distance: Controls point clustering (lower = tighter)"
        )
        info.setFont(QFont("Arial", 9))
        info.setWordWrap(True)
        layout.addWidget(info)

        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

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
                QSpinBox, QDoubleSpinBox {
                    background-color: #1e1e1e;
                    color: #e0e0e0;
                    border: 1px solid #424242;
                    border-radius: 4px;
                    padding: 5px;
                }
                QSpinBox:focus, QDoubleSpinBox:focus {
                    border-color: #2196F3;
                }
                QPushButton {
                    background-color: #424242;
                    color: #e0e0e0;
                    border: 1px solid #616161;
                    border-radius: 4px;
                    padding: 8px 16px;
                    min-width: 80px;
                }
                QPushButton:hover {
                    background-color: #525252;
                }
            """)
        else:
            self.setStyleSheet("""
                QDialog {
                    background-color: white;
                    color: #333;
                }
                QLabel {
                    color: #333;
                }
                QSpinBox, QDoubleSpinBox {
                    background-color: white;
                    color: #333;
                    border: 1px solid #ddd;
                    border-radius: 4px;
                    padding: 5px;
                }
                QSpinBox:focus, QDoubleSpinBox:focus {
                    border-color: #2196F3;
                }
                QPushButton {
                    background-color: #f5f5f5;
                    color: #333;
                    border: 1px solid #ddd;
                    border-radius: 4px;
                    padding: 8px 16px;
                    min-width: 80px;
                }
                QPushButton:hover {
                    background-color: #e0e0e0;
                }
            """)

    def get_parameters(self):
        """Get selected parameters"""
        return {
            "n_components": self.n_components_spin.value(),
            "n_neighbors": self.n_neighbors_spin.value(),
            "min_dist": self.min_dist_spin.value(),
        }


def main():
    app = QApplication(sys.argv)

    # Set application-wide font (system-appropriate)
    import platform

    if platform.system() == "Darwin":  # macOS
        font = QFont("SF Pro Text", 12)
    elif platform.system() == "Windows":
        font = QFont("Segoe UI", 12)
    else:  # Linux
        font = QFont("Ubuntu", 12)
    app.setFont(font)

    viewer = SpatialViewer()
    viewer.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
