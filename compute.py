"""
compute.py - Computational functions for Spatial Viewer
"""

from PyQt5.QtCore import QThread, pyqtSignal
import scanpy as sc
import traceback


class UMAPWorker(QThread):
    """Worker thread for UMAP computation"""

    finished = pyqtSignal(object)  # Emits adata when done
    error = pyqtSignal(str)  # Emits error message
    progress = pyqtSignal(str)  # Emits progress messages

    def __init__(self, adata, n_components=2, n_neighbors=15, min_dist=0.5):
        super().__init__()
        self.adata = adata.copy()  # Work on a copy
        self.n_components = n_components
        self.n_neighbors = n_neighbors
        self.min_dist = min_dist

    def run(self):
        """Run UMAP computation"""
        try:
            # Step 1: Check if neighbors are computed properly
            # Need to check both .uns and .obsp for proper neighbor graph
            needs_neighbors = False

            if "neighbors" not in self.adata.uns:
                needs_neighbors = True
            elif "connectivities" not in self.adata.obsp:
                needs_neighbors = True
            elif "distances" not in self.adata.obsp:
                needs_neighbors = True

            if needs_neighbors:
                self.progress.emit("Computing neighbor graph...")
                sc.pp.neighbors(self.adata, n_neighbors=self.n_neighbors)

            # Step 2: Compute UMAP
            self.progress.emit(f"Computing UMAP ({self.n_components}D)...")
            sc.tl.umap(self.adata, n_components=self.n_components, min_dist=self.min_dist)

            self.progress.emit("UMAP completed!")
            self.finished.emit(self.adata)

        except Exception as e:
            error_details = traceback.format_exc()
            self.error.emit(f"UMAP computation failed:\n{str(e)}\n\nDetails:\n{error_details}")


class ComputeManager:
    """Manages computational tasks"""

    @staticmethod
    def has_umap(adata):
        """Check if UMAP coordinates exist"""
        if adata is None:
            return False

        umap_keys = ["X_umap", "umap"]
        for key in umap_keys:
            if key in adata.obsm:
                return True

        # Check for partial matches
        for key in adata.obsm.keys():
            if "umap" in key.lower():
                return True

        return False

    @staticmethod
    def get_umap_dimensions(adata):
        """Get number of UMAP dimensions if exists"""
        if not ComputeManager.has_umap(adata):
            return 0

        umap_keys = ["X_umap", "umap"]
        for key in umap_keys:
            if key in adata.obsm:
                return adata.obsm[key].shape[1]

        # Check for partial matches
        for key in adata.obsm.keys():
            if "umap" in key.lower():
                return adata.obsm[key].shape[1]

        return 0

    @staticmethod
    def has_neighbors(adata):
        """Check if neighbor graph exists"""
        if adata is None:
            return False

        # Check both .uns and .obsp
        has_uns = "neighbors" in adata.uns
        has_connectivities = "connectivities" in adata.obsp
        has_distances = "distances" in adata.obsp

        return has_uns and has_connectivities and has_distances
