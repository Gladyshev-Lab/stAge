"""
plotting.py - Plotting functions for Spatial Viewer
"""


class PlotManager:
    """Manages plot generation for different view modes"""

    @staticmethod
    def find_key(dict_like, possible_keys):
        """Find first matching key in dictionary"""
        for key in possible_keys:
            if key in dict_like:
                return key

        # Check for partial matches
        for key in dict_like.keys():
            for pkey in possible_keys:
                if pkey.lower() in key.lower():
                    return key

        return None

    @staticmethod
    def plot_spatial(figure, adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors):
        """Plot spatial distribution"""
        figure.clear()
        ax = figure.add_subplot(111)

        # Set background
        ax.set_facecolor(plot_colors["bg_color"])
        text_color = plot_colors["text_color"]

        # Find spatial coordinates
        spatial_key = PlotManager.find_key(adata.obsm, ["spatial", "X_spatial"])

        if spatial_key is None:
            raise ValueError(f"No spatial coordinates found.\nAvailable: {list(adata.obsm.keys())}")

        coords = adata.obsm[spatial_key]

        if color_vals is not None:
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1], c=color_vals, s=size, alpha=0.7, cmap=cmap, edgecolors="none"
            )
            if show_colorbar:
                cbar = figure.colorbar(scatter, ax=ax, label=color_label)
                cbar.ax.yaxis.label.set_color(text_color)
                cbar.ax.tick_params(colors=text_color)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=size, alpha=0.6, c="steelblue", edgecolors="none")

        ax.set_xlabel("X Coordinate", fontsize=12, color=text_color, fontweight="bold")
        ax.set_ylabel("Y Coordinate", fontsize=12, color=text_color, fontweight="bold")
        ax.set_title(
            f"Spatial Distribution • {adata.n_obs:,} cells", fontsize=14, fontweight="bold", color=text_color, pad=15
        )
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.tick_params(colors=text_color)

        for spine in ax.spines.values():
            spine.set_edgecolor(plot_colors["spine_color"])

        figure.tight_layout()

    @staticmethod
    def plot_umap_2d(figure, adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors):
        """Plot UMAP in 2D"""
        figure.clear()
        ax = figure.add_subplot(111)

        # Set background
        ax.set_facecolor(plot_colors["bg_color"])
        text_color = plot_colors["text_color"]

        umap_key = PlotManager.find_key(adata.obsm, ["X_umap", "umap"])

        if umap_key is None:
            ax.text(
                0.5,
                0.5,
                "UMAP coordinates not found\nPlease run UMAP first",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=14,
                color=text_color,
                fontweight="bold",
            )
            return

        coords = adata.obsm[umap_key]

        if color_vals is not None:
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1], c=color_vals, s=size, alpha=0.7, cmap=cmap, edgecolors="none"
            )
            if show_colorbar:
                cbar = figure.colorbar(scatter, ax=ax, label=color_label)
                cbar.ax.yaxis.label.set_color(text_color)
                cbar.ax.tick_params(colors=text_color)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=size, alpha=0.6, c="steelblue", edgecolors="none")

        ax.set_xlabel("UMAP 1", fontsize=12, color=text_color, fontweight="bold")
        ax.set_ylabel("UMAP 2", fontsize=12, color=text_color, fontweight="bold")
        ax.set_title(
            f"UMAP 2D Projection • {adata.n_obs:,} cells", fontsize=14, fontweight="bold", color=text_color, pad=15
        )
        ax.tick_params(colors=text_color)

        for spine in ax.spines.values():
            spine.set_edgecolor(plot_colors["spine_color"])

        figure.tight_layout()

    @staticmethod
    def plot_umap_3d(figure, adata, color_vals, color_label, size, cmap, show_colorbar, plot_colors):
        """Plot UMAP in 3D"""
        figure.clear()
        ax = figure.add_subplot(111, projection="3d")

        # Set background
        ax.set_facecolor(plot_colors["bg_color"])
        text_color = plot_colors["text_color"]

        if plot_colors["bg_color"] == "#2b2b2b":  # Dark mode
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
            ax.xaxis.pane.set_edgecolor(plot_colors["spine_color"])
            ax.yaxis.pane.set_edgecolor(plot_colors["spine_color"])
            ax.zaxis.pane.set_edgecolor(plot_colors["spine_color"])

        umap_key = PlotManager.find_key(adata.obsm, ["X_umap", "umap"])

        if umap_key is None:
            figure.text(
                0.5,
                0.5,
                "UMAP coordinates not found\nPlease run UMAP first",
                ha="center",
                va="center",
                fontsize=14,
                color=text_color,
                fontweight="bold",
            )
            return

        coords = adata.obsm[umap_key]

        if coords.shape[1] < 3:
            figure.text(
                0.5,
                0.5,
                "UMAP 3D requires 3 components\nCurrent UMAP has only 2 dimensions",
                ha="center",
                va="center",
                fontsize=14,
                color=text_color,
                fontweight="bold",
            )
            return

        if color_vals is not None:
            scatter = ax.scatter(
                coords[:, 0], coords[:, 1], coords[:, 2], c=color_vals, s=size, alpha=0.6, cmap=cmap, edgecolors="none"
            )
            if show_colorbar:
                cbar = figure.colorbar(scatter, ax=ax, label=color_label, shrink=0.6)
                cbar.ax.yaxis.label.set_color(text_color)
                cbar.ax.tick_params(colors=text_color)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=size, alpha=0.6, c="steelblue", edgecolors="none")

        ax.set_xlabel("UMAP 1", fontsize=11, color=text_color, fontweight="bold")
        ax.set_ylabel("UMAP 2", fontsize=11, color=text_color, fontweight="bold")
        ax.set_zlabel("UMAP 3", fontsize=11, color=text_color, fontweight="bold")
        ax.set_title(
            f"UMAP 3D Projection • {adata.n_obs:,} cells", fontsize=14, fontweight="bold", color=text_color, pad=20
        )
        ax.tick_params(colors=text_color)

        figure.tight_layout()
