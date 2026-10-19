#!/usr/bin/env python3
# Author: Khatereh Kashmari
"""Read LAMMPS ave/chunk output and plot mean density profiles along each axis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# Defaults for running without command-line arguments.

input_DIR = Path(
    "/Users/khatereh/Desktop/Research/PE/Lammps_Replicate_Method/PCFF_Class2/rep2/300K"
)
file_name = "PEdens_rep2_densityprofile"
output_DIR = "/Users/khatereh/Desktop/Research/Polymer_CoScientist/generate_uff_lammps"


def read_ave_chunk(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return fixed bin coordinates, density blocks, and their timesteps."""
    if not path.is_file():
        raise FileNotFoundError(f"Density-profile file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    coordinates: list[np.ndarray] = []
    densities: list[np.ndarray] = []
    timesteps: list[int] = []
    index = 0

    while index < len(lines):
        text = lines[index].strip()
        index += 1
        if not text or text.startswith("#"):
            continue
        fields = text.split()
        if len(fields) != 3:
            continue
        try:
            timestep, number_of_chunks = int(fields[0]), int(fields[1])
            float(fields[2])
        except ValueError:
            continue

        block_coordinates: list[float] = []
        block_densities: list[float] = []
        for _ in range(number_of_chunks):
            if index >= len(lines):
                raise ValueError(f"Incomplete timestep block in {path}")
            row = lines[index].split()
            index += 1
            if len(row) < 4:
                raise ValueError(f"Malformed chunk row in {path}: {' '.join(row)}")
            block_coordinates.append(float(row[1]))
            block_densities.append(float(row[3]))

        coordinates.append(np.asarray(block_coordinates))
        densities.append(np.asarray(block_densities))
        timesteps.append(timestep)

    if not densities:
        raise ValueError(f"No ave/chunk density blocks found in {path}")
    if any(len(block) != len(densities[0]) for block in densities):
        raise ValueError(f"The number of spatial bins changes in {path}")
    if not all(np.allclose(block, coordinates[0]) for block in coordinates[1:]):
        raise ValueError(f"Bin coordinates change between timestep blocks in {path}")

    return coordinates[0], np.vstack(densities), np.asarray(timesteps)


def plot_profiles(input_directory: Path, prefix: str, output_directory: Path) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    colors = {"x": "#D55E00", "y": "#0072B2", "z": "#009E73"}
    line_styles = {"x": "-", "y": "-", "z": "-"}
    profiles: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for axis in "xyz":
        path = input_directory / f"{prefix}_{axis}.txt"
        profiles[axis] = read_ave_chunk(path)

    plot_styles = {
        "font.family": "Times New Roman",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "axes.labelweight": "bold",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "lines.linewidth": 1.25,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
    plt.rcParams.update(plot_styles)
    figure, axes = plt.subplots(
        3, 1, figsize=(7.1, 6.6), sharex=True, sharey=True, constrained_layout=True
    )

    for panel, (subplot, axis) in enumerate(zip(axes, "xyz")):
        coordinate, density_blocks, timesteps = profiles[axis]
        mean = density_blocks.mean(axis=0)
        subplot.plot(
            coordinate,
            mean,
            color=colors[axis],
            linestyle=line_styles[axis],
            label=f"{axis.upper()} direction",
        )
        subplot.axhline(
            mean.mean(),
            color="0.25",
            linestyle="-",
            linewidth=0.9,
            label=f"Spatial mean: {mean.mean():.3f} g cm$^{{-3}}$",
        )
        subplot.set_ylabel(r"Density (g cm$^{-3}$)")
        subplot.text(
            0.012,
            0.93,
            f"({chr(97 + panel)})",
            transform=subplot.transAxes,
            fontweight="bold",
            va="top",
        )
        subplot.text(
            0.985,
            0.93,
            f"{axis.upper()} direction",
            transform=subplot.transAxes,
            ha="right",
            va="top",
        )
        subplot.tick_params(which="both", direction="in", top=True, right=True)
        subplot.tick_params(which="major", length=4, width=0.8)
        subplot.tick_params(which="minor", length=2, width=0.6)
        subplot.xaxis.set_minor_locator(AutoMinorLocator(2))
        subplot.yaxis.set_minor_locator(AutoMinorLocator(2))
        subplot.grid(False)
        subplot.legend(
            loc="lower right", frameon=False, handlelength=2.5, borderaxespad=0.3
        )
    axes[-1].set_xlabel(r"Position ($\mathrm{\AA}$)")

    panel_stem = output_directory / f"{prefix}_xyz_panels"
    for extension, dpi in (("pdf", None), ("png", 600)):
        options = (
            {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
        )
        figure.savefig(
            panel_stem.with_suffix(f".{extension}"),
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
            **options,
        )
    plt.close(figure)
    figure, axis_plot = plt.subplots(figsize=(3.5, 2.75), constrained_layout=True)
    for axis in "xyz":
        coordinate, density_blocks, _ = profiles[axis]
        axis_plot.plot(
            coordinate,
            density_blocks.mean(axis=0),
            color=colors[axis],
            linestyle=line_styles[axis],
            label=axis.upper(),
        )
    axis_plot.set_xlabel(r"Position ($\mathrm{\AA}$)")
    axis_plot.set_ylabel(r"Density (g cm$^{-3}$)")
    axis_plot.tick_params(which="both", direction="in", top=True, right=True)
    axis_plot.tick_params(which="major", length=4, width=0.8)
    axis_plot.tick_params(which="minor", length=2, width=0.6)
    axis_plot.xaxis.set_minor_locator(AutoMinorLocator(2))
    axis_plot.yaxis.set_minor_locator(AutoMinorLocator(2))
    axis_plot.legend(
        frameon=False, ncol=3, loc="upper center", handlelength=2.2, columnspacing=1.2
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input-directory", type=Path, default=input_DIR)
    parser.add_argument("--prefix", default=file_name)
    parser.add_argument("-o", "--output-directory", type=Path, default=output_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    plot_profiles(
        args.input_directory.expanduser().resolve(),
        args.prefix,
        args.output_directory.expanduser().resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
