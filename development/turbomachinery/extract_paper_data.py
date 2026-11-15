"""Extract Figures 1-7 from the DTU PDF's vector paths using Poppler.

This development tool needs pdftocairo on PATH, not a Python PDF dependency.
It never changes production equations or regenerates references during tests.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "paper_reference/curves.json"
PDF_URL = "https://backend.orbit.dtu.dk/ws/portalfiles/portal/406859073/1-s2.0-S0196890425006090-main.pdf"
PDF_HASH = "b66ae80aa6c8e51470149723690ba7372017885face2966226da775adccf5622"
NUMBER = r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?"

# SVG coordinates are in PDF points from the top left. These frames were
# checked against rendered pages and axis ticks, not fitted to model outputs.
FIGURES = {
    1: dict(page=7, frame=[82.628906, 56.171875, 271.867143, 213.871094],
            limits=[0, .15, .55, .95], legend=[135, 160, 268, 211],
            curves=[(81, "aungier", "Aungier baseline", True),
                    (82, "casey", "Casey correlation", False)]),
    2: dict(page=7, frame=[351.316407, 56.171875, 540.554643, 213.871094],
            limits=[0, .15, .45, .85], legend=[355, 170, 472, 211],
            curves=[(428, "proposed", "Proposed work coefficient", True)]),
    3: dict(page=7, frame=[355.007813, 282.335938, 540.695269, 437.074219],
            limits=[0, .15, .675, .925], legend=None,
            curves=[(781, "pr_1.1", "PR = 1.1", True),
                    (782, "pr_3", "PR = 3", True),
                    (783, "pr_5", "PR = 5", True),
                    (784, "pr_10", "PR = 10", True)]),
    4: dict(page=8, frame=[82.628906, 59.988281, 271.867143, 217.6875],
            limits=[0, .15, .4, .7], legend=[87, 64, 203, 105],
            curves=[(98, "proposed", "Published velocity-ratio curve", True)]),
    5: dict(page=8, frame=[351.316407, 59.355469, 540.554643, 217.054688],
            limits=[0, .15, 1, 5], legend=[355, 63, 474, 105],
            curves=[(420, "proposed", "Proposed diameter/chord ratio", True)]),
    6: dict(page=9, frame=[356.21875, 112.21875, 547.54683, 271.660156],
            limits=[.25, .95, .775, .975], legend=[456, 225, 544, 269],
            curves=[(129, "kofskey", "Kofskey correlation", False),
                    (130, "aungier", "Aungier correlation", False),
                    (131, "proposed", "Proposed total-to-total", True)]),
    7: dict(page=10, frame=[356.21875, 56.175781, 547.54683, 215.617188],
            limits=[.25, .95, .775, .975], legend=[456, 59, 544, 103],
            curves=[(122, "kofskey", "Kofskey correlation", False),
                    (123, "aungier", "Aungier correlation", False),
                    (124, "proposed", "Proposed total-to-static", True)]),
}


def svg_paths(filename):
    """Read visible SVG M/L/C/Z paths and transform their coordinates to the page.

    Glyph definitions are excluded. Curve endpoints and all control points
    are retained separately so polylines and marker centers can be recovered
    without treating Bezier control points as measured observations.
    """
    root = ET.parse(filename).getroot()
    records = []

    def visit(element, matrix):
        """Traverse drawing groups, composing local affine transforms in order."""
        if element.tag.rsplit("}", 1)[-1] == "defs":
            return
        transform = element.get("transform", "")
        if transform:
            if not transform.startswith("matrix("):
                raise ValueError(f"Unsupported SVG transform: {transform}")
            a, b, c, d, e, f = map(float, re.findall(NUMBER, transform))
            matrix = matrix @ np.array([[a, c, e], [b, d, f], [0, 0, 1]])
        if element.tag.endswith("path"):
            tokens = re.findall(r"[MLCZmlcz]|" + NUMBER, element.get("d", ""))
            vertices, endpoints, commands = [], [], []
            cursor = 0
            while cursor < len(tokens):
                command = tokens[cursor]
                cursor += 1
                commands.append(command)
                count = {"M": 2, "L": 2, "C": 6, "Z": 0}[command]
                values = list(map(float, tokens[cursor:cursor + count]))
                cursor += count
                for offset in range(0, count, 2):
                    vertices.append((matrix @ [values[offset], values[offset + 1], 1])[:2])
                if count:
                    endpoints.append(vertices[-1])
            if vertices:
                array = np.asarray(vertices)
                records.append(dict(style=element.get("style", ""), commands="".join(commands),
                                    bounds=[*array.min(axis=0), *array.max(axis=0)],
                                    endpoints=np.asarray(endpoints)))
        for child in element:
            visit(child, matrix)

    visit(root, np.eye(3))
    return records


def axis_coordinates(points, specification):
    """Map page coordinates to the linear data axes using independently read ticks."""
    left, top, right, bottom = specification["frame"]
    x_min, x_max, y_min, y_max = specification["limits"]
    points = np.asarray(points)
    x = x_min + (points[:, 0] - left) / (right - left) * (x_max - x_min)
    y = y_min + (bottom - points[:, 1]) / (bottom - top) * (y_max - y_min)
    return x, y


def marker_groups(records, figure, specification):
    """Recover visible mean-line marker centers, excluding legend samples.

    Circles and diamonds use symmetric bounds; triangle centers use the
    centroid of their three vertices. Points clipped out of the published
    axes cannot be recovered and are not invented or extrapolated.
    """
    if figure == 3:
        return []
    colors = [("30.979919%", "theta_0", "Mean-line, 0 deg", "o"),
              ("70.979309%", "theta_20", "Mean-line, 20 deg", "^"),
              ("98.429871%", "theta_40", "Mean-line, 40 deg", "d")]
    if figure >= 6:
        colors = [("98.429871%", "meanline", "Mean-line", "o")]
    groups = []
    for color, key, label, marker in colors:
        points = []
        for record in records:
            style = record["style"]
            if color not in style:
                continue
            left, top, right, bottom = record["bounds"]
            if not (.8 < right - left < 5 and .8 < bottom - top < 5):
                continue
            if marker == "o" and "stroke:none" not in style:
                continue
            center = np.array([(left + right) / 2, (top + bottom) / 2])
            if marker == "^":
                center = record["endpoints"][:3].mean(axis=0)
            x0, y0, x1, y1 = specification["frame"]
            if not (x0 - .05 <= center[0] <= x1 + .05 and y0 <= center[1] <= y1):
                continue
            legend = specification["legend"]
            if legend and legend[0] <= center[0] <= legend[2] and legend[1] <= center[1] <= legend[3]:
                continue
            points.append(center)
        if not points:
            raise ValueError(f"No mean-line points found for Figure {figure}, {key}")
        x, y = axis_coordinates(points, specification)
        order = np.argsort(x)
        groups.append(dict(key=key, label=label, kind="published_meanline", marker=marker,
                           x=x[order].tolist(), y=y[order].tolist()))
    return groups


def extract(pdf, work_directory):
    """Extract all seven component figures, preserving source hash and calibration.

    This configuration targets the DTU open-access copy with its cover page.
    Refuse a different file rather than silently applying wrong page/path
    indices. The output is publication-vector digitization, not raw solver or
    experimental data, and carries the paper's CC BY 4.0 attribution.
    """
    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if digest != PDF_HASH:
        raise ValueError("PDF differs from the reviewed DTU copy; review page and axis calibration first.")
    program = shutil.which("pdftocairo")
    if program is None:
        raise RuntimeError("pdftocairo is required only to re-extract the PDF references.")
    work_directory.mkdir(parents=True, exist_ok=True)
    pages = {}
    for page in sorted({item["page"] for item in FIGURES.values()}):
        svg = work_directory / f"page_{page:02d}.svg"
        subprocess.run([program, "-f", str(page), "-l", str(page), "-svg", str(pdf), str(svg)],
                       check=True, capture_output=True)
        pages[page] = svg_paths(svg)
    figures = {}
    for number, specification in FIGURES.items():
        records = pages[specification["page"]]
        series = []
        for index, key, label, compared in specification["curves"]:
            record = records[index]
            if "fill:none;stroke-width:" not in record["style"]:
                raise ValueError(f"Unexpected source path for Figure {number}: {index}")
            x, y = axis_coordinates(record["endpoints"], specification)
            if not np.all(np.diff(x) > 0):
                raise ValueError(f"Source curve is not ordered for Figure {number}, {key}")
            series.append(dict(key=key, label=label, kind="published_correlation", compared=compared,
                               svg_path_index=index, x=x.tolist(), y=y.tolist()))
        series.extend(marker_groups(records, number, specification))
        figures[str(number)] = dict(pdf_page=specification["page"], axis_frame_pdf_points=specification["frame"],
                                    axis_limits=specification["limits"], series=series)
    return dict(source=dict(title="Advancing non-dimensional models for radial turbomachinery applied to pumped thermal energy storage",
                            authors="Parisi, Agromayor, La Seta, Nielsen and Haglind (2025)",
                            doi="10.1016/j.enconman.2025.120085", url=PDF_URL, sha256=digest,
                            license="CC BY 4.0", license_url="https://creativecommons.org/licenses/by/4.0/"),
                method="PDF to SVG with Poppler; affine transforms and linear-axis calibration; no fitting to thermopt.",
                limitations=["Published curves and marker positions, not original mean-line or experimental datasets.",
                             "Coordinates have PDF rounding and clipping/interpolation error.",
                             "Only visible marker centers are retained; hidden/clipped points are not reconstructed."],
                figures=figures)


def main():
    """Extract reviewed paper references only with an explicit source PDF and overwrite flag."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--work-directory", type=Path, default=ROOT / ".test_artifacts/paper_vectors")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        parser.error("Reference exists; review extraction changes before using --overwrite.")
    data = extract(args.pdf, args.work_directory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Extracted {len(data['figures'])} figures to {args.output}")
    for number, figure in data["figures"].items():
        print(f"Figure {number}: " + ", ".join(f"{s['key']}={len(s['x'])}" for s in figure["series"]))


if __name__ == "__main__":
    main()
