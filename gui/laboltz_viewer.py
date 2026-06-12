from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass, field
from pathlib import Path
from tkinter import (
    BOTH,
    BOTTOM,
    END,
    LEFT,
    RIGHT,
    TOP,
    X,
    Y,
    Canvas,
    filedialog,
    messagebox,
    StringVar,
    Tk,
)
from tkinter import ttk


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs"


@dataclass
class VtkData:
    path: Path
    dimensions: tuple[int, int, int]
    scalars: dict[str, list[float]] = field(default_factory=dict)
    vectors: dict[str, list[tuple[float, float, float]]] = field(default_factory=dict)

    @property
    def cell_count(self) -> int:
        nx, ny, nz = self.dimensions
        return nx * ny * nz


@dataclass
class CaseData:
    output_dir: Path
    case_name: str = "case"
    history: list[dict[str, float]] = field(default_factory=list)
    profile: list[dict[str, float]] = field(default_factory=list)
    vtk: VtkData | None = None


def read_float_csv(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        return []

    with path.open(newline="") as handle:
        rows: list[dict[str, float]] = []
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append({key: float(value) for key, value in row.items() if key is not None})
        return rows


def parse_legacy_vtk(path: Path) -> VtkData:
    lines = path.read_text(encoding="utf-8").splitlines()
    dimensions: tuple[int, int, int] | None = None
    point_count: int | None = None
    scalars: dict[str, list[float]] = {}
    vectors: dict[str, list[tuple[float, float, float]]] = {}

    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
            continue

        if parts[0] == "DIMENSIONS":
            dimensions = (int(parts[1]), int(parts[2]), int(parts[3]))
            index += 1
            continue

        if parts[0] == "POINT_DATA":
            point_count = int(parts[1])
            index += 1
            continue

        if parts[0] == "SCALARS":
            if point_count is None:
                raise ValueError("VTK SCALARS section appeared before POINT_DATA")
            name = parts[1]
            index += 1
            if index < len(lines) and lines[index].startswith("LOOKUP_TABLE"):
                index += 1
            values = [float(lines[index + offset].strip()) for offset in range(point_count)]
            scalars[name] = values
            index += point_count
            continue

        if parts[0] == "VECTORS":
            if point_count is None:
                raise ValueError("VTK VECTORS section appeared before POINT_DATA")
            name = parts[1]
            values: list[tuple[float, float, float]] = []
            index += 1
            for offset in range(point_count):
                vx, vy, vz = (float(value) for value in lines[index + offset].split()[:3])
                values.append((vx, vy, vz))
            vectors[name] = values
            index += point_count
            continue

        index += 1

    if dimensions is None:
        raise ValueError("Could not find VTK DIMENSIONS")
    vtk = VtkData(path=path, dimensions=dimensions, scalars=scalars, vectors=vectors)
    if point_count is not None and vtk.cell_count != point_count:
        raise ValueError("VTK POINT_DATA count does not match DIMENSIONS")
    return vtk


def case_label(case_name: str) -> str:
    return case_name.replace("_", " ").title()


def discover_cases(output_dir: Path) -> list[str]:
    output_dir = output_dir.resolve()
    cases: set[str] = set()
    for path in output_dir.glob("*_history.csv"):
        cases.add(path.name.removesuffix("_history.csv"))
    for path in output_dir.glob("*_profile.csv"):
        cases.add(path.name.removesuffix("_profile.csv"))
    for path in output_dir.glob("*_final.vtk"):
        cases.add(path.name.removesuffix("_final.vtk"))
    return sorted(cases)


def latest_case_file_time(output_dir: Path, case_name: str | None = None) -> float:
    newest = -1.0
    prefix = case_name if case_name else "*"
    for pattern in [f"{prefix}_history.csv", f"{prefix}_profile.csv", f"{prefix}_final.vtk", f"{prefix}_step_*.vtk"]:
        for path in output_dir.glob(pattern):
            if path.is_file():
                newest = max(newest, path.stat().st_mtime)
    return newest


def resolve_output_dir(output_dir: Path, case_name: str | None = None) -> Path:
    output_dir = output_dir.resolve()
    candidates: list[tuple[float, Path]] = []

    direct_time = latest_case_file_time(output_dir, case_name)
    if direct_time >= 0.0:
        candidates.append((direct_time, output_dir))

    for child in output_dir.iterdir() if output_dir.exists() else []:
        if child.is_dir():
            modified = latest_case_file_time(child, case_name)
            if modified >= 0.0:
                candidates.append((modified, child))

    if not candidates:
        return output_dir

    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def latest_case(output_dir: Path, cases: list[str]) -> str | None:
    newest_case_name: str | None = None
    newest_time = -1.0
    for name in cases:
        paths = [
            output_dir / f"{name}_history.csv",
            output_dir / f"{name}_profile.csv",
            output_dir / f"{name}_final.vtk",
        ]
        paths.extend(output_dir.glob(f"{name}_step_*.vtk"))
        for path in paths:
            if path.exists():
                modified = path.stat().st_mtime
                if modified > newest_time:
                    newest_time = modified
                    newest_case_name = name
    return newest_case_name


def load_case(output_dir: Path, case_name: str | None = None) -> CaseData:
    output_dir = resolve_output_dir(output_dir, case_name)
    cases = discover_cases(output_dir)
    if case_name is None:
        case_name = latest_case(output_dir, cases) or (cases[0] if cases else "case")

    final_vtk = output_dir / f"{case_name}_final.vtk"
    if not final_vtk.exists():
        vtk_files = sorted(output_dir.glob(f"{case_name}_*.vtk"))
        final_vtk = vtk_files[-1] if vtk_files else final_vtk

    return CaseData(
        output_dir=output_dir,
        case_name=case_name,
        history=read_float_csv(output_dir / f"{case_name}_history.csv"),
        profile=read_float_csv(output_dir / f"{case_name}_profile.csv"),
        vtk=parse_legacy_vtk(final_vtk) if final_vtk.exists() else None,
    )


def finite_range(values: list[float]) -> tuple[float, float]:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return 0.0, 1.0
    minimum = min(finite)
    maximum = max(finite)
    if minimum == maximum:
        padding = abs(minimum) * 0.1 if minimum else 1.0
        return minimum - padding, maximum + padding
    padding = 0.05 * (maximum - minimum)
    return minimum - padding, maximum + padding


def row_columns(rows: list[dict[str, float]]) -> list[str]:
    if not rows:
        return []
    return list(rows[0].keys())


def preferred_x_column(rows: list[dict[str, float]]) -> str:
    columns = row_columns(rows)
    for name in ["step", "y_from_lower_wall", "coordinate", "x", "y"]:
        if name in columns:
            return name
    return columns[0] if columns else ""


def plottable_columns(rows: list[dict[str, float]], x_column: str) -> list[str]:
    return [name for name in row_columns(rows) if name != x_column]


def color_ramp(value: float, minimum: float, maximum: float) -> str:
    if maximum <= minimum:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - minimum) / (maximum - minimum)))

    stops = [
        (0.00, (35, 51, 94)),
        (0.35, (42, 157, 143)),
        (0.70, (239, 201, 76)),
        (1.00, (204, 76, 57)),
    ]
    for left, right in zip(stops, stops[1:]):
        if left[0] <= t <= right[0]:
            local = (t - left[0]) / (right[0] - left[0])
            red = round(left[1][0] + local * (right[1][0] - left[1][0]))
            green = round(left[1][1] + local * (right[1][1] - left[1][1]))
            blue = round(left[1][2] + local * (right[1][2] - left[1][2]))
            return f"#{red:02x}{green:02x}{blue:02x}"
    red, green, blue = stops[-1][1]
    return f"#{red:02x}{green:02x}{blue:02x}"


class PlotCanvas(Canvas):
    def __init__(self, parent, **kwargs):
        super().__init__(parent, background="#f8fafc", highlightthickness=0, **kwargs)

    def plot_lines(
        self,
        x_values: list[float],
        series: list[tuple[str, list[float], str]],
        title: str,
        x_label: str,
        y_label: str,
    ) -> None:
        self.delete("all")
        width = max(self.winfo_width(), 640)
        height = max(self.winfo_height(), 360)
        left_margin = 72
        right_margin = 28
        top_margin = 48
        bottom_margin = 56
        plot_width = width - left_margin - right_margin
        plot_height = height - top_margin - bottom_margin

        y_values = [value for _, values, _ in series for value in values]
        x_min, x_max = finite_range(x_values)
        y_min, y_max = finite_range(y_values)

        def project_x(value: float) -> float:
            return left_margin + (value - x_min) / (x_max - x_min) * plot_width

        def project_y(value: float) -> float:
            return top_margin + (y_max - value) / (y_max - y_min) * plot_height

        self.create_text(left_margin, 22, text=title, anchor="w", fill="#172033", font=("Segoe UI", 13, "bold"))
        self.create_rectangle(left_margin, top_margin, left_margin + plot_width, top_margin + plot_height, outline="#c8d0dc")

        for tick in range(6):
            fraction = tick / 5
            x = left_margin + fraction * plot_width
            y = top_margin + fraction * plot_height
            x_value = x_min + fraction * (x_max - x_min)
            y_value = y_max - fraction * (y_max - y_min)
            self.create_line(x, top_margin, x, top_margin + plot_height, fill="#edf1f5")
            self.create_line(left_margin, y, left_margin + plot_width, y, fill="#edf1f5")
            self.create_text(x, top_margin + plot_height + 18, text=f"{x_value:.3g}", fill="#42526b", font=("Segoe UI", 9))
            self.create_text(left_margin - 10, y, text=f"{y_value:.3g}", anchor="e", fill="#42526b", font=("Segoe UI", 9))

        self.create_text(left_margin + plot_width / 2, height - 18, text=x_label, fill="#172033", font=("Segoe UI", 10))
        self.create_text(18, top_margin + plot_height / 2, text=y_label, fill="#172033", font=("Segoe UI", 10), angle=90)

        legend_x = left_margin + 12
        legend_y = top_margin + 12
        for label, values, color in series:
            points: list[float] = []
            for x_value, y_value in zip(x_values, values):
                points.extend([project_x(x_value), project_y(y_value)])
            if len(points) >= 4:
                self.create_line(*points, fill=color, width=2)
            for x_value, y_value in zip(x_values, values):
                self.create_oval(project_x(x_value) - 3, project_y(y_value) - 3, project_x(x_value) + 3, project_y(y_value) + 3, fill=color, outline="")
            self.create_rectangle(legend_x, legend_y - 5, legend_x + 14, legend_y + 5, fill=color, outline="")
            self.create_text(legend_x + 20, legend_y, text=label, anchor="w", fill="#172033", font=("Segoe UI", 9))
            legend_y += 18

    def plot_heatmap(
        self,
        matrix: list[list[float]],
        title: str,
        x_label: str,
        y_label: str,
    ) -> None:
        self.delete("all")
        width = max(self.winfo_width(), 640)
        height = max(self.winfo_height(), 360)
        if not matrix or not matrix[0]:
            self.create_text(width / 2, height / 2, text="No field data loaded", fill="#42526b", font=("Segoe UI", 12))
            return

        rows = len(matrix)
        cols = len(matrix[0])
        values = [value for row in matrix for value in row]
        minimum, maximum = finite_range(values)
        left_margin = 72
        right_margin = 96
        top_margin = 48
        bottom_margin = 56
        plot_width = width - left_margin - right_margin
        plot_height = height - top_margin - bottom_margin
        cell_width = plot_width / cols
        cell_height = plot_height / rows

        self.create_text(left_margin, 22, text=title, anchor="w", fill="#172033", font=("Segoe UI", 13, "bold"))
        for row_index, row in enumerate(matrix):
            for col_index, value in enumerate(row):
                x0 = left_margin + col_index * cell_width
                y0 = top_margin + (rows - 1 - row_index) * cell_height
                self.create_rectangle(
                    x0,
                    y0,
                    x0 + cell_width + 1,
                    y0 + cell_height + 1,
                    fill=color_ramp(value, minimum, maximum),
                    outline="",
                )

        self.create_rectangle(left_margin, top_margin, left_margin + plot_width, top_margin + plot_height, outline="#c8d0dc")
        self.create_text(left_margin + plot_width / 2, height - 18, text=x_label, fill="#172033", font=("Segoe UI", 10))
        self.create_text(18, top_margin + plot_height / 2, text=y_label, fill="#172033", font=("Segoe UI", 10), angle=90)

        legend_x = left_margin + plot_width + 28
        legend_y = top_margin
        legend_height = plot_height
        for step in range(80):
            fraction0 = step / 80
            fraction1 = (step + 1) / 80
            y0 = legend_y + (1.0 - fraction1) * legend_height
            y1 = legend_y + (1.0 - fraction0) * legend_height
            value = minimum + fraction0 * (maximum - minimum)
            self.create_rectangle(legend_x, y0, legend_x + 18, y1, fill=color_ramp(value, minimum, maximum), outline="")
        self.create_rectangle(legend_x, legend_y, legend_x + 18, legend_y + legend_height, outline="#c8d0dc")
        self.create_text(legend_x + 28, legend_y, text=f"{maximum:.3g}", anchor="w", fill="#42526b", font=("Segoe UI", 9))
        self.create_text(legend_x + 28, legend_y + legend_height, text=f"{minimum:.3g}", anchor="w", fill="#42526b", font=("Segoe UI", 9))


def vtk_field_names(vtk: VtkData | None) -> list[str]:
    if vtk is None:
        return []
    names = list(vtk.scalars)
    if "velocity" in vtk.vectors:
        names.extend(["velocity_x", "velocity_y", "velocity_z", "velocity_magnitude"])
    return names


def vtk_field_values(vtk: VtkData, field_name: str) -> list[float]:
    if field_name in vtk.scalars:
        return vtk.scalars[field_name]
    if "velocity" not in vtk.vectors:
        return []
    velocity = vtk.vectors["velocity"]
    if field_name == "velocity_x":
        return [value[0] for value in velocity]
    if field_name == "velocity_y":
        return [value[1] for value in velocity]
    if field_name == "velocity_z":
        return [value[2] for value in velocity]
    if field_name == "velocity_magnitude":
        return [math.sqrt(value[0] ** 2 + value[1] ** 2 + value[2] ** 2) for value in velocity]
    return []


def vtk_xy_slice(vtk: VtkData, field_name: str, z_index: int) -> list[list[float]]:
    nx, ny, nz = vtk.dimensions
    z_index = max(0, min(nz - 1, z_index))
    values = vtk_field_values(vtk, field_name)
    if not values:
        return []
    matrix: list[list[float]] = []
    for y in range(ny):
        row: list[float] = []
        for x in range(nx):
            row.append(values[(z_index * ny + y) * nx + x])
        matrix.append(row)
    return matrix


class LaboltzViewer:
    def __init__(self, root: Tk, initial_dir: Path) -> None:
        self.root = root
        self.root.title("LaBoltz Research Viewer")
        self.root.geometry("1180x760")
        self.root.minsize(960, 640)
        self.output_dir = initial_dir
        self.case = CaseData(output_dir=initial_dir)
        self.available_cases: list[str] = []
        self.case_var = StringVar(value="")
        self.metric_var = StringVar(value="relative_l2_error")
        self.field_var = StringVar(value="velocity_x")
        self.slice_var = StringVar(value="0")

        self._configure_style()
        self._build_layout()
        self.load_output_dir(initial_dir)

    def _configure_style(self) -> None:
        style = ttk.Style()
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("TFrame", background="#eef2f6")
        style.configure("Panel.TFrame", background="#ffffff")
        style.configure("TLabel", background="#eef2f6", foreground="#172033", font=("Segoe UI", 10))
        style.configure("Panel.TLabel", background="#ffffff", foreground="#172033", font=("Segoe UI", 10))
        style.configure("Title.TLabel", background="#eef2f6", foreground="#172033", font=("Segoe UI", 16, "bold"))
        style.configure("TButton", font=("Segoe UI", 10), padding=(10, 6))
        style.configure("TNotebook", background="#eef2f6")
        style.configure("TNotebook.Tab", font=("Segoe UI", 10), padding=(12, 6))

    def _build_layout(self) -> None:
        main = ttk.Frame(self.root, padding=12)
        main.pack(fill=BOTH, expand=True)

        top_bar = ttk.Frame(main)
        top_bar.pack(side=TOP, fill=X, pady=(0, 10))
        ttk.Label(top_bar, text="LaBoltz Research Viewer", style="Title.TLabel").pack(side=LEFT)
        ttk.Button(top_bar, text="Open Output Folder", command=self.choose_output_dir).pack(side=RIGHT, padx=(6, 0))
        ttk.Button(top_bar, text="Export HTML Report", command=self.export_report).pack(side=RIGHT, padx=(6, 0))
        ttk.Button(top_bar, text="Refresh", command=lambda: self.load_output_dir(self.output_dir)).pack(side=RIGHT, padx=(6, 0))
        self.case_box = ttk.Combobox(top_bar, textvariable=self.case_var, state="readonly", width=18)
        self.case_box.pack(side=RIGHT, padx=(6, 12))
        self.case_box.bind("<<ComboboxSelected>>", lambda _event: self.load_output_dir(self.output_dir, self.case_var.get()))
        ttk.Label(top_bar, text="Case").pack(side=RIGHT)

        self.status = ttk.Label(main, text="", anchor="w")
        self.status.pack(side=BOTTOM, fill=X, pady=(8, 0))

        self.notebook = ttk.Notebook(main)
        self.notebook.pack(fill=BOTH, expand=True)

        self.summary_frame = ttk.Frame(self.notebook, style="Panel.TFrame", padding=14)
        self.convergence_frame = ttk.Frame(self.notebook, style="Panel.TFrame", padding=10)
        self.profile_frame = ttk.Frame(self.notebook, style="Panel.TFrame", padding=10)
        self.field_frame = ttk.Frame(self.notebook, style="Panel.TFrame", padding=10)
        self.table_frame = ttk.Frame(self.notebook, style="Panel.TFrame", padding=10)

        self.notebook.add(self.summary_frame, text="Summary")
        self.notebook.add(self.convergence_frame, text="Convergence")
        self.notebook.add(self.profile_frame, text="Profile")
        self.notebook.add(self.field_frame, text="Field Slice")
        self.notebook.add(self.table_frame, text="History Table")

        self.summary_text = ttk.Label(self.summary_frame, text="", style="Panel.TLabel", justify=LEFT)
        self.summary_text.pack(anchor="nw")

        convergence_controls = ttk.Frame(self.convergence_frame, style="Panel.TFrame")
        convergence_controls.pack(side=TOP, fill=X, pady=(0, 8))
        ttk.Label(convergence_controls, text="Metric", style="Panel.TLabel").pack(side=LEFT, padx=(0, 6))
        self.metric_box = ttk.Combobox(
            convergence_controls,
            textvariable=self.metric_var,
            values=["max_ux", "average_ux", "relative_l2_error", "relative_max_change", "fluid_mass"],
            state="readonly",
            width=24,
        )
        self.metric_box.pack(side=LEFT)
        self.metric_box.bind("<<ComboboxSelected>>", lambda _event: self.refresh_plots())
        self.convergence_canvas = PlotCanvas(self.convergence_frame)
        self.convergence_canvas.pack(fill=BOTH, expand=True)
        self.convergence_canvas.bind("<Configure>", lambda _event: self.refresh_convergence_plot())

        self.profile_canvas = PlotCanvas(self.profile_frame)
        self.profile_canvas.pack(fill=BOTH, expand=True)
        self.profile_canvas.bind("<Configure>", lambda _event: self.refresh_profile_plot())

        field_controls = ttk.Frame(self.field_frame, style="Panel.TFrame")
        field_controls.pack(side=TOP, fill=X, pady=(0, 8))
        ttk.Label(field_controls, text="Field", style="Panel.TLabel").pack(side=LEFT, padx=(0, 6))
        self.field_box = ttk.Combobox(field_controls, textvariable=self.field_var, state="readonly", width=24)
        self.field_box.pack(side=LEFT, padx=(0, 14))
        self.field_box.bind("<<ComboboxSelected>>", lambda _event: self.refresh_field_plot())
        ttk.Label(field_controls, text="Z slice", style="Panel.TLabel").pack(side=LEFT, padx=(0, 6))
        self.slice_box = ttk.Spinbox(field_controls, from_=0, to=0, textvariable=self.slice_var, width=6, command=self.refresh_field_plot)
        self.slice_box.pack(side=LEFT)
        self.slice_box.bind("<Return>", lambda _event: self.refresh_field_plot())
        self.field_canvas = PlotCanvas(self.field_frame)
        self.field_canvas.pack(fill=BOTH, expand=True)
        self.field_canvas.bind("<Configure>", lambda _event: self.refresh_field_plot())

        self.history_table = ttk.Treeview(self.table_frame, show="headings", height=18)
        table_scroll = ttk.Scrollbar(self.table_frame, orient="vertical", command=self.history_table.yview)
        self.history_table.configure(yscrollcommand=table_scroll.set)
        self.history_table.pack(side=LEFT, fill=BOTH, expand=True)
        table_scroll.pack(side=RIGHT, fill=Y)

    def choose_output_dir(self) -> None:
        directory = filedialog.askdirectory(initialdir=self.output_dir)
        if directory:
            self.load_output_dir(Path(directory))

    def load_output_dir(self, output_dir: Path, case_name: str | None = None) -> None:
        try:
            self.case = load_case(output_dir, case_name)
            self.output_dir = self.case.output_dir
            self.available_cases = discover_cases(self.output_dir)
            self.case_var.set(self.case.case_name)
            self.case_box.configure(values=self.available_cases)
            self.status.configure(text=f"Loaded {case_label(self.case.case_name)} from {self.output_dir}")
            self.refresh_all()
        except Exception as error:
            messagebox.showerror("Could not load output folder", str(error))
            self.status.configure(text=f"Load failed: {error}")

    def refresh_all(self) -> None:
        self.refresh_summary()
        self.refresh_metric_controls()
        self.refresh_field_controls()
        self.refresh_plots()
        self.refresh_history_table()

    def refresh_plots(self) -> None:
        self.refresh_convergence_plot()
        self.refresh_profile_plot()
        self.refresh_field_plot()

    def refresh_summary(self) -> None:
        history = self.case.history
        profile = self.case.profile
        vtk = self.case.vtk
        lines = [f"Output folder: {self.case.output_dir}"]
        lines.append(f"Case: {case_label(self.case.case_name)}")

        if history:
            final = history[-1]
            lines.extend(["", f"{case_label(self.case.case_name)} convergence"])
            for key, value in final.items():
                if key == "step":
                    lines.append(f"  final step: {value:.0f}")
                else:
                    lines.append(f"  {key}: {value:.6e}")
        else:
            lines.extend(["", f"No {self.case.case_name}_history.csv found."])

        if profile:
            lines.extend(["", f"Profile rows: {len(profile)}"])
            if "ux_mean" in profile[0] and "ux_analytical" in profile[0]:
                max_error = max(abs(row["ux_mean"] - row["ux_analytical"]) for row in profile)
                lines.append(f"maximum profile error: {max_error:.6e}")

        if vtk is not None:
            lines.extend(
                [
                    "",
                    f"VTK file: {vtk.path.name}",
                    f"dimensions: {vtk.dimensions[0]} x {vtk.dimensions[1]} x {vtk.dimensions[2]}",
                    f"scalar fields: {', '.join(vtk.scalars) if vtk.scalars else 'none'}",
                    f"vector fields: {', '.join(vtk.vectors) if vtk.vectors else 'none'}",
                ]
            )
        else:
            lines.extend(["", "No VTK file found."])

        self.summary_text.configure(text="\n".join(lines))

    def refresh_metric_controls(self) -> None:
        metrics = plottable_columns(self.case.history, "step")
        if not metrics:
            metrics = ["max_ux", "average_ux", "relative_l2_error", "relative_max_change", "fluid_mass"]
        self.metric_box.configure(values=metrics)
        if self.metric_var.get() not in metrics:
            self.metric_var.set(metrics[0])

    def refresh_convergence_plot(self) -> None:
        history = self.case.history
        if not history:
            self.convergence_canvas.delete("all")
            self.convergence_canvas.create_text(320, 180, text="No convergence history loaded", fill="#42526b", font=("Segoe UI", 12))
            return

        metric = self.metric_var.get()
        x_values = [row["step"] for row in history]
        y_values = [row.get(metric, 0.0) for row in history]
        self.convergence_canvas.plot_lines(
            x_values,
            [(metric, y_values, "#2a7fbb")],
            f"Convergence: {metric}",
            "step",
            metric,
        )

    def refresh_profile_plot(self) -> None:
        profile = self.case.profile
        if not profile:
            self.profile_canvas.delete("all")
            self.profile_canvas.create_text(320, 180, text="No profile CSV loaded", fill="#42526b", font=("Segoe UI", 12))
            return

        if "y_from_lower_wall" in profile[0]:
            x_column = "y_from_lower_wall"
        else:
            x_column = preferred_x_column(profile)
        x_values = [row[x_column] for row in profile]

        colors = ["#2a7fbb", "#c94f3d", "#2a9d8f", "#8b5fbf", "#e9a03f"]
        if "ux_mean" in profile[0] and "ux_analytical" in profile[0]:
            series = [
                ("LBM ux", [row["ux_mean"] for row in profile], colors[0]),
                ("analytical ux", [row["ux_analytical"] for row in profile], colors[1]),
            ]
        else:
            series = [
                (column, [row[column] for row in profile], colors[index % len(colors)])
                for index, column in enumerate(plottable_columns(profile, x_column))
            ]

        self.profile_canvas.plot_lines(
            x_values,
            series,
            f"{case_label(self.case.case_name)} Profile",
            x_column,
            "velocity",
        )

    def refresh_field_controls(self) -> None:
        names = vtk_field_names(self.case.vtk)
        self.field_box.configure(values=names)
        if names and self.field_var.get() not in names:
            self.field_var.set(names[0])
        if self.case.vtk is not None:
            max_slice = self.case.vtk.dimensions[2] - 1
            self.slice_box.configure(to=max_slice)
            if int_or_zero(self.slice_var.get()) > max_slice:
                self.slice_var.set(str(max_slice))

    def refresh_field_plot(self) -> None:
        vtk = self.case.vtk
        if vtk is None:
            self.field_canvas.delete("all")
            self.field_canvas.create_text(320, 180, text="No VTK field loaded", fill="#42526b", font=("Segoe UI", 12))
            return
        field_name = self.field_var.get()
        z_index = int_or_zero(self.slice_var.get())
        matrix = vtk_xy_slice(vtk, field_name, z_index)
        self.field_canvas.plot_heatmap(matrix, f"{field_name} at z={z_index}", "x", "y")

    def refresh_history_table(self) -> None:
        for column in self.history_table["columns"]:
            self.history_table.heading(column, text="")
        self.history_table.delete(*self.history_table.get_children())
        history = self.case.history
        if not history:
            self.history_table.configure(columns=[])
            return

        columns = list(history[0].keys())
        self.history_table.configure(columns=columns)
        for column in columns:
            self.history_table.heading(column, text=column)
            self.history_table.column(column, width=140, anchor="e")
        for row in history:
            self.history_table.insert("", END, values=[f"{row[column]:.8g}" for column in columns])

    def export_report(self) -> None:
        try:
            report_path = self.case.output_dir / f"{self.case.case_name}_report.html"
            report_path.write_text(render_html_report(self.case), encoding="utf-8")
            self.status.configure(text=f"Wrote {report_path}")
            messagebox.showinfo("Report exported", f"Wrote {report_path}")
        except Exception as error:
            messagebox.showerror("Could not export report", str(error))


def int_or_zero(value: str) -> int:
    try:
        return int(value)
    except ValueError:
        return 0


def render_svg_line_plot(
    x_values: list[float],
    series: list[tuple[str, list[float], str]],
    title: str,
    x_label: str,
    y_label: str,
) -> str:
    width = 860
    height = 420
    left = 76
    right = 32
    top = 54
    bottom = 58
    plot_width = width - left - right
    plot_height = height - top - bottom
    y_values = [value for _, values, _ in series for value in values]
    x_min, x_max = finite_range(x_values)
    y_min, y_max = finite_range(y_values)

    def px(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def py(value: float) -> float:
        return top + (y_max - value) / (y_max - y_min) * plot_height

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title)}">',
        '<rect width="100%" height="100%" fill="#f8fafc"/>',
        f'<text x="{left}" y="30" font-family="Segoe UI, Arial" font-size="18" font-weight="700" fill="#172033">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="#ffffff" stroke="#c8d0dc"/>',
    ]
    for tick in range(6):
        fraction = tick / 5
        x = left + fraction * plot_width
        y = top + fraction * plot_height
        x_value = x_min + fraction * (x_max - x_min)
        y_value = y_max - fraction * (y_max - y_min)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_height}" stroke="#edf1f5"/>')
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_width}" y2="{y:.2f}" stroke="#edf1f5"/>')
        parts.append(f'<text x="{x:.2f}" y="{top + plot_height + 24}" text-anchor="middle" font-size="12" fill="#42526b">{x_value:.3g}</text>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-size="12" fill="#42526b">{y_value:.3g}</text>')
    parts.append(f'<text x="{left + plot_width / 2}" y="{height - 16}" text-anchor="middle" font-size="13" fill="#172033">{html.escape(x_label)}</text>')
    parts.append(f'<text x="18" y="{top + plot_height / 2}" transform="rotate(-90 18 {top + plot_height / 2})" text-anchor="middle" font-size="13" fill="#172033">{html.escape(y_label)}</text>')

    legend_y = top + 18
    for label, values, color in series:
        points = " ".join(f"{px(x):.2f},{py(y):.2f}" for x, y in zip(x_values, values))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        parts.append(f'<rect x="{left + 14}" y="{legend_y - 9}" width="16" height="10" fill="{color}"/>')
        parts.append(f'<text x="{left + 38}" y="{legend_y}" font-size="12" fill="#172033">{html.escape(label)}</text>')
        legend_y += 20
    parts.append("</svg>")
    return "\n".join(parts)


def render_html_report(case: CaseData) -> str:
    history = case.history
    profile = case.profile
    final = history[-1] if history else {}
    sections = [
        "<!doctype html>",
        "<html lang=\"en\">",
        "<head>",
        "<meta charset=\"utf-8\">",
        f"<title>LaBoltz {html.escape(case_label(case.case_name))} Report</title>",
        "<style>",
        "body{font-family:Segoe UI,Arial,sans-serif;margin:32px;color:#172033;background:#eef2f6}",
        "main{max-width:980px;margin:auto}",
        "section{background:white;padding:22px;margin:18px 0;border:1px solid #d8dee8}",
        "h1,h2{margin-top:0}",
        ".metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}",
        ".metric{background:#f8fafc;border:1px solid #d8dee8;padding:12px}",
        ".label{color:#5d6b82;font-size:12px}.value{font-size:18px;font-weight:700}",
        "</style>",
        "</head>",
        "<body><main>",
        f"<h1>LaBoltz {html.escape(case_label(case.case_name))} Report</h1>",
        f"<p>Output folder: {html.escape(str(case.output_dir))}</p>",
        "<section><h2>Final Metrics</h2><div class=\"metrics\">",
    ]

    metric_keys = list(final.keys()) if final else ["step"]
    for key in metric_keys:
        value = final.get(key)
        shown = f"{value:.6g}" if isinstance(value, float) else "n/a"
        sections.append(f"<div class=\"metric\"><div class=\"label\">{html.escape(key)}</div><div class=\"value\">{shown}</div></div>")
    sections.append("</div></section>")

    if history:
        x_values = [row["step"] for row in history]
        metrics = plottable_columns(history, "step")
        primary_metric = "relative_l2_error" if "relative_l2_error" in metrics else metrics[0]
        sections.append("<section>")
        sections.append(render_svg_line_plot(x_values, [(primary_metric, [row[primary_metric] for row in history], "#2a7fbb")], case_label(primary_metric), "step", primary_metric))
        sections.append("</section>")
        comparison_metrics = [name for name in ["max_ux", "average_ux", "max_speed", "average_speed"] if name in metrics]
        if len(comparison_metrics) >= 2:
            colors = ["#2a7fbb", "#c94f3d", "#2a9d8f", "#8b5fbf"]
            sections.append("<section>")
            sections.append(render_svg_line_plot(
                x_values,
                [(name, [row[name] for row in history], colors[index]) for index, name in enumerate(comparison_metrics[:4])],
                "Velocity Convergence",
                "step",
                "velocity",
            ))
            sections.append("</section>")

    if profile:
        x_column = "y_from_lower_wall" if "y_from_lower_wall" in profile[0] else preferred_x_column(profile)
        x_values = [row[x_column] for row in profile]
        if "ux_mean" in profile[0] and "ux_analytical" in profile[0]:
            series = [
                ("LBM ux", [row["ux_mean"] for row in profile], "#2a7fbb"),
                ("analytical ux", [row["ux_analytical"] for row in profile], "#c94f3d"),
            ]
        else:
            colors = ["#2a7fbb", "#c94f3d", "#2a9d8f", "#8b5fbf"]
            series = [
                (column, [row[column] for row in profile], colors[index % len(colors)])
                for index, column in enumerate(plottable_columns(profile, x_column)[:4])
            ]
        sections.append("<section>")
        sections.append(render_svg_line_plot(x_values, series, "Velocity Profile", x_column, "velocity"))
        sections.append("</section>")

    if case.vtk is not None:
        sections.append("<section>")
        sections.append(f"<h2>Field Data</h2><p>VTK file: {html.escape(case.vtk.path.name)}</p>")
        sections.append(f"<p>Dimensions: {case.vtk.dimensions[0]} x {case.vtk.dimensions[1]} x {case.vtk.dimensions[2]}</p>")
        sections.append("</section>")

    sections.extend(["</main></body></html>"])
    return "\n".join(sections)


def run_check(output_dir: Path, case_name: str | None = None) -> int:
    case = load_case(output_dir, case_name)
    print(f"Loaded output folder: {case.output_dir}")
    print(f"Case: {case.case_name}")
    print(f"History rows: {len(case.history)}")
    print(f"Profile rows: {len(case.profile)}")
    if case.history:
        final = case.history[-1]
        for key, value in final.items():
            if key == "step":
                print(f"Final step: {value:.0f}")
            else:
                print(f"Final {key}: {value:.6e}")
    if case.vtk is not None:
        print(f"VTK: {case.vtk.path.name}")
        print(f"Dimensions: {case.vtk.dimensions}")
        print(f"Fields: {', '.join(vtk_field_names(case.vtk))}")
    report = render_html_report(case)
    print(f"Report preview bytes: {len(report.encode('utf-8'))}")
    return 0


def export_report(output_dir: Path, case_name: str | None = None, report_path: Path | None = None) -> int:
    case = load_case(output_dir, case_name)
    if report_path is None:
        report_path = case.output_dir / f"{case.case_name}_report.html"
    report_path.write_text(render_html_report(case), encoding="utf-8")
    print(f"Wrote {report_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="LaBoltz research output viewer")
    parser.add_argument("output_dir", nargs="?", default=str(DEFAULT_OUTPUT_DIR), help="Directory containing solver outputs")
    parser.add_argument("--case", help="Case prefix to load, for example poiseuille or couette")
    parser.add_argument("--check", action="store_true", help="Load data and exit without opening the GUI")
    parser.add_argument("--export-report", action="store_true", help="Write an HTML report and exit without opening the GUI")
    parser.add_argument("--report-path", help="Optional output path for --export-report")
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    if args.check:
        return run_check(output_dir, args.case)
    if args.export_report:
        return export_report(output_dir, args.case, Path(args.report_path) if args.report_path else None)

    root = Tk()
    LaboltzViewer(root, output_dir)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
