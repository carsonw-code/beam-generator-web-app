
import io
import random
from dataclasses import dataclass, field
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from matplotlib.patches import Circle, Polygon


st.set_page_config(page_title="Random Beam Practice", layout="wide")

BEAM_YLIM = (-1.7, 4.1)
PIN_LABEL_Y = -1.10
TYPE_LABEL_Y = -1.35
POINT_LABEL_LEVELS = [3.05, 3.35, 3.65]
DIST_LABEL_LEVELS = [2.35, 2.62]
DIST_END_LABEL_Y = 1.82
MOMENT_LABEL_LEVELS = [-0.98, -1.25]

SUPPORT_TRI_HALF_WIDTH = 0.020
SUPPORT_TRI_HEIGHT = 0.120
SUPPORT_WHEEL_RADIUS = 0.0085
SUPPORT_WHEEL_GAP = 0.013
SUPPORT_BASE_EXTRA = 0.020
SUPPORT_BASE_HALF_WIDTH = 0.030
ROLLER_WHEEL_OFFSETS = (-0.018, 0.0, 0.018)

TEMPLATES = {
    10.0: [{"dist": (2.0, 6.0)}, {"dist": (4.0, 8.0)}, {"dist": (2.0, 8.0)}],
    12.0: [{"dist": (2.0, 6.0)}, {"dist": (4.0, 10.0)}, {"dist": (2.0, 8.0)}, {"dist": (4.0, 8.0)}],
}


@dataclass
class PointLoad:
    x: float
    P: float


@dataclass
class DistLoad:
    a: float
    b: float
    w1: float
    w2: float

    @property
    def length(self) -> float:
        return self.b - self.a

    @property
    def slope(self) -> float:
        return (self.w2 - self.w1) / self.length

    @property
    def intercept(self) -> float:
        return self.w1 - self.slope * self.a

    def resultant(self) -> Tuple[float, float]:
        W = 0.5 * (self.w1 + self.w2) * self.length
        xbar_from_a = self.length * (self.w1 + 2 * self.w2) / (3 * (self.w1 + self.w2))
        return W, self.a + xbar_from_a

    def intensity(self, x: np.ndarray) -> np.ndarray:
        q = np.zeros_like(x, dtype=float)
        mask = (x >= self.a) & (x <= self.b)
        q[mask] = self.slope * x[mask] + self.intercept
        return q

    def area_left_of(self, x: float) -> float:
        if x <= self.a:
            return 0.0
        xr = min(x, self.b)
        dx = xr - self.a
        return self.w1 * dx + 0.5 * self.slope * dx**2

    def moment_about_cut(self, x: float) -> float:
        if x <= self.a:
            return 0.0
        xr = min(x, self.b)
        dx = xr - self.a
        return 0.5 * self.w1 * dx**2 + (self.slope / 6.0) * dx**3


@dataclass
class FreeMoment:
    x: float
    M: float


@dataclass
class Section:
    start: float
    end: float
    v_coeffs: np.ndarray
    m_coeffs: np.ndarray
    v_text: str
    m_text: str
    range_text: str
    v_latex: str
    m_latex: str
    range_latex: str


@dataclass
class BeamProblem:
    L: float
    point_load: PointLoad
    dist_load: DistLoad
    free_moment: FreeMoment
    Ay: float
    By: float
    x: np.ndarray
    V: np.ndarray
    M: np.ndarray
    sections: List[Section] = field(default_factory=list)
    reaction_text: str = ""
    reaction_latex: str = ""
    dist_label: str = ""
    point_label: str = ""
    moment_label: str = ""


def fmt_num(value: float, digits: int = 2) -> str:
    if abs(value) < 1e-10:
        value = 0.0
    text = f"{value:.{digits}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def poly_text(coeffs: np.ndarray, start_power: int, digits: int = 2) -> str:
    parts = []
    degree = start_power
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if abs(coeff) < 1e-10:
            continue
        mag = fmt_num(abs(coeff), digits)
        if power == 0:
            term = mag
        elif power == 1:
            term = "x" if abs(abs(coeff) - 1.0) < 1e-10 else f"{mag}x"
        else:
            term = f"x^{power}" if abs(abs(coeff) - 1.0) < 1e-10 else f"{mag}x^{power}"

        if not parts:
            parts.append(term if coeff > 0 else f"-{term}")
        else:
            parts.append(f" {'+' if coeff > 0 else '-'} {term}")
    return "".join(parts) if parts else "0"


def poly_latex(coeffs: np.ndarray, start_power: int, digits: int = 2) -> str:
    parts = []
    degree = start_power
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if abs(coeff) < 1e-10:
            continue

        mag_value = abs(coeff)
        mag = fmt_num(mag_value, digits)

        if power == 0:
            term = mag
        elif power == 1:
            term = "x" if abs(mag_value - 1.0) < 1e-10 else f"{mag}x"
        else:
            term = rf"x^{{{power}}}" if abs(mag_value - 1.0) < 1e-10 else rf"{mag}x^{{{power}}}"

        if not parts:
            parts.append(term if coeff > 0 else f"-{term}")
        else:
            parts.append(f" + {term}" if coeff > 0 else f" - {term}")

    return "".join(parts) if parts else "0"


def label_lane(x_value: float, used_x: List[List[float]], levels: List[float], min_sep: float) -> float:
    for xs, level in zip(used_x, levels):
        if all(abs(x_value - old) >= min_sep for old in xs):
            xs.append(x_value)
            return level
    used_x[-1].append(x_value)
    return levels[-1]


def choose_station(L: float, forbidden: List[float], min_sep: float = 2.0):
    candidates = [float(x) for x in range(1, int(L)) if all(abs(float(x) - f) >= min_sep for f in forbidden)]
    return random.choice(candidates) if candidates else None


def section_start_values(beam: BeamProblem, s: float) -> Tuple[float, float, bool]:
    point = beam.point_load
    dist = beam.dist_load
    free = beam.free_moment

    V0 = beam.Ay - dist.area_left_of(s)
    if point.x <= s + 1e-10:
        V0 -= point.P

    M0 = beam.Ay * s - dist.moment_about_cut(s)
    if point.x <= s + 1e-10:
        M0 -= point.P * (s - point.x)
    if free.x <= s + 1e-10:
        M0 -= free.M

    active = dist.a <= (s + 1e-8) < dist.b
    return V0, M0, active


def build_section(beam: BeamProblem, start: float, end: float) -> Section:
    V0, M0, active = section_start_values(beam, start)

    if active:
        m = beam.dist_load.slope
        c = beam.dist_load.intercept
        a2 = -0.5 * m
        a1 = -c
        a0 = V0 + 0.5 * m * start**2 + c * start
    else:
        a2, a1, a0 = 0.0, 0.0, V0

    b3 = a2 / 3.0
    b2 = a1 / 2.0
    b1 = a0
    b0 = M0 - (b3 * start**3 + b2 * start**2 + b1 * start)

    v_coeffs = np.array([a2, a1, a0], dtype=float)
    m_coeffs = np.array([b3, b2, b1, b0], dtype=float)

    range_text = f"{fmt_num(start,1)} ≤ x ≤ {fmt_num(end,1)}"
    range_latex = rf"{fmt_num(start,1)} \le x \le {fmt_num(end,1)}"

    return Section(
        start=start,
        end=end,
        v_coeffs=v_coeffs,
        m_coeffs=m_coeffs,
        v_text=f"V(x) = {poly_text(v_coeffs, 2)}",
        m_text=f"M(x) = {poly_text(m_coeffs, 3)}",
        range_text=range_text,
        v_latex=rf"V(x) = {poly_latex(v_coeffs, 2)}",
        m_latex=rf"M(x) = {poly_latex(m_coeffs, 3)}",
        range_latex=range_latex,
    )


def precompute_display_data(beam: BeamProblem) -> None:
    points = sorted([0.0, beam.L, beam.point_load.x, beam.dist_load.a, beam.dist_load.b, beam.free_moment.x])
    beam.sections = [build_section(beam, a, b) for a, b in zip(points[:-1], points[1:]) if b - a > 1e-10]

    beam.reaction_text = f"Ay = {fmt_num(beam.Ay,2)} kN, By = {fmt_num(beam.By,2)} kN"
    beam.reaction_latex = rf"A_y = {fmt_num(beam.Ay,2)}\ \mathrm{{kN}}, \qquad B_y = {fmt_num(beam.By,2)}\ \mathrm{{kN}}"

    beam.point_label = rf"$P_1={fmt_num(beam.point_load.P,1)}\ \mathrm{{kN}}$"
    beam.moment_label = rf"$M_1={fmt_num(abs(beam.free_moment.M),1)}\ \mathrm{{kN\cdot m}}\ ({'CCW' if beam.free_moment.M >= 0 else 'CW'})$"

    d = beam.dist_load
    expr = fmt_num(d.w1, 1) if abs(d.w1 - d.w2) < 1e-10 else poly_text(np.array([d.slope, d.intercept]), 1)
    beam.dist_label = rf"$w(x)={expr}\ \mathrm{{kN/m}},\ {fmt_num(d.a,1)}\leq x \leq {fmt_num(d.b,1)}$"


def build_exact_curves(beam: BeamProblem, n: int = 2401) -> None:
    x = np.linspace(0.0, beam.L, n)
    V = np.zeros_like(x)
    M = np.zeros_like(x)

    for i, sec in enumerate(beam.sections):
        mask = (x >= sec.start) & (x < sec.end) if i < len(beam.sections) - 1 else (x >= sec.start) & (x <= sec.end)
        V[mask] = np.polyval(sec.v_coeffs, x[mask])
        M[mask] = np.polyval(sec.m_coeffs, x[mask])

    beam.x = x
    beam.V = V
    beam.M = M


def generate_problem() -> BeamProblem:
    while True:
        L = random.choice([10.0, 12.0])
        template = random.choice(TEMPLATES[L])

        dist = DistLoad(
            a=template["dist"][0],
            b=template["dist"][1],
            w1=random.choice([0.0, 5.0, 10.0, 15.0, 20.0]),
            w2=random.choice([0.0, 5.0, 10.0, 15.0, 20.0]),
        )
        if abs(dist.w1 - dist.w2) < 1e-10 and dist.w1 == 0.0:
            dist.w1 = dist.w2 = random.choice([5.0, 10.0, 15.0, 20.0])

        point_x = choose_station(L, [dist.a, dist.b], min_sep=2.0)
        moment_x = choose_station(L, [dist.a, dist.b, point_x] if point_x is not None else [dist.a, dist.b], min_sep=2.0)
        if point_x is None or moment_x is None:
            continue

        point = PointLoad(x=point_x, P=random.choice([10.0, 15.0, 20.0, 25.0, 30.0]))
        free = FreeMoment(x=moment_x, M=random.choice([-40.0, -20.0, 20.0, 40.0]))

        W, xbar = dist.resultant()
        By = (point.P * point.x + W * xbar - free.M) / L
        Ay = point.P + W - By
        if Ay <= 0 or By <= 0:
            continue

        beam = BeamProblem(
            L=L,
            point_load=point,
            dist_load=dist,
            free_moment=free,
            Ay=Ay,
            By=By,
            x=np.array([]),
            V=np.array([]),
            M=np.array([]),
        )
        precompute_display_data(beam)
        build_exact_curves(beam)

        if abs(beam.M[-1]) > 1.0:
            continue
        return beam


def draw_supports(ax, L: float) -> None:
    y_beam = (0.0 - BEAM_YLIM[0]) / (BEAM_YLIM[1] - BEAM_YLIM[0])
    trans = ax.transAxes

    for x_ax in (0.0, 1.0):
        triangle = Polygon(
            [(x_ax - SUPPORT_TRI_HALF_WIDTH, y_beam - SUPPORT_TRI_HEIGHT),
             (x_ax + SUPPORT_TRI_HALF_WIDTH, y_beam - SUPPORT_TRI_HEIGHT),
             (x_ax, y_beam)],
            closed=True, fill=False, linewidth=1.5, transform=trans, clip_on=False, zorder=6,
        )
        ax.add_patch(triangle)

    wheel_y = y_beam - SUPPORT_TRI_HEIGHT - SUPPORT_WHEEL_GAP
    for off in ROLLER_WHEEL_OFFSETS:
        ax.add_patch(Circle((1.0 + off, wheel_y), SUPPORT_WHEEL_RADIUS,
                            fill=False, linewidth=1.2, transform=trans, clip_on=False, zorder=6))
    base_y = wheel_y - SUPPORT_WHEEL_RADIUS - SUPPORT_BASE_EXTRA
    ax.plot([1.0 - SUPPORT_BASE_HALF_WIDTH, 1.0 + SUPPORT_BASE_HALF_WIDTH], [base_y, base_y],
            color="black", linewidth=1.4, transform=trans, clip_on=False, zorder=6)

    ax.text(0.0, PIN_LABEL_Y, "A", ha="center", fontsize=10, transform=ax.transData)
    ax.text(0.0, TYPE_LABEL_Y, "pin", ha="center", fontsize=9, transform=ax.transData)
    ax.text(L, PIN_LABEL_Y, "B", ha="center", fontsize=10, transform=ax.transData)
    ax.text(L, TYPE_LABEL_Y, "roller", ha="center", fontsize=9, transform=ax.transData)


def make_loading_figure(beam: BeamProblem):
    L = beam.L
    d = beam.dist_load

    fig, ax = plt.subplots(figsize=(11, 3.7))
    ax.set_xlim(0.0, L)
    ax.set_ylim(*BEAM_YLIM)
    ax.margins(x=0)

    ax.plot([0.0, L], [0.0, 0.0], lw=5, color="#3a4f7a", solid_capstyle="round", zorder=5)
    draw_supports(ax, L)

    x_load = np.linspace(d.a, d.b, 19)
    y_top = 0.55 + 1.55 * d.intensity(x_load) / max(d.w1, d.w2, 10.0)
    ax.plot([d.a, d.b], [y_top[0], y_top[-1]], lw=1.7, color="#c03fb5", zorder=3)
    for xi, yi in zip(x_load, y_top):
        ax.annotate("", xy=(xi, 0.05), xytext=(xi, yi),
                    arrowprops=dict(arrowstyle="->", lw=1.15, color="#c03fb5"), zorder=3)

    point_used = [[] for _ in POINT_LABEL_LEVELS]
    dist_used = [[] for _ in DIST_LABEL_LEVELS]
    moment_used = [[] for _ in MOMENT_LABEL_LEVELS]

    dist_y = label_lane(0.5 * (d.a + d.b), dist_used, DIST_LABEL_LEVELS, 2.5)
    ax.text(0.5 * (d.a + d.b), dist_y, beam.dist_label, ha="center", fontsize=9.6,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.8", alpha=0.96), zorder=10)
    ax.text(d.a, DIST_END_LABEL_Y, fmt_num(d.w1,1), ha="center", va="bottom", fontsize=8.8, color="#8a1484",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88), zorder=10)
    ax.text(d.b, DIST_END_LABEL_Y, fmt_num(d.w2,1), ha="center", va="bottom", fontsize=8.8, color="#8a1484",
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.88), zorder=10)

    p = beam.point_load
    ax.annotate("", xy=(p.x, 0.08), xytext=(p.x, 2.25),
                arrowprops=dict(arrowstyle="->", lw=2.0, color="#cc2222"), zorder=4)
    p_level = label_lane(p.x, point_used, POINT_LABEL_LEVELS, 1.7)
    ax.text(p.x, p_level, beam.point_label, ha="center", fontsize=9.8,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.8", alpha=0.96), zorder=11)

    m = beam.free_moment
    r = 0.33
    theta = np.linspace(np.deg2rad(35), np.deg2rad(325), 220) if m.M >= 0 else np.linspace(np.deg2rad(325), np.deg2rad(35), 220)
    cx, cy = m.x, 0.92
    ax.plot(cx + r * np.cos(theta), cy + r * np.sin(theta), lw=1.8, color="#0a7f7f", zorder=4)
    ax.annotate("", xy=(cx + r*np.cos(theta[-1]), cy + r*np.sin(theta[-1])),
                xytext=(cx + r*np.cos(theta[-8]), cy + r*np.sin(theta[-8])),
                arrowprops=dict(arrowstyle="->", lw=1.8, color="#0a7f7f"), zorder=4)
    m_level = label_lane(m.x, moment_used, MOMENT_LABEL_LEVELS, 2.0)
    ax.text(m.x, m_level, beam.moment_label, ha="center", fontsize=9.2,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.8", alpha=0.96), zorder=11)

    ax.set_xticks(np.arange(0.0, L + 0.001, 1.0))
    ax.grid(True, axis="x", alpha=0.25)
    ax.set_yticks([])
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Beam Loading Diagram", fontsize=13, pad=12)
    fig.tight_layout()
    return fig


def make_solution_figure(beam: BeamProblem):
    fig, (ax_v, ax_m) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True)

    x = beam.x
    ax_v.axhline(0.0, lw=1.0, color="black")
    ax_v.plot(x, beam.V, lw=2.2)
    ax_v.fill_between(x, 0.0, beam.V, alpha=0.15)
    ax_v.set_ylabel("V (kN)")
    ax_v.set_title("Shear Diagram")
    ax_v.grid(True, axis="y", alpha=0.3)
    ax_v.grid(True, axis="x", alpha=0.25)

    ax_m.axhline(0.0, lw=1.0, color="black")
    ax_m.plot(x, beam.M, lw=2.2)
    ax_m.fill_between(x, 0.0, beam.M, alpha=0.15)
    ax_m.set_ylabel("M (kN·m)")
    ax_m.set_title("Moment Diagram")
    ax_m.set_xlabel("x (m)")
    ax_m.grid(True, axis="y", alpha=0.3)
    ax_m.grid(True, axis="x", alpha=0.25)

    for xi, style in (
        (beam.point_load.x, "--"),
        (beam.free_moment.x, ":"),
        (beam.dist_load.a, "-."),
        (beam.dist_load.b, "-."),
    ):
        ax_v.axvline(xi, ls=style, lw=0.8, alpha=0.5)
        ax_m.axvline(xi, ls=style, lw=0.8, alpha=0.5)

    ax_m.set_xlim(0.0, beam.L)
    ax_m.set_xticks(np.arange(0.0, beam.L + 0.001, 1.0))
    fig.tight_layout()
    return fig


def fig_to_png_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def prepare_problem_bundle() -> None:
    beam = generate_problem()
    st.session_state.beam = beam
    st.session_state.loading_png = fig_to_png_bytes(make_loading_figure(beam))
    st.session_state.solution_png = fig_to_png_bytes(make_solution_figure(beam))
    st.session_state.show_solution = False
    st.session_state.show_reactions = False


if "beam" not in st.session_state:
    prepare_problem_bundle()

if "show_solution" not in st.session_state:
    st.session_state.show_solution = False

if "show_reactions" not in st.session_state:
    st.session_state.show_reactions = False

st.title("Random Beam Practice Generator")

c1, c2, c3 = st.columns([1, 1, 1])
with c1:
    st.button("New Problem", use_container_width=True, on_click=prepare_problem_bundle)
with c2:
    st.toggle("Show support reactions", key="show_reactions")
with c3:
    st.toggle("Show solution", key="show_solution")

beam = st.session_state.beam

st.image(st.session_state.loading_png, use_container_width=True)

reaction_slot = st.container()
with reaction_slot:
    st.markdown("### Support reactions")
    if st.session_state.show_reactions:
        st.latex(beam.reaction_latex)
    else:
        st.caption("Toggle on to reveal the support reactions.")

if st.session_state.show_solution:
    st.image(st.session_state.solution_png, use_container_width=True)

    st.subheader("Piecewise equations")
    for i, sec in enumerate(beam.sections, start=1):
        with st.container(border=True):
            st.markdown(f"**Section {i}:** {sec.range_text}")
            st.latex(sec.v_latex)
            st.latex(sec.m_latex)
else:
    st.info("Use **Show support reactions** and **Show solution** independently, depending on what you want to practice.")
