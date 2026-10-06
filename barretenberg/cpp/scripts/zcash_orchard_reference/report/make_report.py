#!/usr/bin/env python3
"""Builds the Orchard benchmark report (HTML with inline SVG charts) from the result files."""
import json
import math
import sys
from pathlib import Path

RESULTS = Path(sys.argv[1])
OUT = Path(sys.argv[2])

SYSTEMS = [
    # key, label, short description, css class, dashed
    ("zcash", "Zcash halo2 (production)", "orchard 0.16 · halo2_proofs 0.4 · PLONKish · Vesta IPA", "s-zcash", False),
    ("halo2-honk-pasta", "halo2 gates · Honk · Pasta", "custom Honk flavor with Orchard's halo2 gates · Vesta · halo2 IPA", "s-h2p", False),
    ("halo2-honk-bn254", "halo2 gates · Honk · BN254", "same gates ported to BN254/Grumpkin · KZG", "s-h2b", True),
    ("ultra-zk-pasta", "UltraHonk · Pasta", "hand-written Ultra gates over F_p · Vesta · halo2 IPA", "s-ulp", False),
    ("ultra-zk-bn254", "UltraHonk · BN254", "barretenberg stdlib over BN254/Grumpkin · KZG", "s-ulb", True),
]


def load_jsonl(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def load():
    data = {}
    for mode in ("mt", "st"):
        rows = {}
        z = RESULTS / f"zcash_{mode}.json"
        if z.exists():
            for r in json.loads(z.read_text()):
                rows.setdefault("zcash", {})[r["actions"]] = {
                    "prove": r["prove_s"],
                    "verify": r["verify_s"],
                    "size": r["proof_bytes"],
                    "keygen": r["keygen_s"],
                }
        for key, *_ in SYSTEMS[1:]:
            for r in load_jsonl(RESULTS / f"bb_{key}_{mode}.jsonl"):
                assert r["verified"], (key, r)
                rows.setdefault(key, {})[r["actions"]] = {
                    "prove": r["prove_total_s"],
                    "verify": r["verify_s"],
                    "size": r["proof_bytes_compressed"],
                    "size_raw": r["proof_bytes"],
                    "keygen": r["keygen_s"],
                    "rows": r["rows"],
                    "log_n": r["log_n"],
                    "synth": r["synthesize_s"],
                }
        data[mode] = rows
    return data


def fmt_time(s):
    if s < 1e-2:
        return f"{s * 1e3:.1f} ms"
    if s < 1:
        return f"{s * 1e3:.0f} ms"
    return f"{s:.2f} s"


def fmt_bytes(b):
    return f"{b / 1024:.1f} KiB"


SHORT_LABELS = {
    "zcash": "Zcash halo2",
    "halo2-honk-pasta": "halo2 gates · Pasta",
    "halo2-honk-bn254": "halo2 gates · BN254",
    "ultra-zk-pasta": "UltraHonk · Pasta",
    "ultra-zk-bn254": "UltraHonk · BN254",
}


def tick_time(s):
    if s < 1:
        ms = s * 1e3
        return f"{ms:g} ms"
    return f"{s:g} s"


def chart(series, title, log_y, counts, ylabel, tick_fmt, value_fmt):
    """series: list of (key, label, cls, dashed, {actions: value}). One line per system, labelled at its last point."""
    W, H = 960, 440
    L, R, T, B = 110, 190, 20, 58
    xs = counts
    vals = [v for *_, d in series for a, v in d.items() if a in xs]
    if not vals:
        return ""
    lo, hi = min(vals), max(vals)
    if log_y:
        nice = [m * 10**e for e in range(-4, 3) for m in (1, 2, 5)]
        ticks = [t for t in nice if lo / 1.6 <= t <= hi * 1.6]
        y_lo, y_hi = min(ticks[0], lo / 1.25), max(ticks[-1], hi * 1.25)
        ymap = lambda v: T + (H - T - B) * (math.log10(y_hi) - math.log10(v)) / (math.log10(y_hi) - math.log10(y_lo))
    else:
        top = hi * 1.08
        step = 10 ** math.floor(math.log10(top / 4))
        for m in (1, 2, 5, 10):
            if top / (m * step) <= 6:
                step *= m
                break
        ticks = [i * step for i in range(int(math.ceil(top / step)) + 1)]
        y_hi = ticks[-1]
        ymap = lambda v: T + (H - T - B) * (1 - v / y_hi)
    span = math.log2(xs[-1]) if xs[-1] > 1 else 1
    xmap = lambda a: L + 18 + (W - L - R - 36) * math.log2(a) / span
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{title}">']
    for t in ticks:
        y = ymap(t)
        out.append(f'<line class="grid" x1="{L}" x2="{W - R + 8}" y1="{y:.1f}" y2="{y:.1f}"/>')
        out.append(f'<text class="tick" x="{L - 10}" y="{y + 4.5:.1f}" text-anchor="end">{tick_fmt(t)}</text>')
    out.append(f'<line class="baseline" x1="{L}" x2="{W - R + 8}" y1="{H - B}" y2="{H - B}"/>')
    for a in xs:
        x = xmap(a)
        out.append(f'<line class="baseline" x1="{x:.1f}" x2="{x:.1f}" y1="{H - B}" y2="{H - B + 6}"/>')
        out.append(f'<text class="tick" x="{x:.1f}" y="{H - B + 24}" text-anchor="middle">{a}</text>')
    out.append(f'<text class="axis" x="{(L + W - R) / 2}" y="{H - 8}" text-anchor="middle">Actions per proof</text>')
    out.append(
        f'<text class="axis" transform="translate(20 {(T + H - B) / 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>'
    )
    # End-of-line labels, spread apart vertically so they never overlap.
    ends = []
    for key, label, cls, dashed, d in series:
        pts = [(a, xmap(a), ymap(d[a])) for a in xs if a in d]
        if not pts:
            continue
        path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (_, x, y) in enumerate(pts))
        out.append(f'<path class="line {cls}{" dashed" if dashed else ""}" d="{path}"><title>{label}</title></path>')
        for a, x, y in pts:
            out.append(
                f'<circle class="dot {cls}" cx="{x:.1f}" cy="{y:.1f}" r="4.5"><title>{label}: {value_fmt(d[a])} at {a} Actions</title></circle>'
            )
        a_last, x_last, y_last = pts[-1]
        ends.append([y_last, x_last, key, cls, value_fmt(d[a_last])])
    ends.sort()
    gap = 34
    for i in range(1, len(ends)):
        ends[i][0] = max(ends[i][0], ends[i - 1][0] + gap)
    overflow = ends[-1][0] - (H - B) if ends else 0
    if overflow > 0:
        for e in ends:
            e[0] -= overflow
        for i in range(len(ends) - 2, -1, -1):
            ends[i][0] = min(ends[i][0], ends[i + 1][0] - gap)
    for y, x, key, cls, v in ends:
        lx = W - R + 18
        out.append(f'<text class="endlabel {cls}" x="{lx}" y="{y - 2:.1f}">{SHORT_LABELS[key]}</text>')
        out.append(f'<text class="endvalue" x="{lx}" y="{y + 14:.1f}">{v}</text>')
    out.append("</svg>")
    return "\n".join(out)


def legend():
    items = []
    for key, label, desc, cls, dashed in SYSTEMS:
        items.append(
            f'<li><svg viewBox="0 0 28 10" aria-hidden="true"><line class="line {cls}{" dashed" if dashed else ""}" x1="1" x2="27" y1="5" y2="5"/></svg>'
            f"<span><b>{label}</b><small>{desc}</small></span></li>"
        )
    return '<ul class="legend">' + "".join(items) + "</ul>"


def table(data, metric, fmt, counts, highlight_min=True, extra=None):
    head = "".join(f"<th>{a}</th>" for a in counts)
    rows = []
    for key, label, desc, cls, dashed in SYSTEMS:
        d = data.get(key, {})
        cells = []
        best = {a: min((data[k][a][metric] for k in data if a in data[k]), default=None) for a in counts}
        for a in counts:
            if a in d:
                v = d[a][metric]
                mark = ' class="best"' if highlight_min and best[a] == v else ""
                cells.append(f"<td{mark}>{fmt(v)}</td>")
            else:
                cells.append("<td class='na'>—</td>")
        rows.append(f'<tr><th scope="row"><i class="sw {cls}"></i>{label}</th>{"".join(cells)}</tr>')
    return (
        '<div class="tablewrap"><table><thead><tr><th scope="col">Actions</th>'
        + head
        + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
    )


def size_table(data, counts):
    head = "".join(f"<th>{a}</th>" for a in counts)
    rows = []
    for key, label, desc, cls, dashed in SYSTEMS[1:]:
        d = data.get(key, {})
        cells = []
        for a in counts:
            if a in d:
                cells.append(f"<td>{d[a]['rows']:,}<small>2<sup>{d[a]['log_n']}</sup></small></td>")
            else:
                cells.append("<td class='na'>—</td>")
        rows.append(f'<tr><th scope="row"><i class="sw {cls}"></i>{label}</th>{"".join(cells)}</tr>')
    return (
        '<div class="tablewrap"><table class="sizes"><thead><tr><th scope="col">Actions</th>'
        + head
        + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
    )


def main():
    data = load()
    mt, st = data["mt"], data["st"]
    counts_mt = sorted({a for d in mt.values() for a in d})
    counts_st = sorted({a for d in st.values() for a in d})

    def series(d, metric):
        return [(k, lab, cls, dashed, {a: v[metric] for a, v in d.get(k, {}).items()}) for k, lab, _, cls, dashed in SYSTEMS]

    tmpl = (Path(__file__).parent / "report_template.html").read_text()
    html = tmpl
    html = html.replace("{{LEGEND}}", legend())
    def in_kib(d):
        return [(k, lab, cls, dashed, {a: v / 1024 for a, v in vals.items()}) for k, lab, cls, dashed, vals in d]
    html = html.replace(
        "{{CHART_PROVE}}", chart(series(mt, "prove"), "Prover time", True, counts_mt, "Prover time (log scale)", tick_time, fmt_time)
    )
    html = html.replace(
        "{{CHART_VERIFY}}", chart(series(mt, "verify"), "Verifier time", True, counts_mt, "Verifier time (log scale)", tick_time, fmt_time)
    )
    html = html.replace(
        "{{CHART_SIZE}}", chart(in_kib(series(mt, "size")), "Proof size", False, counts_mt, "Proof size", lambda t: f"{t:g} KiB", lambda v: f"{v:.1f} KiB")
    )
    html = html.replace(
        "{{CHART_PROVE_ST}}",
        chart(series(st, "prove"), "Prover time, one thread", True, counts_st, "Prover time (log scale)", tick_time, fmt_time),
    )
    html = html.replace("{{TABLE_PROVE}}", table(mt, "prove", fmt_time, counts_mt))
    html = html.replace("{{TABLE_VERIFY}}", table(mt, "verify", fmt_time, counts_mt))
    html = html.replace("{{TABLE_SIZE}}", table(mt, "size", fmt_bytes, counts_mt))
    html = html.replace("{{TABLE_PROVE_ST}}", table(st, "prove", fmt_time, counts_st))
    html = html.replace("{{TABLE_KEYGEN}}", table(mt, "keygen", fmt_time, counts_mt))
    html = html.replace("{{TABLE_ROWS}}", size_table(mt, counts_mt))
    for key, *_ in SYSTEMS:
        for a in counts_mt:
            for m in ("prove", "verify", "size"):
                v = mt.get(key, {}).get(a, {}).get(m)
                if v is not None:
                    html = html.replace("{{" + f"{key}:{a}:{m}" + "}}", fmt_bytes(v) if m == "size" else fmt_time(v))
            r = mt.get(key, {}).get(a, {}).get("rows")
            if r is not None:
                html = html.replace("{{" + f"{key}:{a}:rows_n" + "}}", f"{r:,}")
        for a in counts_st:
            v = st.get(key, {}).get(a, {}).get("prove")
            if v is not None:
                html = html.replace("{{" + f"{key}:{a}:prove_st" + "}}", fmt_time(v))
    def ratio(k1, k2, a, m, d=mt):
        return f"{d[k1][a][m] / d[k2][a][m]:.1f}×"
    import re
    html = re.sub(
        r"\{\{ratio:([a-z0-9-]+):([a-z0-9-]+):(\d+):([a-z_]+)\}\}",
        lambda mm: ratio(mm.group(1), mm.group(2), int(mm.group(3)), mm.group(4).replace("_st", ""), st if mm.group(4).endswith("_st") else mt),
        html,
    )
    left = re.findall(r"\{\{[^}]*\}\}", html)
    assert not left, left
    OUT.write_text(html)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
