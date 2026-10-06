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


def chart(series, metric, title, unit_fmt, log_y, counts, ylabel):
    """series: list of (key, label, cls, dashed, {actions: value})."""
    W, H = 640, 330
    L, R, T, B = 64, 20, 18, 46
    xs = counts
    vals = [v for *_, d in series for a, v in d.items() if a in xs]
    if not vals:
        return ""
    lo, hi = min(vals), max(vals)
    if log_y:
        lo_e = math.floor(math.log10(lo) * 2) / 2
        hi_e = math.ceil(math.log10(hi) * 2) / 2
        ticks = []
        e = math.floor(lo_e)
        while e <= hi_e + 1e-9:
            for m in (1, 2, 5):
                t = m * 10**e
                if 10**lo_e * 0.999 <= t <= 10**hi_e * 1.001:
                    ticks.append(t)
            e += 1
        ymap = lambda v: T + (H - T - B) * (hi_e - math.log10(v)) / (hi_e - lo_e)
    else:
        top = hi * 1.1
        step = 10 ** math.floor(math.log10(top / 4))
        for m in (1, 2, 5, 10):
            if top / (m * step) <= 6:
                step *= m
                break
        ticks = [i * step for i in range(int(top / step) + 1)]
        top = ticks[-1] if ticks[-1] >= hi else ticks[-1] + step
        ymap = lambda v: T + (H - T - B) * (1 - v / top)
    xmap = lambda a: L + (W - L - R) * (math.log2(a) / math.log2(xs[-1]) if xs[-1] > 1 else 0)
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{title}">']
    for t in ticks:
        y = ymap(t)
        out.append(f'<line class="grid" x1="{L}" x2="{W - R}" y1="{y:.1f}" y2="{y:.1f}"/>')
        out.append(f'<text class="tick" x="{L - 8}" y="{y + 4:.1f}" text-anchor="end">{unit_fmt(t)}</text>')
    for a in xs:
        x = xmap(a)
        out.append(f'<line class="grid faint" x1="{x:.1f}" x2="{x:.1f}" y1="{T}" y2="{H - B}"/>')
        out.append(f'<text class="tick" x="{x:.1f}" y="{H - B + 18}" text-anchor="middle">{a}</text>')
    out.append(f'<text class="axis" x="{(L + W - R) / 2}" y="{H - 6}" text-anchor="middle">Actions per proof</text>')
    out.append(
        f'<text class="axis" transform="translate(14 {(T + H - B) / 2}) rotate(-90)" text-anchor="middle">{ylabel}</text>'
    )
    for key, label, cls, dashed, d in series:
        pts = [(xmap(a), ymap(d[a])) for a in xs if a in d]
        if not pts:
            continue
        path = " ".join(f"{'M' if i == 0 else 'L'}{x:.1f},{y:.1f}" for i, (x, y) in enumerate(pts))
        out.append(f'<path class="line {cls}{" dashed" if dashed else ""}" d="{path}"><title>{label}</title></path>')
        for (x, y), a in zip(pts, [a for a in xs if a in d]):
            out.append(f'<circle class="dot {cls}" cx="{x:.1f}" cy="{y:.1f}" r="3.5"><title>{label}: {unit_fmt(d[a])} at {a}</title></circle>')
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
    html = html.replace(
        "{{CHART_PROVE}}", chart(series(mt, "prove"), "prove", "Prover time", fmt_time, True, counts_mt, "Prover time")
    )
    html = html.replace(
        "{{CHART_VERIFY}}", chart(series(mt, "verify"), "verify", "Verifier time", fmt_time, True, counts_mt, "Verifier time")
    )
    html = html.replace(
        "{{CHART_SIZE}}",
        chart(series(mt, "size"), "size", "Proof size", lambda b: f"{b / 1024:.0f} KiB", False, counts_mt, "Proof size"),
    )
    html = html.replace(
        "{{CHART_PROVE_ST}}",
        chart(series(st, "prove"), "prove", "Prover time, one thread", fmt_time, True, counts_st, "Prover time"),
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
