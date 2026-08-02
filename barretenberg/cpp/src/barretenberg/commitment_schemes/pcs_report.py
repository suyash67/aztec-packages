#!/usr/bin/env python3
"""Render pcs_report.html from a google-benchmark JSON dump of whir_honk_bench.

    ./build-arm64/bin/whir_honk_bench --benchmark_min_time=1x \
        --benchmark_format=json --benchmark_out=sweep.json
    python3 pcs_report.py sweep.json pcs_bench_results.json pcs_report.html

Backends are faceted by family so no chart carries more than six categorical hues.
"""

import json
import math
import re
import sys

# Display name, family, and the hue slot each backend keeps across every chart.
BACKENDS = [
    ("ultra_honk_kzg", None, "KZG (Gemini+Shplonk)", "pairing", 0),
    ("mercury_honk", None, "Mercury", "pairing", 1),
    ("vela_honk", None, "Vela", "pairing", 2),
    ("chopin_honk", None, "CHOPIN", "pairing", 3),
    ("kzh_honk", None, "KZH2", "pairing", 4),
    ("kzh3_honk", None, "KZH3", "pairing", 5),
    ("whir_honk", "Blake3sMerkleHasher", "WHIR", "hash", 0),
    ("whir_honk_repaired", "Blake3sMerkleHasher", "WHIR (repaired)", "hash", 1),
    ("ligero_honk", "Blake3sMerkleHasher", "Ligero", "hash", 2),
    ("switchfold_honk", "Blake3sMerkleHasher", "SwitchFold", "hash", 3),
    ("ipa_honk", None, "IPA", "dl", 0),
    ("hyrax_honk", None, "Hyrax", "dl", 1),
    ("dory_honk", None, "Dory", "dl", 2),
]

FAMILIES = [
    ("pairing", "Pairing-based (structured setup)"),
    ("hash", "Hash-based (transparent)"),
    ("dl", "Discrete-log (transparent)"),
]

METRICS = [
    ("prove", "Prover time", "ms", "Full proof construction; circuit building and key setup excluded."),
    ("verify", "Verifier time", "ms", "Full verification of the proof produced at the same size."),
    ("proof", "Proof size", "KiB", "Whole UltraHonk proof: commitments, sumcheck and the opening."),
]

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]

NAME_RE = re.compile(r"^([a-z0-9_]+)(?:<([A-Za-z0-9]+)>)?/(\d+)$")


def parse(path):
    """{(key, hasher): {metric: {log_n: value}}} plus the benchmark context."""
    with open(path) as handle:
        blob = json.load(handle)
    out = {}
    # Only the prove benchmarks report the realised log_n; verify rows carry just the target size,
    # so resolve them through the target -> realised map the prove rows establish.
    realised = {}
    for entry in blob.get("benchmarks", []):
        match = NAME_RE.match(entry["name"])
        if match and entry.get("log_n"):
            realised[int(match.group(3))] = int(entry["log_n"])

    for entry in blob.get("benchmarks", []):
        match = NAME_RE.match(entry["name"])
        if not match:
            continue
        stem, hasher, target = match.groups()
        for suffix, metric in (("_prove", "prove"), ("_verify", "verify")):
            if stem.endswith(suffix):
                key = (stem[: -len(suffix)], hasher)
                log_n = int(entry.get("log_n") or realised.get(int(target), 0))
                if not log_n:
                    continue
                series = out.setdefault(key, {})
                series.setdefault(metric, {})[log_n] = entry["real_time"]
                if metric == "prove" and entry.get("proof_KiB"):
                    series.setdefault("proof", {})[log_n] = entry["proof_KiB"]
                break
    return out, blob.get("context", {})


def nice_ticks(lo, hi):
    """Decade ticks spanning [lo, hi] on a log scale."""
    start = math.floor(math.log10(lo))
    stop = math.ceil(math.log10(hi))
    return [10.0**e for e in range(int(start), int(stop) + 1)]


def fmt(value):
    if value >= 1000:
        return f"{value:,.0f}"
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


def chart(series, title, unit, sizes):
    """One log-y line chart. `series` is [(label, colour, {log_n: value})]."""
    width, height = 430, 268
    left, right, top, bottom = 52, 104, 14, 34
    plot_w, plot_h = width - left - right, height - top - bottom

    values = [v for _, _, points in series for v in points.values() if v > 0]
    if not values:
        return ""
    lo, hi = min(values), max(values)
    ticks = nice_ticks(lo, hi)
    ymin, ymax = math.log10(ticks[0]), math.log10(ticks[-1])
    span = max(ymax - ymin, 1e-9)

    def sx(log_n):
        if len(sizes) == 1:
            return left + plot_w / 2
        return left + plot_w * (log_n - sizes[0]) / (sizes[-1] - sizes[0])

    def sy(value):
        return top + plot_h * (1 - (math.log10(value) - ymin) / span)

    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{title}">']
    for tick in ticks:
        y = sy(tick)
        parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}"/>')
        label = fmt(tick) if tick >= 1 else f"{tick:g}"
        parts.append(f'<text class="tick" x="{left - 7}" y="{y + 3.5:.1f}" text-anchor="end">{label}</text>')
    parts.append(
        f'<line class="axis" x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}"/>'
    )
    for log_n in sizes:
        x = sx(log_n)
        parts.append(
            f'<text class="tick" x="{x:.1f}" y="{top + plot_h + 15}" text-anchor="middle">2^{log_n}</text>'
        )
    parts.append(
        f'<text class="axislabel" x="{left + plot_w / 2:.0f}" y="{height - 4}" '
        f'text-anchor="middle">circuit size</text>'
    )
    parts.append(
        f'<text class="axislabel" x="6" y="{top + 8}" text-anchor="start">{unit}</text>'
    )

    labels = []
    for label, colour, points in series:
        ordered = [(n, points[n]) for n in sizes if n in points and points[n] > 0]
        if not ordered:
            continue
        path = " ".join(
            f"{'M' if i == 0 else 'L'}{sx(n):.1f},{sy(v):.1f}" for i, (n, v) in enumerate(ordered)
        )
        parts.append(f'<path d="{path}" fill="none" stroke="{colour}" stroke-width="2" stroke-linejoin="round"/>')
        for n, v in ordered:
            parts.append(
                f'<circle class="mark" cx="{sx(n):.1f}" cy="{sy(v):.1f}" r="4.2" fill="{colour}" '
                f'stroke="var(--surface)" stroke-width="2" '
                f'data-tip="{label} — 2^{n}: {fmt(v)} {unit}"/>'
            )
        last_n, last_v = ordered[-1]
        labels.append([sy(last_v), sx(last_n), label, colour])

    # Series whose endpoints nearly coincide would otherwise print on top of each other.
    labels.sort()
    for i in range(1, len(labels)):
        labels[i][0] = max(labels[i][0], labels[i - 1][0] + 12.5)
    for y, x, label, colour in labels:
        parts.append(
            f'<text class="dlabel" x="{x + 9:.1f}" y="{y + 4:.1f}" fill="{colour}">{label}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def main():
    sweep_path, json_out, html_out = sys.argv[1], sys.argv[2], sys.argv[3]
    data, context = parse(sweep_path)

    sizes = sorted({n for series in data.values() for pts in series.values() for n in pts})
    resolved = [(key, hasher, label, family, slot) for key, hasher, label, family, slot in BACKENDS
                if (key, hasher) in data]

    # Raw data alongside the report, keyed by display name.
    dump = {
        "context": {k: context.get(k) for k in ("date", "host_name", "num_cpus", "mhz_per_cpu", "library_build_type")},
        "sizes": sizes,
        "backends": {
            label: {metric: {str(n): v for n, v in sorted(data[(key, hasher)].get(metric, {}).items())}
                    for metric in ("prove", "verify", "proof")}
            for key, hasher, label, _family, _slot in resolved
        },
        "all_benchmarks": {f"{key}{'<' + hasher + '>' if hasher else ''}": {
            metric: {str(n): v for n, v in sorted(points.items())} for metric, points in series.items()
        } for (key, hasher), series in sorted(data.items(), key=lambda kv: (kv[0][0], kv[0][1] or ""))},
    }
    with open(json_out, "w") as handle:
        json.dump(dump, handle, indent=1, sort_keys=True)
        handle.write("\n")

    def best(metric, at):
        pool = [(data[(k, h)].get(metric, {}).get(at), lbl)
                for k, h, lbl, _f, _s in resolved if data[(k, h)].get(metric, {}).get(at)]
        return min(pool) if pool else None

    big = sizes[-1]
    tiles = []
    for metric, unit, caption in (("proof", "KiB", "smallest proof"), ("verify", "ms", "fastest verifier"),
                                  ("prove", "ms", "fastest prover")):
        winner = best(metric, big)
        if winner is None:
            continue
        value, label = winner
        tiles.append((f"{fmt(value)} {unit}", f"{caption} at 2^{big} — {label}"))
    tiles.append((str(len(resolved)), f"backends compared, 2^{sizes[0]}–2^{big}"))

    sections = []
    for metric, mlabel, unit, blurb in METRICS:
        figures = []
        for family, fam_label in FAMILIES:
            series = [(label, PALETTE[slot], data[(key, hasher)].get(metric, {}))
                      for key, hasher, label, fam, slot in resolved if fam == family]
            svg = chart(series, f"{mlabel}, {fam_label}", unit, sizes)
            if svg:
                figures.append(
                    f'<figure><figcaption><strong>{fam_label}</strong></figcaption>{svg}</figure>'
                )
        sections.append(
            f'<h2>{mlabel} <span class="unit">({unit}, log scale)</span></h2>'
            f'<p class="sub">{blurb}</p><div class="grid2">{"".join(figures)}</div>'
        )

    header = "".join(f"<th>2^{n}</th>" for n in sizes)
    rows = []
    for key, hasher, label, _family, _slot in resolved:
        for i, (metric, mlabel, unit, _blurb) in enumerate(METRICS):
            points = data[(key, hasher)].get(metric, {})
            cells = "".join(f"<td>{fmt(points[n]) if n in points else '—'}</td>" for n in sizes)
            name = f'<td class="bname">{label}</td>' if i == 0 else "<td></td>"
            cls = ' class="p2"' if i else ""
            rows.append(f'<tr{cls}>{name}<td class="metric">{mlabel} ({unit})</td>{cells}</tr>')

    host = context.get("host_name", "unknown host")
    cpus = context.get("num_cpus", "?")
    date = context.get("date", "")
    build = context.get("library_build_type", "")

    tile_html = "".join(f'<div class="tile"><div class="tv">{v}</div><div class="tl">{l}</div></div>'
                        for v, l in tiles)

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>UltraHonk PCS backends — consolidated comparison</title>
<style>
  .viz-root {{
    color-scheme: light;
    --surface:#fcfcfb; --plane:#f9f9f7; --ink1:#0b0b0b; --ink2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --axis:#c3c2b7; --ring:rgba(11,11,11,.10);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif; color:var(--ink1);
    background:var(--plane); margin:0; padding:28px 20px 60px;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:where(:not([data-theme="light"])) .viz-root {{
      color-scheme: dark;
      --surface:#1a1a19; --plane:#0d0d0d; --ink1:#ffffff; --ink2:#c3c2b7; --muted:#898781;
      --grid:#2c2c2a; --axis:#383835; --ring:rgba(255,255,255,.10);
    }}
  }}
  :root[data-theme="dark"] .viz-root {{
    color-scheme: dark;
    --surface:#1a1a19; --plane:#0d0d0d; --ink1:#ffffff; --ink2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --axis:#383835; --ring:rgba(255,255,255,.10);
  }}
  .wrap {{ max-width: 1320px; margin: 0 auto; }}
  h1 {{ font-size: 26px; margin: 0 0 6px; }}
  h2 {{ font-size: 19px; margin: 40px 0 4px; }}
  h2 .unit {{ font-weight: 400; color: var(--ink2); font-size: 15px; }}
  .sub {{ color: var(--ink2); margin: 0 0 18px; max-width: 78ch; font-size: 13.5px; }}
  .lede {{ color: var(--ink2); margin-bottom: 22px; max-width: 78ch; }}
  .tiles {{ display:flex; flex-wrap:wrap; gap:12px; margin: 20px 0 8px; }}
  .tile {{ background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:14px 20px; min-width:150px; }}
  .tv {{ font-size:24px; font-weight:650; }}
  .tl {{ font-size:12.5px; color:var(--ink2); margin-top:2px; }}
  .grid2 {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(430px, 1fr)); gap:18px; }}
  figure {{ background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:16px 16px 8px; margin:0; }}
  figcaption {{ font-size:13px; color:var(--ink2); margin-bottom:8px; }}
  figcaption strong {{ color:var(--ink1); font-size:14px; }}
  svg {{ width:100%; height:auto; display:block; overflow:visible; }}
  .grid {{ stroke:var(--grid); stroke-width:1; }}
  .axis {{ stroke:var(--axis); stroke-width:1; }}
  .tick {{ fill:var(--muted); font-size:11px; font-variant-numeric: tabular-nums; }}
  .axislabel {{ fill:var(--ink2); font-size:11.5px; }}
  .dlabel {{ font-size:11px; }}
  .mark {{ cursor:pointer; }}
  .mark:hover {{ stroke:var(--ink1); }}
  #tip {{ position:fixed; display:none; background:var(--ink1); color:var(--surface); font-size:12px;
         padding:6px 9px; border-radius:6px; pointer-events:none; z-index:10; max-width:280px; }}
  .tablewrap {{ overflow-x:auto; background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:8px; }}
  table {{ border-collapse:collapse; width:100%; font-size:13px; }}
  th, td {{ text-align:right; padding:6px 12px; border-bottom:1px solid var(--grid); font-variant-numeric: tabular-nums; }}
  th {{ color:var(--ink2); font-weight:600; }}
  td.bname, th:first-child {{ text-align:left; font-weight:600; white-space:nowrap; }}
  td.metric {{ text-align:left; color:var(--ink2); }}
  tr.p2 td {{ color:var(--ink2); }}
  footer {{ color:var(--muted); font-size:12px; margin-top:34px; }}
</style>
</head>
<body>
<div class="viz-root"><div class="wrap">
<h1>UltraHonk PCS backends — consolidated comparison</h1>
<p class="lede">Thirteen polynomial commitment backends proving and verifying the <em>same</em> UltraHonk
circuits through the shared <code>TransparentHonk</code> shell. Arithmetization, trace layout, relations
and sumcheck are identical across every row; only the commitment scheme and its opening argument change.
Charts are faceted by family so that no chart carries more than six series; hues are fixed per backend.
Hash-based rows use the Blake3s Merkle hasher — the Poseidon2 variants are in the raw JSON.</p>
<div class="tiles">{tile_html}</div>
{"".join(sections)}
<h2>All measurements</h2>
<p class="sub">Every backend at every size. Values are single-run wall clock, so treat differences under
a few percent as noise.</p>
<div class="tablewrap"><table>
<thead><tr><th>Backend</th><th>Metric</th>{header}</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table></div>
<footer>{host} · {cpus} CPUs · {build} build · {date}. Generated by <code>pcs_report.py</code> from
<code>whir_honk_bench</code>; raw data in <code>pcs_bench_results.json</code>.</footer>
</div></div>
<div id="tip"></div>
<script>
(function () {{
  var tip = document.getElementById('tip');
  document.addEventListener('mouseover', function (e) {{
    var t = e.target.getAttribute && e.target.getAttribute('data-tip');
    if (!t) return;
    tip.textContent = t; tip.style.display = 'block';
  }});
  document.addEventListener('mousemove', function (e) {{
    if (tip.style.display !== 'block') return;
    tip.style.left = Math.min(e.clientX + 14, window.innerWidth - 290) + 'px';
    tip.style.top = (e.clientY + 16) + 'px';
  }});
  document.addEventListener('mouseout', function (e) {{
    if (e.target.getAttribute && e.target.getAttribute('data-tip')) tip.style.display = 'none';
  }});
}})();
</script>
</body>
</html>
"""
    with open(html_out, "w") as handle:
        handle.write(html)
    print(f"{len(resolved)} backends, sizes {sizes} -> {html_out}, {json_out}")


if __name__ == "__main__":
    main()
