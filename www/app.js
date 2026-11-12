"use strict";
/* FCC-ee readout rates: static page, no build step. Data is fetched from data/*.json next to the page.

   Model (same as tools/calc_ref.py):
     Gbps/layer = occupancy * rate[MHz] * bits * multiplier / 1000 * n_cells * SF
     links      = Gbps / (line rate * encoding efficiency)
     installed  = fibers_per_module * n_cells / cells_per_module                                   */

const $ = (s, r) => (r || document).querySelector(s);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

const SRC = { ipc: "IPC", z: "Z decays", inj: "Injection" };
const SRC_VAR = { ipc: "--s-ipc", z: "--s-z", inj: "--s-inj" };
const UNIT_NAME = { fiber_5g12_8b10b: "fibers", gbe_poe_1g: "ports", sipho_100g: "links" };
const FLAG_TEXT = {
  partial_layers: "Only some of the layers are included.",
  no_geometry: "No geometry yet, so no rate can be computed.",
  multiplier_inherited_from_sheet: "Multiplier (3) inherited from the pixel cluster size in the old sheet; to be set for this detector type.",
  no_stats: "Single-event sample: no statistics.",
  capacity_unreliable: "Installed capacity is not defined yet (the sheet assumed 9 cells per module for every detector), so no capacity is shown.",
};

const TIPS = {
  mean: "<b>Mean</b><br>Occupancy per event averaged over all simulated events: the typical load.",
  p95: "<b>95th percentile</b><br>Occupancy that 95% of events stay below. Sizes for busy events instead of the typical one.",
  peak: "<b>Peak</b><br>Continuous sources plus injection at its peak rate, as if they all coincide. The conservative choice for sizing links.",
  det: "<b>Detector</b><br>Detector design to show. ALLEGRO is the main one; IDEA only has the vertex so far.",
  ipcsample: "<b>IPC sample</b><br>Simulated IPC sample (beam lattice version) the occupancies come from. LCC_V106 is the current default.",
  rateipc: "<b>IPC rate</b><br>Events per second fed to the IPC sample, in MHz. Default is the bunch-crossing rate of the working point (52 MHz at the Z-pole): one event per crossing.",
  wp: "<b>Working point</b><br>Collider energy the samples were simulated for. Background and physics rates depend on it. Only the Z-pole has samples so far.",
  zsample: "<b>Z sample</b><br>Which simulated Z decays the occupancy comes from: the sheet's mixture of Zbb and Zmumu proxies (weights below 1 because invisible decays contribute nothing), or one sample alone counted as 100% of Z decays.",
  injsample: "<b>Injection sample</b><br>Simulated injection-loss sample. Core and halo are separate components once their samples exist.",
  useipc: "<b>IPC</b><br>Include or exclude the incoherent-pair-creation background. Choosing \"No sample\" in the list does the same.",
  fpm: "<b>Fibers per module</b><br>Readout fibers (links) installed per module. Empty keeps the value of each layer from the data. Layers can also be set individually in the Layers tab.",
  cpm: "<b>Cells per module</b><br>Cells read out by one module (9 in the old sheet for every subdetector). Installed fibers = fibers/module × cells ÷ cells/module.",
  usez: "<b>Z decays</b><br>Add the Z-decay physics signal (Z-pole working point) on top of the background. It is physics, not background, so it can be left out of the total.",
  ratez: "<b>Z decay rate</b><br>Z event rate in MHz (0.1 MHz = 100 kHz).",
  useinj: "<b>Injection</b><br>Add top-up injection background. Off by default because its scale factors are not known yet.",
  injpeak: "<b>Peak scale (example)</b><br>Effective event rate in MHz that turns the one-event injection sample into a peak rate. The average is this times the duty. Only an example value for now.",
  sf: "<b>Safety factor</b><br>Multiplies every rate to leave margin for uncertainties in the simulation and the design.",
  tech: "<b>Link technology</b><br>Link type used to turn bandwidth into a number of links. Entries marked unverified are placeholders.",
  tline: "<b>Line rate</b><br>Raw bit rate of one link, in Gbps.",
  teff: "<b>Encoding efficiency</b><br>Share of the line rate left for data after encoding overhead (8b10b gives 0.8).",
  tpow: "<b>Power per link</b><br>Power drawn by one link in watts. Used by the Power unit; unknown for now.",
  bits: "<b>Bits per hit</b><br>Bits sent for one hit of this subdetector.",
  mult: "<b>Multiplier</b><br>Factor on the hits: the cluster size for pixel detectors. The inherited value 3 still has to be set per detector type.",
  logscale: "<b>Log scale</b><br>Bar length on a logarithmic axis so small subdetectors stay visible. Segments then show shares of each bar.",
  avg: () => "<b>Average</b><br>Continuous sources plus injection weighted by its duty (" + (duty() * 100).toFixed(1) + "% of the time). The relevant one for average throughput and storage.",
};

const tipAttr = (k) => 'data-tip="' + esc(typeof TIPS[k] === "function" ? TIPS[k]() : TIPS[k]) + '"';

let D = null; // all data
let S = null; // UI state

/* ---------------------------------------------------------------- data + state */
async function loadData() {
  const get = (p) => fetch(p).then((r) => { if (!r.ok) throw new Error(p + ": " + r.status); return r.json(); });
  const [index, settings, sources, technologies, working_points] = await Promise.all([
    get("data/index.json"), get("data/settings.json"), get("data/sources.json"), get("data/technologies.json"), get("data/working_points.json"),
  ]);
  const detectors = {};
  for (const d of index.detectors) {
    const [detector, samples] = await Promise.all([get(`data/${d.id}/detector.json`), get(`data/${d.id}/samples.json`)]);
    detectors[d.id] = { detector, samples };
  }
  let facts = [];
  try { facts = (await fetch("data/penguin_facts.txt").then((r) => r.text())).split("\n").map((l) => l.trim()).filter(Boolean); } catch (e) { /* optional */ }
  return { index, settings, sources, technologies, working_points, detectors, facts };
}

function detData() { return D.detectors[S.det]; }
function samplesOf(source) { return detData().samples.samples.filter((s) => s.source === source); }
function zParts() {
  return S.zSel === "mixture" || !sampleById(S.zSel) ? D.sources.z.mixture : [{ sample: S.zSel, weight: 1 }];
}
function sampleDesc(s) {
  if (!s) return "";
  const p = [];
  if (s.process) p.push(s.process);
  if (s.lattice) p.push("lattice " + s.lattice);
  if (s.source === "injection" || s.source === "sr") p.push(s.component ? "component " + s.component : "component not specified (core/halo TBC)");
  if (s.n_events != null) p.push(s.n_events + " event" + (s.n_events === 1 ? "" : "s"));
  return p.join(" · ");
}
function sampleById(id) { return detData().samples.samples.find((s) => s.id === id); }

function freshState(detId, keepUi) {
  const st = {
    det: detId,
    tab: keepUi ? keepUi.tab : "overview",
    sub: null,
    unit: keepUi ? keepUi.unit : "gbps",
    log: keepUi ? keepUi.log : false,
    stat: "mean", basis: "peak",
    sf: D.settings.safety_factor,
    techId: D.settings.default_technology,
    tech: {},
    ipcSample: null, includeIpc: true, zSel: D.sources.z.default_selection || "mixture", injSample: null,
    includeZ: true, includeInj: false,
    rates: { ipc: wpData().f_bc_MHz, z: D.sources.z.rate_MHz },
    injPeak: 52,
    enabled: new Set(),
    model: {},
    layerFib: {},
    sort: { col: "total", dir: -1 },
    filter: "",
    csv: null,
    openSec: keepUi ? keepUi.openSec : {},
  };
  for (const [id, t] of Object.entries(D.technologies)) {
    st.tech[id] = { line: t.line_rate_Gbps, eff: t.encoding_efficiency, power: t.power_W_per_link };
  }
  const dd = D.detectors[detId];
  const ipc = dd.samples.samples.filter((s) => s.source === "ipc");
  st.ipcSample = (ipc.find((s) => s.default) || ipc[0] || {}).id || null;
  st.injSample = (dd.samples.samples.find((s) => s.source === "injection") || {}).id || null;
  for (const sub of dd.detector.subdetectors) {
    st.enabled.add(sub.id);
    st.model[sub.id] = { bits: sub.readout.bits_per_unit, mult: sub.readout.multiplier, fpm: null, cpm: sub.cells_per_module };
  }
  st.sub = dd.detector.subdetectors[0].id;
  return st;
}

/* ---------------------------------------------------------------- model */
function techP() { return S.tech[S.techId]; }
function linkCap() { const t = techP(); return t.line * t.eff; }
function wpData() { return D.working_points[D.settings.working_point]; }
function duty() {
  const p = D.sources.injection.pulse;
  return Math.min(1, p.turns_after_injection / wpData().f_rev_Hz / p.injection_period_s);
}
function occ(sample, subId, layerId, stat) {
  const o = sample && sample.occupancy[subId] && sample.occupancy[subId][String(layerId)];
  if (!o) return null;
  return stat === "p95" && o.p95 != null ? o.p95 : o.mean;
}
function toGbps(o, rate, m, cells) { return (o * rate * m.bits * m.mult) / 1000 * cells * S.sf; }

function compute() {
  const dd = detData();
  const cap = linkCap();
  const ipcS = S.includeIpc ? sampleById(S.ipcSample) : null;
  const hasZ = samplesOf("z").length > 0, hasInj = samplesOf("injection").length > 0;
  const useZ = S.includeZ && hasZ, useInj = S.includeInj && hasInj;
  const injS = sampleById(S.injSample);
  const injRate = S.basis === "peak" ? S.injPeak : S.injPeak * duty();
  const res = { subs: [], none: [], total: 0, links: 0, linksCeil: 0, installed: 0, over: 0, layers: 0, noIpc: 0 };

  for (const sub of dd.detector.subdetectors) {
    if (!S.enabled.has(sub.id)) continue;
    const m = S.model[sub.id];
    const rec = { sub, rows: [], parts: { ipc: 0, z: 0, inj: 0 }, total: 0, allowed: 0, over: 0, noIpc: 0 };
    // a user-set fiber count makes the capacity of a flagged subdetector meaningful
    const userCap = m.fpm != null || m.cpm !== sub.cells_per_module || sub.layers.some((l) => S.layerFib[sub.id + ":" + l.id] != null);
    for (const layer of sub.layers) {
      if (layer.n_cells == null) continue;
      const oi = occ(ipcS, sub.id, layer.id, S.stat);
      let oz = null;
      if (useZ) {
        let acc = 0, seen = false;
        for (const mx of zParts()) {
          const v = occ(sampleById(mx.sample), sub.id, layer.id, S.stat);
          if (v != null) { acc += mx.weight * v; seen = true; }
        }
        if (seen) oz = acc;
      }
      const oj = useInj ? occ(injS, sub.id, layer.id, "mean") : null;
      const parts = {
        ipc: oi == null ? 0 : toGbps(oi, S.rates.ipc, m, layer.n_cells),
        z: oz == null ? 0 : toGbps(oz, S.rates.z, m, layer.n_cells),
        inj: oj == null ? 0 : toGbps(oj, injRate, m, layer.n_cells),
      };
      const total = parts.ipc + parts.z + parts.inj;
      const capOk = sub.flags.indexOf("capacity_unreliable") < 0 || userCap;
      const fpm = S.layerFib[sub.id + ":" + layer.id] ?? m.fpm ?? layer.fibers_per_module;
      const inst = capOk ? (fpm * layer.n_cells) / m.cpm : 0;
      const allowed = inst * cap;
      const row = {
        sub, layer, parts, total, allowed, installed: inst, occ: { ipc: oi, z: oz, inj: oj },
        channels: layer.n_cells * sub.sensors_per_cell * sub.channels_per_sensor,
        over: capOk && total > allowed * 1.0000001, noIpc: oi == null && S.includeIpc,
      };
      rec.rows.push(row);
      for (const k of ["ipc", "z", "inj"]) rec.parts[k] += parts[k];
      rec.total += total; rec.allowed += allowed;
      if (row.over) rec.over++;
      if (row.noIpc) rec.noIpc++;
      res.links += total / cap; res.linksCeil += Math.ceil(total / cap - 1e-9);
      res.installed += inst; res.layers++;
    }
    const anyData = rec.rows.some((r) => r.total > 0);
    if (!rec.rows.length || !anyData) { res.none.push(rec); continue; }
    res.subs.push(rec);
    res.total += rec.total; res.over += rec.over; res.noIpc += rec.noIpc;
  }
  return res;
}

/* ---------------------------------------------------------------- formatting */
function short(x) {
  if (x === 0) return "0";
  const a = Math.abs(x);
  let s = a >= 1000 ? Math.round(x).toLocaleString("en-US") : a >= 100 ? x.toFixed(0) : a >= 10 ? x.toFixed(1) : a >= 1 ? x.toFixed(2) : x.toPrecision(2);
  return s.indexOf(".") >= 0 && s.indexOf("e") < 0 ? s.replace(/\.?0+$/, "") : s;
}
function unitName() { return UNIT_NAME[S.techId] || "links"; }
function powerOk() { return techP().power != null && techP().power > 0; }
/* Gbps -> display value in the selected unit (null when not computable) */
function conv(g) {
  if (S.unit === "links") return g / linkCap();
  if (S.unit === "power") return powerOk() ? (g / linkCap()) * techP().power : null;
  return g;
}
function fmtV(v) {
  if (v == null) return "set W/link";
  if (S.unit === "links") return short(v) + " " + unitName();
  if (S.unit === "power") return v >= 1000 ? short(v / 1000) + " kW" : short(v) + " W";
  return v >= 1000 ? short(v / 1000) + " Tbps" : short(v) + " Gbps";
}
const fmtG = (g) => fmtV(conv(g));
const pct = (x) => (x == null ? "–" : (x * 100 < 0.1 && x > 0 ? String(Number((x * 100).toPrecision(2))) : (x * 100).toFixed(2).replace(/\.?0+$/, "")) + "%");

/* ---------------------------------------------------------------- charts */
function niceTicks(max) {
  if (!(max > 0)) return { max: 1, ticks: [0, 0.5, 1] };
  const raw = max / 4, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  const step = (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * p;
  const n = Math.ceil(max / step - 1e-9);
  return { max: n * step, ticks: Array.from({ length: n + 1 }, (_, i) => i * step) };
}

function scaleFor(values) {
  const pos = values.filter((v) => v > 0);
  const vmax = Math.max(...pos, 0);
  if (S.log && pos.length) {
    const lo = Math.pow(10, Math.floor(Math.log10(Math.min(...pos))));
    let hi = Math.pow(10, Math.ceil(Math.log10(vmax)));
    if (hi <= lo) hi = lo * 10;
    const ticks = [];
    for (let v = lo; v <= hi * 1.0001; v *= 10) ticks.push(v);
    return { log: true, pos: (v) => (v <= 0 ? 0 : Math.max(0.015, Math.log(v / lo) / Math.log(hi / lo))), ticks };
  }
  const t = niceTicks(vmax);
  return { log: false, pos: (v) => Math.max(0, v / t.max), ticks: t.ticks };
}

function axisDiv(maxTick) {
  if (S.unit === "gbps") return maxTick >= 2000 ? { d: 1000, u: "Tbps" } : { d: 1, u: "Gbps" };
  if (S.unit === "power") return maxTick >= 2000 ? { d: 1000, u: "kW" } : { d: 1, u: "W" };
  return { d: 1, u: unitName() };
}

function tipFor(row, k) {
  const share = row.total > 0 ? Math.round((row.parts[k] / row.total) * 100) : 0;
  return `<b>${esc(row.label)}</b><br>${SRC[k]}: ${fmtG(row.parts[k])} (${share}%)<br>Total: ${fmtG(row.total)}`
    + (row.allowed ? `<br>Installed: ${fmtG(row.allowed)}` : "");
}

/* rows: {key,label,sub,parts,total,allowed,clickable} in Gbps */
function barChart(rows, opt) {
  opt = opt || {};
  if (S.unit === "power" && !powerOk()) {
    return `<div class="empty">Power needs a W/link value for ${esc(D.technologies[S.techId].label)}. Enter it in the sidebar under Link technology.</div>`;
  }
  const vals = rows.map((r) => conv(r.total));
  const sc = scaleFor(vals);
  const tmax = sc.ticks[sc.ticks.length - 1];
  const ad = axisDiv(tmax);
  const grid = '<div class="grid-lines">' + sc.ticks.map((t) => `<u style="left:${(sc.pos(t) * 100).toFixed(2)}%"></u>`).join("") + "</div>";
  const keys = ["ipc", "z", "inj"].filter((k) => rows.some((r) => r.parts[k] > 0));
  let h = '<div class="legend">'
    + keys.map((k) => `<span><i class="sw" style="background:var(${SRC_VAR[k]})"></i>${SRC[k]}</span>`).join("")
    + (opt.cap ? '<span><i class="tick"></i>Installed capacity</span>' : "")
    + (sc.log ? '<span class="muted">log scale: bar length is log, segments are shares</span>' : "") + "</div>";
  h += '<div class="chart">';
  for (const r of rows) {
    const v = conv(r.total);
    const len = sc.pos(v);
    const bar = r.total > 0
      ? `<div class="bar" style="width:${(len * 100).toFixed(2)}%">`
        + keys.filter((k) => r.parts[k] > 0).map((k) => `<i data-tip="${esc(tipFor(r, k))}" style="flex:${r.parts[k]} 1 0;background:var(${SRC_VAR[k]})"></i>`).join("")
        + "</div>"
      : "";
    let cap = "";
    if (opt.cap && r.allowed > 0) {
      const cp = sc.pos(conv(r.allowed));
      cap = cp <= 1
        ? `<span class="cap" style="left:calc(${(cp * 100).toFixed(2)}% - 1px)" data-tip="Installed capacity: ${esc(fmtG(r.allowed))}"></span>`
        : `<span class="cap off" data-tip="Installed capacity: ${esc(fmtG(r.allowed))}">capacity ${esc(fmtG(r.allowed))} ›</span>`;
    }
    let pill = "";
    if (opt.cap && r.allowed > 0 && r.total > 0) {
      const u = r.total / r.allowed;
      pill = r.overN > 0 ? `<span class="pill over">▲ ${r.overN} layer${r.overN > 1 ? "s" : ""} over capacity</span>`
        : u > 1 ? `<span class="pill over">▲ ${short(u)}× capacity</span>`
        : u >= 0.8 ? `<span class="pill near">◆ ${Math.round(u * 100)}% of capacity</span>`
        : `<span class="pill">${u < 0.01 ? "<1" : Math.round(u * 100)}% of capacity</span>`;
    }
    const inner = `<div class="lab"><b>${esc(r.label)}</b>${r.sub ? `<small>${esc(r.sub)}</small>` : ""}</div>`
      + `<div class="track">${grid}${bar}${cap}</div><div class="val">${esc(fmtV(v))}</div>${pill ? pill : "<span></span>"}`;
    h += r.clickable
      ? `<button class="crow" data-go="${esc(r.key)}">${inner}</button>`
      : `<div class="crow">${inner}</div>`;
  }
  h += '<div class="axis-row"><span></span><div class="axis">'
    + sc.ticks.map((t) => `<span style="left:${(sc.pos(t) * 100).toFixed(2)}%">${short(t / ad.d)}</span>`).join("")
    + `</div><span class="axis-cap">${esc(ad.u)}</span><span></span></div></div>`;
  return h;
}

/* ---------------------------------------------------------------- views */
function overviewView(res) {
  const rows = res.subs.map((r) => ({
    key: r.sub.id, label: r.sub.label, sub: r.sub.group, parts: r.parts, total: r.total, allowed: r.allowed, overN: r.over, clickable: true,
  })).sort((a, b) => b.total - a.total);
  let h = `<div class="panel-head"><h2>Readout by subdetector</h2><span class="spacer"></span>`
    + `<label class="muted" ${tipAttr("logscale")}><input type="checkbox" id="logscale" ${S.log ? "checked" : ""}> log scale</label></div>`;
  h += rows.length ? barChart(rows, { cap: true }) : '<div class="empty">No subdetector selected, or none has data for the chosen sample.</div>';
  h += '<p class="help muted" style="margin-top:8px">Select a row to see its layers.</p>';
  if (res.none.length) {
    h += `<div class="nodata">No rate available: ${res.none.map((r) => esc(r.sub.label)).join(", ")} (missing geometry or occupancy for the chosen sample).</div>`;
  }
  return h;
}

function drillView(res) {
  const subs = detData().detector.subdetectors.filter((s) => S.enabled.has(s.id));
  if (!subs.find((s) => s.id === S.sub) && subs.length) S.sub = subs[0].id;
  const rec = res.subs.concat(res.none).find((r) => r.sub.id === S.sub);
  let h = `<div class="panel-head"><h2>Layers</h2><div class="controls">`
    + `<select id="subsel" aria-label="Subdetector">${subs.map((s) => `<option value="${esc(s.id)}" ${s.id === S.sub ? "selected" : ""}>${esc(s.label)}</option>`).join("")}</select>`
    + `<label class="muted" ${tipAttr("logscale")}><input type="checkbox" id="logscale" ${S.log ? "checked" : ""}> log scale</label></div></div>`;
  if (!rec) return h + '<div class="empty">Select at least one subdetector in the sidebar.</div>';
  const rows = rec.rows.map((r) => ({
    key: String(r.layer.id), label: "Layer " + r.layer.id,
    sub: `${r.layer.n_cells.toLocaleString("en-US")} cells · occupancy ${pct(r.occ.ipc)}`,
    parts: r.parts, total: r.total, allowed: r.allowed,
  }));
  h += rows.length && rows.some((r) => r.total > 0) ? barChart(rows, { cap: true }) : '<div class="empty">No rate for this subdetector with the current selection.</div>';
  const notes = [];
  if (rec.sub.notes) notes.push(rec.sub.notes);
  for (const f of rec.sub.flags) if (FLAG_TEXT[f]) notes.push(FLAG_TEXT[f]);
  if (rec.noIpc) notes.push(`${rec.noIpc} layer(s) have no occupancy for the IPC sample ${esc(S.ipcSample)}; their IPC rate counts as zero.`);
  if (rec.rows.some((r) => r.over)) notes.push("A layer above 100% needs more links than installed: its readout does not fit the current link count.");
  h += `<div class="fibtbl"><div class="muted">Fibers per module, by layer (cells per module: ${S.model[rec.sub.id].cpm}; edit here or per subdetector in the sidebar; empty = default)</div><table class="model"><thead><tr><th>Layer</th><th>cells</th><th>fibers/module</th><th>fibers installed</th></tr></thead><tbody>`
    + rec.rows.map((r) => { const k = rec.sub.id + ":" + r.layer.id, m = S.model[rec.sub.id], ov = S.layerFib[k];
      return `<tr><td>${r.layer.id}</td><td class="num">${r.layer.n_cells.toLocaleString("en-US")}</td><td><input type="number" class="lfib" min="0" step="any" data-k="${esc(k)}" value="${ov == null ? "" : ov}" placeholder="${m.fpm ?? r.layer.fibers_per_module}" aria-label="Layer ${r.layer.id} fibers per module"></td><td class="num">${r.installed ? short(r.installed) : "–"}</td></tr>`; }).join("")
    + "</tbody></table></div>";
  h += notes.length ? `<ul class="notes">${notes.map((n) => `<li>${esc(n)}</li>`).join("")}</ul>` : "";
  return h;
}

function compareView() {
  const dd = detData();
  const ipcSamples = samplesOf("ipc");
  const cmpVar = ["--cmp-1", "--cmp-2", "--cmp-3"];
  const groups = [];
  for (const sub of dd.detector.subdetectors) {
    if (!S.enabled.has(sub.id)) continue;
    const m = S.model[sub.id];
    const series = ipcSamples.map((sm) => {
      let tot = 0, n = 0;
      for (const l of sub.layers) {
        if (l.n_cells == null) continue;
        const o = occ(sm, sub.id, l.id, S.stat);
        if (o != null) { tot += toGbps(o, S.rates.ipc, m, l.n_cells); n++; }
      }
      return { sm, tot, n, of: sub.layers.filter((l) => l.n_cells != null).length };
    });
    if (series.some((s) => s.n > 0)) groups.push({ sub, series });
  }
  let h = `<div class="panel-head"><h2>IPC samples side by side</h2><span class="spacer"></span>`
    + `<span class="muted">IPC only, ${S.stat === "p95" ? "95th percentile" : "mean"} occupancy, rate ${short(S.rates.ipc)} MHz</span></div>`;
  if (S.unit === "power" && !powerOk()) return h + `<div class="empty">Power needs a W/link value for ${esc(D.technologies[S.techId].label)}.</div>`;
  if (!groups.length) return h + '<div class="empty">No data.</div>';
  const all = groups.flatMap((g) => g.series.map((s) => conv(s.tot)));
  const sc = scaleFor(all);
  const ad = axisDiv(sc.ticks[sc.ticks.length - 1]);
  const grid = '<div class="grid-lines">' + sc.ticks.map((t) => `<u style="left:${(sc.pos(t) * 100).toFixed(2)}%"></u>`).join("") + "</div>";
  h += '<div class="legend">' + ipcSamples.map((s, i) => `<span><i class="sw" style="background:var(${cmpVar[i % 3]})"></i>${esc(s.label)}</span>`).join("") + "</div><div class=\"chart\">";
  for (const g of groups) {
    h += `<div class="group-head">${esc(g.sub.label)} <span class="muted">· ${esc(g.sub.group)}</span></div>`;
    g.series.forEach((s, i) => {
      const v = conv(s.tot);
      const bar = s.n > 0 && s.tot > 0
        ? `<div class="bar" style="width:${(sc.pos(v) * 100).toFixed(2)}%"><i style="flex:1;background:var(${cmpVar[i % 3]})" data-tip="${esc(s.sm.label + ": " + fmtG(s.tot))}"></i></div>` : "";
      const miss = s.n === 0 ? "no data" : s.n < s.of ? `${s.n}/${s.of} layers` : "";
      h += `<div class="crow"><div class="lab"><b>${esc(s.sm.id)}</b></div><div class="track">${grid}${bar}</div>`
        + `<div class="val">${s.n ? esc(fmtV(v)) : "–"}</div><span class="pill">${esc(miss)}</span></div>`;
    });
  }
  h += '<div class="axis-row"><span></span><div class="axis">' + sc.ticks.map((t) => `<span style="left:${(sc.pos(t) * 100).toFixed(2)}%">${short(t / ad.d)}</span>`).join("")
    + `</div><span class="axis-cap">${esc(ad.u)}</span><span></span></div></div>`;
  h += '<p class="help muted">Log scale applies here too. Layers without an occupancy value for a sample count as zero, so compare only rows without a "layers" note.</p>';
  return h;
}

const COLS = [
  ["sub", "Subdetector", (r) => r.sub.label, false],
  ["layer", "Layer", (r) => r.layer.id, false],
  ["cells", "Cells", (r) => r.layer.n_cells, true],
  ["channels", "Channels", (r) => r.channels, true],
  ["occ", "IPC occ.", (r) => r.occ.ipc, true],
  ["ipc", "IPC Gbps", (r) => r.parts.ipc, true],
  ["z", "Z Gbps", (r) => r.parts.z, true],
  ["inj", "Inj. Gbps", (r) => r.parts.inj, true],
  ["total", "Total Gbps", (r) => r.total, true],
  ["links", "Links req.", (r) => r.total / linkCap(), true],
  ["inst", "Links inst.", (r) => r.installed, true],
  ["util", "Utilisation", (r) => (r.allowed ? r.total / r.allowed : null), true],
];
function tableRows(res) {
  let rows = res.subs.flatMap((r) => r.rows);
  const f = S.filter.trim().toLowerCase();
  if (f) rows = rows.filter((r) => (r.sub.label + " " + r.layer.id).toLowerCase().indexOf(f) >= 0);
  const c = COLS.find((x) => x[0] === S.sort.col) || COLS[8];
  return rows.sort((a, b) => {
    const x = c[2](a), y = c[2](b);
    if (x == null) return 1; if (y == null) return -1;
    return (typeof x === "string" ? x.localeCompare(y) : x - y) * S.sort.dir;
  });
}
function cell(c, r) {
  const v = c[2](r);
  if (v == null) return '<td class="n dim">–</td>';
  switch (c[0]) {
    case "sub": case "layer": return `<td>${esc(v)}</td>`;
    case "cells": case "channels": return `<td class="n">${Number(v).toLocaleString("en-US")}</td>`;
    case "occ": return `<td class="n">${pct(v)}</td>`;
    case "util": return `<td class="n">${Math.round(v * 100)}%</td>`;
    default: return `<td class="n${v === 0 ? " dim" : ""}">${short(v)}</td>`;
  }
}
function tableView(res) {
  const rows = tableRows(res);
  let h = `<div class="panel-head"><h2>All layers</h2><div class="controls">`
    + `<input type="search" id="tfilter" placeholder="Filter subdetector or layer" value="${esc(S.filter)}" aria-label="Filter rows">`
    + `<button class="btn" id="copycsv">Copy CSV</button></div></div>`;
  h += '<div class="tablewrap"><table class="data"><thead><tr>'
    + COLS.map((c) => `<th data-col="${c[0]}" ${S.sort.col === c[0] ? `aria-sort="${S.sort.dir > 0 ? "ascending" : "descending"}"` : ""}>${esc(c[1])}</th>`).join("")
    + "</tr></thead><tbody>"
    + rows.map((r) => `<tr class="${r.over ? "over" : ""}">${COLS.map((c) => cell(c, r)).join("")}</tr>`).join("")
    + "</tbody></table></div>";
  h += `<p class="help muted">${rows.length} layers. Rates are in Gbps and links regardless of the unit chosen above. Rows on a red background need more links than installed.</p>`;
  if (S.csv) h += `<textarea class="csv" id="csvbox" readonly aria-label="CSV">${esc(S.csv)}</textarea>`;
  return h;
}
function toCsv(res) {
  const rows = tableRows(res);
  const head = COLS.map((c) => c[1]).join(",");
  const lines = rows.map((r) => COLS.map((c) => { const v = c[2](r); return v == null ? "" : typeof v === "number" ? +v.toPrecision(6) : '"' + String(v).replace(/"/g, '""') + '"'; }).join(","));
  return [head].concat(lines).join("\n");
}

/* ---------------------------------------------------------------- shell rendering */
function kpis(res) {
  const g = res.total, cap = linkCap();
  const basis = S.basis === "peak" ? "peak (link sizing)" : "average";
  const k = [];
  k.push(`<div class="kpi"><div class="k">Total readout</div><div class="v">${esc(fmtG(g))}</div><div class="d">${basis}, SF ${short(S.sf)}</div></div>`);
  k.push(`<div class="kpi"><div class="k">${esc(D.technologies[S.techId].label)}</div><div class="v">${short(res.linksCeil)} ${unitName()}</div><div class="d">rounded up per layer · ${short(res.links)} unrounded · ${short(cap)} Gbps each</div></div>`);
  k.push(`<div class="kpi"><div class="k">Installed links</div><div class="v">${short(res.installed)}</div><div class="d">${short((res.installed * cap) / 1000)} Tbps, silicon only (others undefined)</div></div>`);
  k.push(`<div class="kpi ${res.over ? "bad" : ""}"><div class="k">Layers over capacity</div><div class="v">${res.over} of ${res.layers}</div><div class="d">${res.noIpc ? res.noIpc + " layers without IPC occupancy" : "all layers have IPC occupancy"}</div></div>`);
  return k.join("");
}

function render() {
  const res = compute();
  S.res = res;
  $("#kpis").innerHTML = kpis(res);
  document.querySelectorAll("#tabs button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.tab === S.tab)));
  document.querySelectorAll("#units button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.unit === S.unit)));
  const view = { overview: overviewView, drill: drillView, compare: compareView, table: tableView }[S.tab];
  $("#view").innerHTML = view(res);
  const dd = detData();
  const smp = S.includeIpc ? sampleById(S.ipcSample) : null;
  $("#foot").innerHTML = `<div>${esc(D.working_points[D.settings.working_point].label)} · ${esc(dd.detector.label)} · ${esc(dd.detector.status)} · IPC ${S.includeIpc ? "sample " + esc(S.ipcSample || "–") : "excluded"}${S.includeZ && samplesOf("z").length ? " · Z " + (S.zSel === "mixture" ? "mixture" : esc(S.zSel)) : ""}${smp && smp.notes ? " · " + esc(smp.notes) : ""}</div>`;
}

/* sections are open unless closed; "model" starts closed */
function op(id, startClosed) {
  if (S.openSec[id] === undefined) S.openSec[id] = !startClosed;
  return S.openSec[id] ? "open" : "";
}

function sw(id, on, enabled, tip) {
  return `<label class="sw-toggle ${on ? "on" : "off"}" ${tip}><input type="checkbox" role="switch" id="${id}" ${on ? "checked" : ""} ${enabled ? "" : "disabled"}><span class="knob"></span><span class="st">${on ? "on" : "off"}</span></label>`;
}

function sidebarHtml() {
  const dd = detData(), t = techP(), tdef = D.technologies[S.techId];
  const hasZ = samplesOf("z").length > 0, hasInj = samplesOf("injection").length > 0;
  const ipcSel = samplesOf("ipc");
  const groups = {};
  for (const s of dd.detector.subdetectors) (groups[s.group] = groups[s.group] || []).push(s);
  const injSample = sampleById(S.injSample);
  const wps = D.working_points;
  const basisOn = S.includeInj && hasInj;
  const basisOff = "<b>Peak / average</b><br>Only the pulsed injection background differs between peak and average. Switch injection on to choose.";
  return `
  <details class="sec" data-sec="det" ${op("det")}><summary>Detector</summary>
    <div class="field stack" ${tipAttr("det")}><select id="det" aria-label="Detector">${D.index.detectors.map((d) => `<option value="${esc(d.id)}" ${d.id === S.det ? "selected" : ""}>${esc(d.label)} (${esc(d.status)})</option>`).join("")}</select></div>
  </details>
  <details class="sec" data-sec="wp" ${op("wp")}><summary>Working point</summary>
    <div class="field stack" ${tipAttr("wp")}><select id="wp" aria-label="Working point">${Object.entries(wps).map(([id, w]) => `<option value="${esc(id)}" ${id === D.settings.working_point ? "selected" : ""} ${id === D.settings.working_point ? "" : "disabled"}>${esc(w.label)}, ${w.sqrt_s_GeV} GeV${id === D.settings.working_point ? "" : " (no samples yet)"}</option>`).join("")}</select></div>
  </details>
  <details class="sec" data-sec="src" ${op("src")}><summary>Sources</summary>
    <div class="src ${S.includeIpc ? "" : "off"}"><div class="src-head"><i class="sw" style="background:var(--s-ipc)"></i><strong>IPC</strong><span class="tag">background · continuous</span>${sw("useipc", S.includeIpc, ipcSel.length > 0, tipAttr("useipc"))}</div>
      <div class="field stack" ${tipAttr("ipcsample")}><select id="ipcsample" aria-label="IPC sample"><option value="" ${S.includeIpc ? "" : "selected"}>No sample (excluded)</option>${ipcSel.map((s) => `<option value="${esc(s.id)}" ${s.id === S.ipcSample && S.includeIpc ? "selected" : ""}>${esc(s.label)}</option>`).join("")}</select></div>
      <div class="field" ${tipAttr("rateipc")}><label for="rate-ipc">Rate (MHz)</label><input type="number" id="rate-ipc" step="any" min="0" value="${S.rates.ipc}"></div></div>
    <div class="src ${S.includeZ && hasZ ? "" : "off"}"><div class="src-head"><i class="sw" style="background:var(--s-z)"></i><strong>Z decays</strong><span class="tag">physics · continuous</span>${sw("usez", S.includeZ && hasZ, hasZ, tipAttr("usez"))}</div>
      ${hasZ ? "" : '<div class="help">No Z sample for this detector.</div>'}
      <div class="field" ${tipAttr("ratez")}><label for="rate-z">Rate (MHz)</label><input type="number" id="rate-z" step="any" min="0" value="${S.rates.z}"></div>
      <div class="field stack" ${tipAttr("zsample")}><select id="zsample" aria-label="Z sample" ${hasZ ? "" : "disabled"}><option value="mixture" ${S.zSel === "mixture" ? "selected" : ""}>Mixture: ${D.sources.z.mixture.map((m) => m.weight + " × " + esc(m.sample)).join(" + ")}</option>${samplesOf("z").map((s) => `<option value="${esc(s.id)}" ${s.id === S.zSel ? "selected" : ""}>${esc(s.id)} only</option>`).join("")}</select></div>
      <div class="help">${S.zSel === "mixture" ? zParts().map((m) => esc(m.sample) + " (" + esc(sampleDesc(sampleById(m.sample))) + ")").join("<br>") : esc(sampleDesc(sampleById(S.zSel))) + ", counted as 100% of Z decays"}</div></div>
    <div class="src ${S.includeInj && hasInj ? "" : "off"}"><div class="src-head"><i class="sw" style="background:var(--s-inj)"></i><strong>Injection</strong><span class="tag">background · pulsed</span><span class="tag warn">example scale</span>${sw("useinj", S.includeInj && hasInj, hasInj, tipAttr("useinj"))}</div>
      ${hasInj ? "" : '<div class="help">No injection sample for this detector.</div>'}
      <div class="field stack" ${tipAttr("injsample")}><select id="injsample" aria-label="Injection sample" ${hasInj ? "" : "disabled"}>${samplesOf("injection").map((s) => `<option value="${esc(s.id)}" ${s.id === S.injSample ? "selected" : ""}>${esc(s.label)}</option>`).join("")}</select></div>
      <div class="help">${injSample ? esc(sampleDesc(injSample)) : ""}</div>
      <div class="field" ${tipAttr("injpeak")}><label for="injpeak">Peak scale (MHz-equiv.)</label><input type="number" id="injpeak" step="any" min="0" value="${S.injPeak}"></div>
      <div class="help">Average = peak × ${(duty() * 100).toFixed(1)}% duty (${D.sources.injection.pulse.turns_after_injection} turns per ${D.sources.injection.pulse.injection_period_s} s). The scale factors from the injection study are not in yet; this value is only an example.${injSample && injSample.flags.indexOf("no_stats") >= 0 ? " Sample has no statistics." : ""}</div></div>
    ${["sr", "beamgas", "touschek"].map((k) => { const c = (D.sources[k].components || []).map((x) => esc(x.label)).join(", ");
      return `<div class="src"><div class="src-head"><i class="sw" style="background:var(--line-strong)"></i><strong>${esc(D.sources[k].label)}</strong><span class="tag">background · continuous</span></div>
      <div class="help">No samples yet${c ? "; one per component (" + c + ", ...) will appear here" : ""}.</div></div>`; }).join("")}
  </details>
  <details class="sec" data-sec="basis" ${op("basis")}><summary>Occupancy and basis</summary>
    <div class="field stack"><div class="seg-ctl" id="stat"><button data-v="mean" aria-pressed="${S.stat === "mean"}" data-tip="${esc(TIPS.mean)}">Mean</button><button data-v="p95" aria-pressed="${S.stat === "p95"}" data-tip="${esc(TIPS.p95)}">95th percentile</button></div></div>
    <div class="field stack"><div class="seg-ctl" id="basis"><button data-v="peak" aria-pressed="${S.basis === "peak"}" ${basisOn ? "" : "disabled"} data-tip="${esc(basisOn ? TIPS.peak : basisOff)}">Peak</button><button data-v="avg" aria-pressed="${S.basis === "avg"}" ${basisOn ? "" : "disabled"} data-tip="${esc(basisOn ? TIPS.avg() : basisOff)}">Average</button></div></div>
    <div class="help">Peak adds the injection window to the continuous sources; average weights it by its duty.${basisOn ? "" : " Greyed out: injection is off, so both give the same result."}</div>
    <div class="field" ${tipAttr("sf")}><label for="sf">Safety factor</label><input type="number" id="sf" step="any" min="0" value="${S.sf}"></div>
  </details>
  <details class="sec" data-sec="tech" ${op("tech")}><summary>Link technology</summary>
    <div class="field stack" ${tipAttr("tech")}><select id="tech" aria-label="Link technology">${Object.entries(D.technologies).map(([id, x]) => `<option value="${esc(id)}" ${id === S.techId ? "selected" : ""}>${esc(x.label)}${x.verified ? "" : " (unverified)"}</option>`).join("")}</select></div>
    <div class="field" ${tipAttr("tline")}><label for="t-line">Line rate (Gbps)</label><input type="number" id="t-line" step="any" min="0" value="${t.line}"></div>
    <div class="field" ${tipAttr("teff")}><label for="t-eff">Encoding efficiency</label><input type="number" id="t-eff" step="any" min="0" max="1" value="${t.eff}"></div>
    <div class="field" ${tipAttr("tpow")}><label for="t-pow">Power per link (W)</label><input type="number" id="t-pow" step="any" min="0" value="${t.power == null ? "" : t.power}" placeholder="unknown"></div>
    <div class="help">Usable: <span class="num">${short(linkCap())}</span> Gbps per link. ${tdef.verified ? "" : "<span class=\"tag warn\">placeholder values</span> "}${esc(tdef.notes || "")}</div>
  </details>
  <details class="sec" data-sec="subs" ${op("subs")}><summary>Subdetectors</summary>
    <div class="mini"><button class="btn" id="all">All</button><button class="btn" id="none">None</button></div>
    ${Object.entries(groups).map(([g, subs]) => `<div class="group-title">${esc(g)}</div>` + subs.map((s) => `<div class="check"><input type="checkbox" class="subchk" id="sub-${esc(s.id)}" data-id="${esc(s.id)}" ${S.enabled.has(s.id) ? "checked" : ""}><label for="sub-${esc(s.id)}">${esc(s.label)}</label></div>`).join("")).join("")}
  </details>
  <details class="sec" data-sec="fib" ${op("fib", true)}><summary>Fibers per subdetector</summary>
    <table class="model"><thead><tr><th>Subdetector</th><th ${tipAttr("fpm")}>fib/module</th><th ${tipAttr("cpm")}>cells/module</th></tr></thead><tbody>
    ${dd.detector.subdetectors.map((s) => `<tr><td>${esc(s.label)}</td><td><input type="number" class="mdl" min="0" step="any" data-id="${esc(s.id)}" data-k="fpm" value="${S.model[s.id].fpm ?? ""}" placeholder="${s.layers[0].fibers_per_module}" aria-label="${esc(s.label)} fibers per module"></td><td><input type="number" class="mdl" min="0" step="any" data-id="${esc(s.id)}" data-k="cpm" value="${S.model[s.id].cpm}" aria-label="${esc(s.label)} cells per module"></td></tr>`).join("")}
    </tbody></table>
    <div class="help">Installed fibers per layer = fibers/module × cells ÷ cells/module. Empty fib/module keeps the layer values from the data. Setting any value for calorimeter, drift chamber or muon subdetectors enables their capacity.</div>
  </details>
  <details class="sec" data-sec="model" ${op("model", true)}><summary>Readout model per subdetector</summary>
    <table class="model"><thead><tr><th>Subdetector</th><th ${tipAttr("bits")}>bits/hit</th><th ${tipAttr("mult")}>mult.</th></tr></thead><tbody>
    ${dd.detector.subdetectors.map((s) => `<tr><td>${esc(s.label)}${s.flags.indexOf("multiplier_inherited_from_sheet") >= 0 ? ' <span class="tag warn" title="' + esc(FLAG_TEXT.multiplier_inherited_from_sheet) + '">inherited</span>' : ""}</td><td><input type="number" class="mdl" data-id="${esc(s.id)}" data-k="bits" value="${S.model[s.id].bits}" aria-label="${esc(s.label)} bits per hit"></td><td><input type="number" class="mdl" data-id="${esc(s.id)}" data-k="mult" step="any" value="${S.model[s.id].mult}" aria-label="${esc(s.label)} multiplier"></td></tr>`).join("")}
    </tbody></table></details>
  <div class="mini" style="margin-top:12px"><button class="btn" id="reset">Reset all inputs</button></div>`;
}

function renderSidebar() { $("#sidebar").innerHTML = sidebarHtml(); }

/* ---------------------------------------------------------------- events */
function num(el, fallback) { const v = parseFloat(el.value); return isFinite(v) ? v : fallback; }

function onSidebarInput(e) {
  const t = e.target;
  switch (t.id) {
    case "det": S = freshState(t.value, S); renderSidebar(); break;
    case "ipcsample":
      if (t.value === "") S.includeIpc = false; else { S.ipcSample = t.value; S.includeIpc = true; }
      renderSidebar(); break;
    case "useipc": S.includeIpc = t.checked; renderSidebar(); break;
    case "rate-ipc": S.rates.ipc = num(t, 0); break;
    case "zsample": S.zSel = t.value; break;
    case "injsample": S.injSample = t.value; break;
    case "rate-z": S.rates.z = num(t, 0); break;
    case "usez": S.includeZ = t.checked; renderSidebar(); break;
    case "useinj": S.includeInj = t.checked; renderSidebar(); break;
    case "injpeak": S.injPeak = num(t, 0); break;
    case "sf": S.sf = num(t, 1); break;
    case "tech": S.techId = t.value; renderSidebar(); break;
    case "t-line": techP().line = num(t, 0); renderLinkHelp(); break;
    case "t-eff": techP().eff = num(t, 1); renderLinkHelp(); break;
    case "t-pow": techP().power = t.value === "" ? null : num(t, null); break;
    default:
      if (t.classList.contains("subchk")) { t.checked ? S.enabled.add(t.dataset.id) : S.enabled.delete(t.dataset.id); }
      else if (t.classList.contains("mdl")) { S.model[t.dataset.id][t.dataset.k] = t.dataset.k === "fpm" && t.value === "" ? null : num(t, 0); }
      else return;
  }
  render();
}
function renderLinkHelp() {
  const el = document.querySelector("#sidebar #t-line");
  if (!el) return;
  const help = el.closest("section").querySelector(".help .num");
  if (help) help.textContent = short(linkCap());
}

function onSidebarClick(e) {
  const b = e.target.closest("button");
  if (!b) return;
  if (b.closest("#stat")) { S.stat = b.dataset.v; renderSidebar(); }
  else if (b.closest("#basis")) { S.basis = b.dataset.v; renderSidebar(); }
  else if (b.id === "all") { detData().detector.subdetectors.forEach((s) => S.enabled.add(s.id)); renderSidebar(); }
  else if (b.id === "none") { S.enabled.clear(); renderSidebar(); }
  else if (b.id === "reset") { S = freshState(S.det, S); renderSidebar(); }
  else return;
  render();
}

function onMainClick(e) {
  const go = e.target.closest("[data-go]");
  if (go) { S.sub = go.dataset.go; S.tab = "drill"; render(); return; }
  const th = e.target.closest("th[data-col]");
  if (th) {
    const c = th.dataset.col;
    S.sort = { col: c, dir: S.sort.col === c ? -S.sort.dir : (c === "sub" || c === "layer" ? 1 : -1) };
    render(); return;
  }
  if (e.target.id === "copycsv") {
    S.csv = toCsv(S.res);
    const done = () => { e.target.textContent = "Copied"; };
    const fallback = () => { render(); const box = $("#csvbox"); if (box) { box.focus(); box.select(); } };
    try { navigator.clipboard.writeText(S.csv).then(done, fallback); } catch (err) { fallback(); }
  }
}
function onMainInput(e) {
  if (e.target.classList.contains("lfib")) {
    const k = e.target.dataset.k, v = e.target.value, pos = e.target.selectionStart;
    if (v === "") delete S.layerFib[k]; else S.layerFib[k] = num(e.target, 0);
    render();
    const f = document.querySelector('#view .lfib[data-k="' + k + '"]'); if (f) { f.focus(); try { f.setSelectionRange(pos, pos); } catch (x) { /* number input */ } }
    return;
  }
  if (e.target.id === "subsel") { S.sub = e.target.value; render(); }
  else if (e.target.id === "logscale") { S.log = e.target.checked; render(); }
  else if (e.target.id === "tfilter") {
    S.filter = e.target.value; const pos = e.target.selectionStart; render();
    const f = $("#tfilter"); if (f) { f.focus(); f.setSelectionRange(pos, pos); }
  }
}

function setupTooltip() {
  const tip = $("#tip");
  document.addEventListener("mousemove", (e) => {
    const t = e.target.closest && e.target.closest("[data-tip]");
    if (!t) { tip.style.display = "none"; return; }
    tip.innerHTML = t.getAttribute("data-tip");
    tip.style.display = "block";
    const w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = Math.min(e.clientX + 14, window.innerWidth - w - 8) + "px";
    tip.style.top = (e.clientY + 16 + h > window.innerHeight ? e.clientY - h - 10 : e.clientY + 16) + "px";
  });
  document.addEventListener("mouseleave", () => { tip.style.display = "none"; });
}

function fillInfo() {
  const p = D.sources.injection.pulse;
  const v = {
    sf: D.settings.safety_factor, ipcrate: wpData().f_bc_MHz, zrate: D.sources.z.rate_MHz,
    zmix: D.sources.z.mixture.map((m) => m.weight + " × " + m.sample).join(" + "),
    turns: p.turns_after_injection, period: p.injection_period_s, duty: (duty() * 100).toFixed(1) + "%",
  };
  document.querySelectorAll("[data-fill]").forEach((e) => { e.textContent = v[e.dataset.fill]; });
}

/* theme: auto (follow the system) -> light -> dark; the choice is kept in the browser and shared by both pages */
function setupTheme(btn) {
  const KEY = "fcc-theme", modes = ["auto", "light", "dark"], icon = { auto: "\u25D0", light: "\u2600", dark: "\u263E" };
  let m = "auto";
  try { const v = localStorage.getItem(KEY); if (modes.includes(v)) m = v; } catch (e) { /* storage blocked */ }
  const apply = () => {
    if (m === "auto") document.documentElement.removeAttribute("data-theme"); else document.documentElement.setAttribute("data-theme", m);
    btn.textContent = icon[m];
    btn.title = "Theme: " + m + " (click to change)";
    btn.setAttribute("aria-label", btn.title);
  };
  apply();
  btn.addEventListener("click", () => {
    m = modes[(modes.indexOf(m) + 1) % modes.length];
    try { localStorage.setItem(KEY, m); } catch (e) { /* storage blocked */ }
    apply();
  });
}

function setupInfo() {
  const modal = $("#modal"), btn = $("#infobtn"), close = $("#modal-close");
  const open = () => { modal.hidden = false; $(".modal-panel").focus(); };
  const shut = () => { modal.hidden = true; btn.focus(); };
  btn.addEventListener("click", open);
  close.addEventListener("click", shut);
  modal.addEventListener("click", (e) => { if (e.target === modal) shut(); });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && !modal.hidden) shut(); });
  fillInfo();
}

function storeGet(k) { try { return localStorage.getItem(k); } catch (e) { return null; } }
function storeSet(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* ignore */ } }

async function main() {
  try { D = await loadData(); } catch (err) {
    $("#view").innerHTML = `<div class="empty">Could not load the data (${esc(err.message)}). Open the single-file build or serve this folder over http.</div>`;
    return;
  }
  S = freshState(D.index.detectors[0].id);
  renderSidebar();
  $("#sidebar").addEventListener("input", onSidebarInput);
  $("#sidebar").addEventListener("click", onSidebarClick);
  $("#sidebar").addEventListener("toggle", (e) => { if (e.target.dataset && e.target.dataset.sec) S.openSec[e.target.dataset.sec] = e.target.open; }, true);
  $("#view").addEventListener("click", onMainClick);
  $("#view").addEventListener("input", onMainInput);
  $("#tabs").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { S.tab = b.dataset.tab; render(); } });
  $("#units").addEventListener("click", (e) => { const b = e.target.closest("button"); if (b) { S.unit = b.dataset.unit; render(); } });
  const lay = $("#layout"), tg = $("#toggle");
  const setCollapsed = (c) => { lay.classList.toggle("collapsed", c); tg.setAttribute("aria-expanded", String(!c)); tg.textContent = c ? "Show controls" : "Hide controls"; storeSet("tdaq-collapsed", c ? "1" : "0"); };
  tg.addEventListener("click", () => setCollapsed(!lay.classList.contains("collapsed")));
  setCollapsed(storeGet("tdaq-collapsed") === "1");
  setupTooltip();
  setupInfo();
  setupTheme(document.querySelector("#themebtn"));
  // penguin disabled (see page.html)
  // if (D.facts && D.facts.length) $("#penguin-bubble").textContent = D.facts[Math.floor(Math.random() * D.facts.length)];
  render();
}
main();
