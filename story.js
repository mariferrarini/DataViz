// One story, told once: a network sorts itself into modules (arcs), the same network
// is shown as an adjacency matrix that sorts into blocks, then a zoom into one module
// reveals its members' sequences, which align around a shared motif. Then it all stops.
(function () {
  // Colours come from styles.css, so the palette lives in one place
  const css = getComputedStyle(document.documentElement);
  const INK = css.getPropertyValue("--viz-ink-rgb").trim(); // as "r, g, b"
  const INK2 = css.getPropertyValue("--viz-ink-2-rgb").trim(); // second layer
  const INK3 = css.getPropertyValue("--viz-ink-3-rgb").trim(); // the protein's second motif
  const AXIS = css.getPropertyValue("--border").trim();

  // ---------- Seeded randomness ----------
  function rng(seed) {
    return function () {
      seed |= 0;
      seed = (seed + 0x6d2b79f5) | 0;
      let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function shuffle(a, r = Math.random) {
    for (let i = a.length - 1; i > 0; i--) {
      const j = Math.floor(r() * (i + 1));
      [a[i], a[j]] = [a[j], a[i]];
    }
    return a;
  }

  // ---------- Network: two layers of connections over the same 8 planted modules ----------
  // Layer A (DNA/RNA: correlated sequence and expression signals) is dense; layer B (protein:
  // shared domains) is sparse and almost always within a module. A is drawn on one side of
  // the arcs and in the matrix's upper triangle, B on the other side and in the lower triangle.
  const rand = rng(1);
  const SIZES = [12, 9, 8, 7, 6, 5, 5, 4];
  const N = SIZES.reduce((a, b) => a + b, 0);
  const planted = SIZES.flatMap((s, g) => Array(s).fill(g));
  function layer(pIn, pOut) {
    const l = [];
    for (let i = 0; i < N; i++)
      for (let j = i + 1; j < N; j++) if (rand() < (planted[i] === planted[j] ? pIn : pOut)) l.push([i, j]);
    return l;
  }
  const layerA = layer(0.65, 0.02);
  const layerB = layer(0.3, 0.003);
  // Communities are detected on both layers together
  const seen = new Set();
  const edges = [...layerA, ...layerB].filter(([i, j]) => !seen.has(i * N + j) && seen.add(i * N + j));
  const neighbors = Array.from({ length: N }, () => []);
  for (const [i, j] of edges) {
    neighbors[i].push(j);
    neighbors[j].push(i);
  }
  const degree = neighbors.map((n) => n.length);

  // ---------- Community detection: Louvain (doesn't see the planted modules) ----------
  // Moves nodes between communities while modularity improves, then merges each
  // community into one weighted node and repeats until nothing moves.
  function louvain(r) {
    const m2 = 2 * edges.length;
    let n = N;
    let adj = neighbors.map((nb) => new Map(nb.map((u) => [u, 1])));
    let member = [...Array(N).keys()]; // original node → current merged node
    for (;;) {
      const k = adj.map((a, v) => [...a].reduce((s, [u, w]) => s + (u === v ? 2 * w : w), 0));
      const com = [...Array(n).keys()];
      const tot = k.slice();
      let moved = false;
      for (let improved = true; improved; ) {
        improved = false;
        for (const v of shuffle([...Array(n).keys()], r)) {
          const cv = com[v];
          const links = new Map();
          for (const [u, w] of adj[v]) if (u !== v) links.set(com[u], (links.get(com[u]) || 0) + w);
          tot[cv] -= k[v];
          let best = cv;
          let bestGain = (links.get(cv) || 0) - (tot[cv] * k[v]) / m2;
          for (const [c, w] of links) {
            const gain = w - (tot[c] * k[v]) / m2;
            if (gain > bestGain + 1e-12) [best, bestGain] = [c, gain];
          }
          tot[best] += k[v];
          if (best !== cv) {
            com[v] = best;
            improved = moved = true;
          }
        }
      }
      if (!moved) break;
      const ids = new Map();
      for (const c of com) if (!ids.has(c)) ids.set(c, ids.size);
      const next = Array.from({ length: ids.size }, () => new Map());
      for (let v = 0; v < n; v++)
        for (const [u, w] of adj[v]) {
          const a = ids.get(com[v]);
          const b = ids.get(com[u]);
          const add = u === v ? w : w / 2; // each edge is visited from both ends
          next[a].set(b, (next[a].get(b) || 0) + add);
          if (a !== b) next[b].set(a, (next[b].get(a) || 0) + add);
        }
      member = member.map((s) => ids.get(com[s]));
      adj = next;
      n = ids.size;
    }
    return member;
  }
  function modularity(com) {
    const m = edges.length;
    let q = edges.filter(([i, j]) => com[i] === com[j]).length / m;
    const tot = new Map();
    for (let v = 0; v < N; v++) tot.set(com[v], (tot.get(com[v]) || 0) + degree[v]);
    for (const t of tot.values()) q -= (t / (2 * m)) ** 2;
    return q;
  }
  function detectCommunities(r, runs = 10) {
    let best = null;
    let bestQ = -Infinity;
    for (let i = 0; i < runs; i++) {
      const c = louvain(r);
      const q = modularity(c);
      if (!best || q > bestQ) [best, bestQ] = [c, q];
    }
    // Renumber communities 0..k-1, largest first
    const sizes = new Map();
    for (const c of best) sizes.set(c, (sizes.get(c) || 0) + 1);
    const rank = new Map([...sizes].sort((a, b) => b[1] - a[1]).map(([c], i) => [c, i]));
    return best.map((c) => rank.get(c));
  }
  const community = detectCommunities(rng(11));
  // Smallest community first, so the largest block ends at the matrix's bottom-right
  // corner, right next to the zoom panel
  const grouped = [...Array(N).keys()].sort(
    (a, b) => community[b] - community[a] || degree[b] - degree[a] || a - b,
  );
  const zoomMembers = grouped.filter((v) => community[v] === 0);

  // Position of each node (in slots) for a given order
  function positions(order) {
    const p = new Array(N);
    order.forEach((v, r) => (p[v] = r));
    return p;
  }

  // ---------- The zoomed module's genes: DNA and protein, each hiding motifs ----------
  // One row per gene in each track. The DNA layer (navy) shares a near-exact motif, AGGAGG
  // (Shine-Dalgarno, the ribosome binding site). The protein layer is insect defensin-like
  // (the antimicrobial peptides insects use to keep bacteria, symbionts included, in check):
  // C-x-C on the last beta-strand (blue-grey) and C-x-x-x-C on the helix (blue), far apart in
  // sequence (the ellipsis stands for the residues between them) but bonded together in 3D.
  const seqRand = rng(5);
  const pick = (alphabet) => alphabet[Math.floor(seqRand() * alphabet.length)];
  const NT = "ACGT";
  const AA = "ADEFGHIKLMNPQRSTVWY"; // no C, so the motifs' cysteines stand out
  // Each row: random letters with a motif block written at a random offset. `motif()` returns
  // the block's letters and, per letter, which motif it belongs to (-1: not fixed).
  function makeRows(len, blockLen, alphabet, motif) {
    return zoomMembers.map(() => {
      const at = Math.floor(seqRand() * (len - blockLen + 1));
      const s = Array.from({ length: len }, () => pick(alphabet));
      const { letters, motifOf } = motif();
      letters.forEach((c, k) => (s[at + k] = c));
      return { s, at, motifOf: s.map((_, k) => (k >= at && k < at + blockLen ? motifOf[k - at] : -1)) };
    });
  }
  const dnaRows = makeRows(14, 6, NT, () => {
    const letters = "AGGAGG".split("");
    if (seqRand() < 0.35) letters[Math.floor(seqRand() * 6)] = pick(NT); // the odd mismatch
    return { letters, motifOf: letters.map((c, k) => (c === "AGGAGG"[k] ? 0 : -1)) };
  });
  // Block "C x x x C … C x C": motif 1 is C-x-x-x-C, motif 0 is C-x-C
  const PROTEIN_MOTIF_OF = [1, -1, -1, -1, 1, -1, 0, -1, 0];
  const proteinRows = makeRows(15, 9, AA, () => ({
    letters: ["C", pick(AA), pick(AA), pick(AA), "C", "…", "C", pick(AA), "C"],
    motifOf: PROTEIN_MOTIF_OF,
  }));
  // The first protein row is real insect defensin A (PDB 1ICA): C16-A-A-H-C20 … C36-V-C38
  proteinRows[0] = { s: "NHSACAAHC…CVCRN".split(""), at: 4, motifOf: [-1, -1, -1, -1, ...PROTEIN_MOTIF_OF, -1, -1] };
  // Motif positions within each track's block, and their colours
  const tracks = [
    { rows: dnaRows, len: 14, rgb: INK, rest: 0.3, motifs: [{ from: 0, len: 6, rgb: INK }] },
    {
      rows: proteinRows,
      len: 15,
      rgb: INK2,
      rest: 0.45,
      motifs: [
        { from: 6, len: 3, rgb: INK2 }, // C-x-C, revealed with the alignment
        { from: 0, len: 5, rgb: INK3 }, // C-x-x-x-C, revealed with the 3D structure
      ],
    },
  ];
  for (const tr of tracks) tr.maxAt = Math.max(...tr.rows.map((r) => r.at));

  // ---------- Layout ----------
  const stage = document.querySelector(".stage");
  const canvas = document.getElementById("stage-canvas");
  const ctx = canvas.getContext("2d");
  let W, H, wide;
  const P = {};

  function layout() {
    W = stage.clientWidth;
    wide = W >= 880;
    if (wide) {
      // arcs | matrix | sequences | structure, every node level with its matrix row
      const ms = Math.min(400, (W - 84) / 2.85);
      const arcW = 0.5 * ms;
      const seqW = 0.6 * ms;
      const structW = 0.75 * ms;
      const left = (W - (arcW + 20 + ms + 32 + seqW + 32 + structW)) / 2; // centre the composition
      P.axis = left + arcW / 2;
      P.matrix = { x: left + arcW + 20, y: 4, w: ms, h: ms };
      P.seq = { x: P.matrix.x + ms + 32, y: 4, w: seqW, h: ms };
      P.struct = { x: P.seq.x + seqW + 32, y: 4, w: structW, h: ms };
      P.arcFlat = (arcW / 2 - 6) / (ms / 2);
      H = ms + 8;
    } else {
      // Stacked: arcs above and below a horizontal axis, every node over its matrix column
      const ms = Math.min(W, 340);
      const arcH = Math.min(70, ms * 0.22);
      const sz = Math.min(W, 260);
      P.axis = arcH + 12;
      P.matrix = { x: (W - ms) / 2, y: P.axis + arcH + 18, w: ms, h: ms };
      P.seq = { x: 0, y: P.matrix.y + ms + 28, w: W, h: 450 };
      P.struct = { x: (W - sz) / 2, y: P.seq.y + P.seq.h + 28, w: sz, h: sz * 1.3 };
      P.arcFlat = arcH / (ms / 2);
      H = P.struct.y + P.struct.h + 4;
    }
    stage.style.height = `${H}px`;
    Object.assign(structEl.style, { left: `${P.struct.x}px`, top: `${P.struct.y}px`, width: `${P.struct.w}px`, height: `${P.struct.h}px` });
    if (!viewer && !viewerFailed) createViewer();
    else if (viewer) {
      viewer.resize();
      viewer.zoomTo();
      viewer.zoom(0.92);
      viewer.render();
    }
    const dpr = window.devicePixelRatio || 1;
    canvas.width = W * dpr;
    canvas.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  // ---------- Drawing ----------
  const lerp = (a, b, t) => a + (b - a) * t;

  // Arcs and matrix share one order, so every node keeps its matrix row/column
  const shuffled = shuffle([...Array(N).keys()]);
  const from = positions(shuffled);
  const to = positions(grouped);

  // Two-sided arcs: layer A on the left (wide) or top (narrow), layer B on the other side.
  // Each side has its own column of dots, a few pixels either side of the axis.
  function drawArcs(e) {
    const m = P.matrix;
    const cell = m.w / N;
    const at = (v) => (wide ? m.y : m.x) + (lerp(from[v], to[v], e) + 0.5) * cell;
    const sides = [
      { edges: layerA, rgb: INK, alpha: 0.3, off: -4 },
      { edges: layerB, rgb: INK2, alpha: 0.45, off: 4 },
    ];
    const dot = Math.min(2.2, cell * 0.38);
    for (const { edges: list, rgb, alpha, off } of sides) {
      const base = P.axis + off;
      const out = Math.sign(off); // -1: left/up, +1: right/down
      ctx.lineWidth = 1;
      ctx.strokeStyle = `rgba(${rgb},${alpha})`;
      ctx.beginPath();
      for (const [i, j] of list) {
        const a = Math.min(at(i), at(j));
        const r = (Math.max(at(i), at(j)) - a) / 2;
        const f = r * P.arcFlat;
        if (wide) {
          // Half-ellipse from node a to node b, bulging left or right of the axis
          ctx.moveTo(base, out < 0 ? a + 2 * r : a);
          ctx.ellipse(base, a + r, f, r, 0, out < 0 ? Math.PI / 2 : -Math.PI / 2, out < 0 ? (3 * Math.PI) / 2 : Math.PI / 2);
        } else {
          ctx.moveTo(out < 0 ? a : a + 2 * r, base);
          ctx.ellipse(a + r, base, r, f, 0, out < 0 ? Math.PI : 0, out < 0 ? 2 * Math.PI : Math.PI);
        }
      }
      ctx.stroke();
      ctx.fillStyle = `rgb(${rgb})`;
      for (let v = 0; v < N; v++) {
        ctx.beginPath();
        if (wide) ctx.arc(base, at(v), dot, 0, 2 * Math.PI);
        else ctx.arc(at(v), base, dot, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
  }

  function drawMatrix(r, e, alpha) {
    const cell = r.w / N;
    const pos = from.map((p, v) => lerp(p, to[v], e));
    ctx.globalAlpha = alpha;
    const put = (row, col) => ctx.fillRect(r.x + col * cell, r.y + row * cell, cell, cell);
    ctx.fillStyle = `rgba(${INK},0.85)`;
    for (let v = 0; v < N; v++) put(pos[v], pos[v]);
    // Layer A above the diagonal, layer B below it
    for (const [i, j] of layerA) put(Math.min(pos[i], pos[j]), Math.max(pos[i], pos[j]));
    ctx.fillStyle = `rgba(${INK2},0.85)`;
    for (const [i, j] of layerB) put(Math.max(pos[i], pos[j]), Math.min(pos[i], pos[j]));
    ctx.strokeStyle = AXIS;
    ctx.lineWidth = 1;
    ctx.strokeRect(r.x + 0.5, r.y + 0.5, r.w - 1, r.h - 1);
    ctx.globalAlpha = 1;
  }

  // Where the sequence text sits: DNA rows on top, a blank row, then protein rows for the
  // same genes. On wide screens it spans the matrix's full height, a true zoom of the block.
  function seqBox(r) {
    const cols = Math.max(...tracks.map((tr) => tr.maxAt + tr.len));
    const cw = Math.min(11, (r.w - 8) / cols);
    const n = zoomMembers.length;
    const rowH = wide ? r.h / (2 * n + 1) : Math.min(18, r.h / (2 * n + 1));
    const x = wide ? r.x + 8 : r.x + (r.w - cols * cw) / 2;
    const y = wide ? r.y : r.y + (r.h - rowH * (2 * n + 1)) / 2;
    return { cw, rowH, x, y, w: cols * cw, h: rowH * (2 * n + 1), trackY: [y, y + (n + 1) * rowH], trackH: n * rowH };
  }

  // Frame the largest module's block, whose genes the sequence panel shows
  function drawZoom(z) {
    const m = P.matrix;
    const size = (zoomMembers.length / N) * m.w;
    const b = { x: m.x + m.w - size, y: m.y + m.h - size, w: size, h: size };
    ctx.strokeStyle = `rgba(${INK},${(0.8 * z).toFixed(3)})`;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(b.x - 1.5, b.y - 1.5, b.w + 3, b.h + 3);
  }

  // Both tracks slide into alignment (e); their first motif darkens as they align, the
  // protein's second motif only when the 3D structure shows it (reveal2)
  function drawSequences(r, e, alpha, reveal2) {
    const { cw, rowH, x: x0, trackY, trackH } = seqBox(r);
    ctx.globalAlpha = alpha;
    ctx.font = `${Math.min(cw * 1.55, rowH * 0.8)}px ui-monospace, Menlo, monospace`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    tracks.forEach((tr, k) => {
      const y = trackY[k];
      const reveal = tr.motifs.map((_, m) => (m === 0 ? e : reveal2));
      tr.motifs.forEach((m, i) => {
        ctx.fillStyle = `rgba(${m.rgb},${(0.08 * reveal[i] * reveal[i]).toFixed(3)})`;
        ctx.fillRect(x0 + (tr.maxAt + m.from) * cw - 2, y, m.len * cw + 4, trackH);
      });
      tr.rows.forEach((row, i) => {
        const shift = (tr.maxAt - row.at) * e;
        row.s.forEach((c, j) => {
          const m = row.motifOf[j];
          const a = m >= 0 ? tr.rest + (1 - tr.rest) * reveal[m] : tr.rest;
          const rgb = m >= 0 && reveal[m] > 0 ? tr.motifs[m].rgb : tr.rgb;
          ctx.fillStyle = `rgba(${rgb},${a.toFixed(3)})`;
          ctx.fillText(c, x0 + (j + shift + 0.5) * cw, y + (i + 0.5) * rowH);
        });
      });
    });
    ctx.textBaseline = "alphabetic";
    ctx.globalAlpha = 1;
  }

  // Frame the protein motif column, the active site of the enzyme shown next
  function drawZoom2(z) {
    const t = seqBox(P.seq);
    const tr = tracks[1];
    const m = tr.motifs[0];
    const b = { x: t.x + (tr.maxAt + m.from) * t.cw - 2, y: t.trackY[1] - 2, w: m.len * t.cw + 4, h: t.trackH + 4 };
    ctx.strokeStyle = `rgba(${INK},${(0.8 * z).toFixed(3)})`;
    ctx.lineWidth = 1.2;
    ctx.strokeRect(b.x, b.y, b.w, b.h);
  }

  // Insect defensin A (PDB 1ICA) in 3D with 3Dmol.js: a thin, pale cartoon of its helix and
  // two-stranded sheet. The strand's C-x-C is in the same blue-grey as in the sequences; the
  // helix's C-x-x-x-C, bonded to it by two disulfides, colours in blue as the phrase types.
  const PROTEIN = window.DEFENSIN_1ICA;
  const hex = (rgb) => "#" + rgb.split(",").map((v) => (+v).toString(16).padStart(2, "0")).join("");
  const PALE = "220, 228, 238";
  const structEl = document.createElement("div");
  structEl.className = "structure";
  stage.insertBefore(structEl, canvas); // under the canvas, so the zoom frames draw on top
  let viewer = null;
  let viewerFailed = false;
  let shownReveal = -1;
  // Created on the first layout, once its container has a size
  function createViewer() {
    try {
      // A solid background in the page colour: with a transparent canvas, see-through
      // cartoon pixels blend with the page and wash out to grey
      const bg = css.getPropertyValue("--bg").trim();
      viewer = $3Dmol.createViewer(structEl, { backgroundColor: bg, backgroundAlpha: 1, antialias: true });
      viewer.addModel(PROTEIN.pdb, "pdb");
      viewer.enableFog(false);
      styleStructure(0);
      viewer.zoomTo();
      viewer.zoom(0.92);
      viewer.render();
    } catch (err) {
      viewer = null; // no WebGL: the rest of the story still plays
      viewerFailed = true;
    }
  }

  // reveal 0 → 1: the helix motif goes from the pale cartoon colour to blue
  function styleStructure(reveal) {
    const mix = PALE.split(",").map((v, i) => Math.round(+v + (+INK3.split(",")[i] - +v) * reveal)).join(",");
    viewer.setStyle({}, { cartoon: { color: "#c3cfdf", thickness: 0.2, arrows: true } });
    viewer.setStyle({ resi: PROTEIN.site }, { cartoon: { color: hex(INK2) } });
    viewer.setStyle({ resi: PROTEIN.site2 }, { cartoon: { color: hex(mix) } });
  }

  function drawStructure(alpha, reveal) {
    if (!viewer) return;
    structEl.style.opacity = alpha;
    const step = Math.round(reveal * 20) / 20; // restyle in 5% steps, not every frame
    if (step !== shownReveal) {
      styleStructure(step);
      viewer.render();
      shownReveal = step;
    }
  }

  // ---------- Timeline: text and visuals are both a function of time ----------
  const PHRASES = [
    "I like to find patterns in complex datasets.",
    "I like to play with data visualization.",
    "I'm passionate about sequences and the motifs they hide.",
    "I love seeing how structure shapes function.",
  ];
  document.getElementById("tagline-text").textContent = PHRASES.join(" ");
  const TYPE_MS = 55;
  const DELETE_MS = 28;
  const HOLD_MS = 2200;
  // Keyframes: type 1 · hold · delete 1 · type 2 · hold · delete 2 · type 3 · hold · delete 3 ·
  // type 4 · hold · delete 4
  const T = [0];
  const add = (ms) => T.push(T[T.length - 1] + ms);
  add(PHRASES[0].length * TYPE_MS);
  add(HOLD_MS);
  add(PHRASES[0].length * DELETE_MS);
  add(PHRASES[1].length * TYPE_MS);
  add(HOLD_MS);
  add(PHRASES[1].length * DELETE_MS);
  add(PHRASES[2].length * TYPE_MS);
  add(HOLD_MS);
  add(PHRASES[2].length * DELETE_MS);
  add(PHRASES[3].length * TYPE_MS);
  add(HOLD_MS);
  add(PHRASES[3].length * DELETE_MS);
  const END = T[12];

  const ease = (t) => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);
  // Eased progress through keyframes a → b, optionally only over fraction f0..f1 of it
  function seg(t, a, b, f0 = 0, f1 = 1) {
    const u = (t - T[a]) / (T[b] - T[a]);
    return ease(Math.max(0, Math.min(1, (u - f0) / (f1 - f0))));
  }

  function textAt(t) {
    const typed = (p, a) => PHRASES[p].slice(0, Math.floor((t - T[a]) / TYPE_MS));
    const deleted = (p, a) => PHRASES[p].slice(0, Math.max(0, PHRASES[p].length - Math.floor((t - T[a]) / DELETE_MS)));
    if (t < T[1]) return typed(0, 0);
    if (t < T[2]) return PHRASES[0];
    if (t < T[3]) return deleted(0, 2);
    if (t < T[4]) return typed(1, 3);
    if (t < T[5]) return PHRASES[1];
    if (t < T[6]) return deleted(1, 5);
    if (t < T[7]) return typed(2, 6);
    if (t < T[8]) return PHRASES[2];
    if (t < T[9]) return deleted(2, 8);
    if (t < T[10]) return typed(3, 9);
    if (t < T[11]) return PHRASES[3];
    if (t < T[12]) return deleted(3, 11);
    return "";
  }

  function draw(t) {
    ctx.clearRect(0, 0, W, H);
    // Arcs sort (1) and scramble again (2→3); the matrix then fades in unordered
    // and sorts itself, with the arcs re-sorting in step so rows always match (3→4)
    const matSort = seg(t, 3, 4, 0.25, 1);
    const arcSort = t < T[2] ? seg(t, 0, 1) : t < T[3] ? 1 - seg(t, 2, 3) : matSort;
    drawArcs(arcSort);
    const matAlpha = seg(t, 3, 4, 0, 0.25);
    if (matAlpha > 0) drawMatrix(P.matrix, matSort, matAlpha);
    const zoom = seg(t, 5, 6);
    if (zoom > 0) {
      drawZoom(zoom);
      drawSequences(P.seq, seg(t, 6, 7), Math.max(0, (zoom - 0.4) / 0.6), seg(t, 9, 10));
    }
    // Zoom from the protein motif into the protein (8→9); its second motif then colours in,
    // in 3D and in the sequences at once (9→10)
    const zoom2 = seg(t, 8, 9);
    if (zoom2 > 0) drawZoom2(zoom2);
    drawStructure(Math.max(0, (zoom2 - 0.4) / 0.6), seg(t, 9, 10));
  }

  const typer = document.getElementById("typewriter");
  const cursor = document.querySelector(".cursor");
  let t0 = null;
  let now = 0;
  let shown = null;

  function frame(ts) {
    if (t0 === null) t0 = ts;
    now = Math.min(END, ts - t0);
    const text = textAt(now);
    if (text !== shown) typer.textContent = shown = text;
    draw(now);
    if (now < END) requestAnimationFrame(frame);
    else setTimeout(() => cursor.classList.add("done"), 1500); // and then it all stops
  }

  layout();
  let resizeId;
  window.addEventListener("resize", () => {
    clearTimeout(resizeId);
    resizeId = setTimeout(() => {
      layout();
      draw(now);
    }, 100);
  });

  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    // No typing or motion: show the finished picture
    now = END;
    cursor.classList.add("done");
    draw(END);
  } else {
    draw(0);
    setTimeout(() => requestAnimationFrame(frame), 800);
  }
})();
