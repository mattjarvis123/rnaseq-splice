// Gene lookup page of the rnaseq-splice report: draws expression, fold change and splicing plots for
// one gene from window.GL_DATA (written in front of this script by child_gene_lookup.Rmd). Uses the Plotly library already loaded by the report.
(function () {
  function start() {
    var D = window.GL_DATA;
    if (!D || !document.getElementById("gl-input") || typeof Plotly === "undefined") return;
    var input = document.getElementById("gl-input");
    var summary = document.getElementById("gl-summary");
    var current = null;

    // Look-ups by symbol and by Ensembl id (case-insensitive)
    var index = new Map();
    D.gene_name.forEach(function (n, i) { var k = String(n).toUpperCase(); if (!index.has(k)) index.set(k, i); });
    D.gene_id.forEach(function (g, i) { index.set(String(g).toUpperCase(), i); });
    var evByGene = new Map();
    if (D.ev) D.ev.gene.forEach(function (g, j) { if (!evByGene.has(g)) evByGene.set(g, []); evByGene.get(g).push(j); });

    var dl = document.getElementById("gl-genes");
    var seen = new Set(), frag = document.createDocumentFragment();
    D.gene_name.forEach(function (n) {
      if (seen.has(n)) return; seen.add(n);
      var o = document.createElement("option"); o.value = n; frag.appendChild(o);
    });
    dl.appendChild(frag);

    var fmtP = function (p) { return p === null ? "NA" : Number(p).toExponential(2); };
    var axis = { showline: true, mirror: true, linecolor: "#333", zeroline: false, gridcolor: "#EBEBEB" };
    var cfg = { responsive: true, displaylogo: false };
    function message(id, text) {
      Plotly.purge(id);
      document.getElementById(id).innerHTML = '<div class="gl-msg">' + text + "</div>";
    }
    function clear(id) { var d = document.getElementById(id); if (d.querySelector(".gl-msg")) d.innerHTML = ""; }

    // i indexes all genes (search list); r = D.row[i] indexes the plot data, null = not significant
    function draw(gi) {
      current = gi;
      var name = D.gene_name[gi], id = D.gene_id[gi], i = D.row[gi];
      input.value = name;
      if (i === null) {
        summary.innerHTML = "<b>" + name + "</b> " + id + " &nbsp;·&nbsp; not significant in any contrast";
        var ns = "Non-significant: " + name + " is not differentially expressed (" + D.de_rule + ") or differentially spliced (" +
                 D.ds_rule + ") in any contrast. Its values are in the full limma and rMATS result files.";
        ["gl-expr", "gl-lfc", "gl-psi"].forEach(function (b) { message(b, ns); });
        document.getElementById("gl-ev-table").innerHTML = "";
        return;
      }
      // summary line
      var parts = [];
      D.contrasts.forEach(function (ct, k) {
        var l = D.lfc[i] ? D.lfc[i][k] : null, f = D.fdr[i] ? D.fdr[i][k] : null;
        if (l !== null) parts.push(ct + ": log2FC " + l.toFixed(2) + " (FDR " + fmtP(f) + ")");
      });
      summary.innerHTML = "<b>" + name + "</b> " + id + (parts.length ? " &nbsp;·&nbsp; " + parts.join(" &nbsp;·&nbsp; ") : " &nbsp;·&nbsp; not tested for expression");

      // 1. expression per sample, grouped by condition
      if (D.expr && D.expr[i] && D.expr[i].some(function (v) { return v !== null; })) {
        clear("gl-expr");
        var tr = D.cond_order.map(function (cd) {
          var idx = []; D.conds.forEach(function (c, s) { if (c === cd) idx.push(s); });
          return { type: "box", name: cd, x: idx.map(function () { return cd; }), y: idx.map(function (s) { return D.expr[i][s]; }),
                   text: idx.map(function (s) { return D.samples[s]; }), boxpoints: "all", jitter: 0.3, pointpos: 0,
                   marker: { color: D.cond_cols[cd], size: 8, line: { width: 1, color: "#333" } }, line: { color: D.cond_cols[cd] },
                   fillcolor: "rgba(0,0,0,0)", hovertemplate: "%{text}<br>log-CPM %{y:.2f}<extra></extra>" };
        });
        Plotly.react("gl-expr", tr, { showlegend: false, margin: { t: 10, b: 90 }, plot_bgcolor: "white",
          xaxis: Object.assign({ title: "", tickangle: -30, type: "category", categoryorder: "array", categoryarray: D.cond_order }, axis),
          yaxis: Object.assign({ title: "voom log-CPM" }, axis) }, cfg);
      } else {
        message("gl-expr", "No normalised expression for " + name + " (filtered out before limma, or not quantified).");
      }

      // 2. logFC per contrast; filled bars are significant
      if (D.lfc[i] && D.lfc[i].some(function (v) { return v !== null; })) {
        clear("gl-lfc");
        var sig = D.contrasts.map(function (ct, k) {
          var l = D.lfc[i][k], f = D.fdr[i][k];
          return l !== null && f !== null && f < D.fdr_de && Math.abs(l) > D.lfc_de;
        });
        var col = D.contrasts.map(function (ct, k) { var l = D.lfc[i][k]; return l === null ? D.col_ns : (l > 0 ? D.col_up : D.col_down); });
        Plotly.react("gl-lfc", [{
          type: "bar", x: D.contrasts, y: D.lfc[i],
          marker: { color: col.map(function (c, k) { return sig[k] ? c : "white"; }), line: { color: col, width: 2 } },
          text: D.contrasts.map(function (ct, k) { return "FDR " + fmtP(D.fdr[i][k]) + (sig[k] ? " (significant)" : ""); }),
          textposition: "none", hovertemplate: "%{x}<br>log2FC %{y:.2f}<br>%{text}<extra></extra>"
        }], { margin: { t: 10, b: 90 }, plot_bgcolor: "white",
          xaxis: Object.assign({ title: "", tickangle: -30, type: "category" }, axis),
          yaxis: Object.assign({ title: "log2 fold change (treatment vs control)", zeroline: true, zerolinecolor: "#333" }, axis) }, cfg);
      } else {
        message("gl-lfc", "No limma result for " + name + ".");
      }

      // 3. PSI per sample of each of the gene's top events, grouped by condition
      var evs = evByGene.get(i) || [];
      var tbl = document.getElementById("gl-ev-table");
      if (!D.ev || !evs.length) {
        message("gl-psi", "No significant splicing event for " + name + " (" + D.ds_rule + ").");
        tbl.innerHTML = "";
        return;
      }
      clear("gl-psi");
      var palette = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"];
      var traces = evs.map(function (j, n) {
        var lab = D.ev.type[j] + " " + D.ev.loc[j];
        return { type: "box", name: (n + 1) + ". " + lab, legendgroup: "e" + n, offsetgroup: "e" + n,
                 x: D.conds, y: D.ev.psi[j], text: D.samples, boxpoints: "all", jitter: 0.3, pointpos: 0,
                 marker: { color: palette[n % palette.length], size: 6 }, line: { color: palette[n % palette.length] },
                 hovertemplate: "%{text}<br>PSI %{y:.3f}<extra>" + (n + 1) + ". " + D.ev.type[j] + "</extra>" };
      });
      Plotly.react("gl-psi", traces, { boxmode: "group", margin: { t: 10, b: 90 }, plot_bgcolor: "white",
        legend: { orientation: "h", y: -0.25, font: { size: 11 } },
        xaxis: Object.assign({ title: "", type: "category", categoryorder: "array", categoryarray: D.cond_order }, axis),
        yaxis: Object.assign({ title: "PSI", range: [-0.05, 1.05] }, axis) }, cfg);

      // 4. table of the events
      var h = "<table><tr><th>#</th><th>Type</th><th>Location</th>" +
        D.contrasts.map(function (ct) { return "<th>" + ct + "</th>"; }).join("") + "</tr>";
      evs.forEach(function (j, n) {
        h += "<tr><td>" + (n + 1) + "</td><td>" + D.ev.type[j] + "</td><td>" + D.ev.loc[j] + "</td>";
        D.contrasts.forEach(function (ct, k) {
          var d = D.ev.dpsi[j][k], f = D.ev.fdr[j][k];
          if (d === null) { h += "<td>not tested</td>"; return; }
          var s = f !== null && f < D.fdr_ds && Math.abs(d) > D.dpsi_ds;
          h += '<td class="' + (s ? "sig" : "") + '">' + d.toFixed(3) + " (" + fmtP(f) + ")</td>";
        });
        h += "</tr>";
      });
      tbl.innerHTML = h + "</table>";
    }

    function lookup() {
      var k = index.get(input.value.trim().toUpperCase());
      if (k === undefined) { summary.innerHTML = '<span style="color:#b00">No gene "' + input.value + '" in the results.</span>'; return; }
      draw(k);
    }
    input.addEventListener("change", lookup);
    input.addEventListener("keydown", function (e) { if (e.key === "Enter") lookup(); });

    // The page is hidden until opened, and Plotly cannot size a hidden plot: draw on first visit
    function onShow() {
      if (!document.getElementById("gl-expr").offsetParent) return;
      if (current === null) draw(D.default_gene);
      else ["gl-expr", "gl-lfc", "gl-psi"].forEach(function (id) { try { Plotly.Plots.resize(id); } catch (e) {} });
    }
    window.addEventListener("hashchange", function () { setTimeout(onShow, 50); });
    setTimeout(onShow, 300);
  }
  if (document.readyState === "complete") start(); else window.addEventListener("load", start);
})();
