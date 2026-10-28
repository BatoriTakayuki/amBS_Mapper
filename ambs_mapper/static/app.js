"use strict";

const PARAMS = [
  ["min_coverage", "アンプリコン被覆率の下限(%)", 1],
  ["min_identity", "一致率の下限(%)（QUMA -p）", 1],
  ["max_mismatch", "ミスマッチ率の上限(%)（QUMA -mは個数）", 1],
  ["min_conversion", "CpH変換率の下限(%)（QUMA -c）", 1],
  ["max_unconverted", "未変換CpHの上限(個)（QUMA -u、空欄で制限なし）", 1],
  ["tag_max_errors", "タグの許容誤差(塩基)", 1],
  ["tag_slack", "タグ探索の余裕(塩基)", 1],
  ["snp_min_quality", "SNP位置の塩基品質の下限", 1],
  ["allele_min_fraction", "アレル判定に必要なSNP一致率(%)", 1],
  ["photo_rows", "写真に表示する最大リード数", 10],
  ["min_length", "最短リード長(bp)", 10],
  ["auto_min_reads", "自動サンプルにする最少合格リード数", 10],
  ["exclude_primer_cpg", "プライマー内のCpGを除外", null],
  ["single_index", "Nextera: 片側のインデックスだけのリードも、サンプルが一意なら使う", null],
];
const STATUS = { queued: "待機中", running: "実行中", completed: "完了", failed: "失敗", cancelled: "中止" };
const MODES = { analyze: "解析", check: "確認", setup: "参照準備", explore: "構造調査" };
const STORE = "ambs-mapper-form";
let health = null;
let inputs = { files: [], folders: [] };
let selected = new Set();
let shownJob = null;

function el(tag, attrs = {}, children = []) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else if (value !== undefined && value !== null && value !== false) node.setAttribute(key, value === true ? "" : value);
  }
  for (const child of [].concat(children)) if (child !== null && child !== undefined) node.append(child);
  return node;
}

async function api(path, options = {}) {
  const response = await fetch(path, options.body === undefined ? options : {
    ...options, method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(options.body) });
  const type = response.headers.get("Content-Type") || "";
  const data = type.includes("json") ? await response.json() : await response.text();
  if (!response.ok) throw new Error((data && data.error) || `HTTP ${response.status}`);
  return data;
}

function size(bytes) {
  if (bytes > 1e9) return (bytes / 1e9).toFixed(1) + " GB";
  if (bytes > 1e6) return (bytes / 1e6).toFixed(1) + " MB";
  return Math.max(1, Math.round(bytes / 1e3)) + " kB";
}

/* ---------- inputs ---------- */
async function loadInputs() {
  try { inputs = await api("/api/inputs"); } catch (error) { inputs = { files: [], folders: [] }; }
  renderInputs();
  refreshFastqOptions();
}

function renderInputs() {
  const box = document.getElementById("inputs");
  const filter = document.getElementById("inputs-filter").value.trim().toLowerCase();
  box.replaceChildren();
  const entries = [...inputs.folders.map(path => ({ path: path + "/", folder: true })), ...inputs.files];
  if (!entries.length) {
    box.append(el("div", { class: "hint", text: "data\\input にファイルがありません。" }));
    return;
  }
  for (const entry of entries) {
    if (filter && !entry.path.toLowerCase().includes(filter)) continue;
    const key = entry.folder ? entry.path.slice(0, -1) : entry.path;
    const box_ = el("input", { type: "checkbox" });
    box_.checked = selected.has(key);
    box_.addEventListener("change", () => { box_.checked ? selected.add(key) : selected.delete(key); save(); });
    box.append(el("label", { class: entry.folder ? "folder" : "" },
      [box_, " " + entry.path + (entry.folder ? "（フォルダ内すべて）" : `  (${size(entry.size)})`)]));
  }
}

/* ---------- tables ---------- */
function amplicons() {
  return [...document.querySelectorAll("#amplicons tbody tr")].map(row => {
    const [name, forward, reverse, max, region] = row.querySelectorAll("input");
    return { name: name.value.trim(), forward: forward.value.trim(), reverse: reverse.value.trim(), max_size: Number(max.value) || 5000,
             region: region.value.trim(), strand: row.querySelector("select").value };
  });
}

function addAmplicon(values = {}) {
  const remove = el("button", { class: "small danger", text: "削除" });
  const strand = el("select", {}, [el("option", { value: "+", text: "+（上鎖）" }), el("option", { value: "-", text: "-（下鎖）" })]);
  strand.value = values.strand || "+";
  strand.addEventListener("change", save);
  const row = el("tr", {}, [
    el("td", {}, el("input", { value: values.name || "", placeholder: "例: IG-DMR" })),
    el("td", { class: "primer" }, el("input", { value: values.forward || "", placeholder: "5'-...-3'" })),
    el("td", { class: "primer" }, el("input", { value: values.reverse || "", placeholder: "5'-...-3'" })),
    el("td", {}, el("input", { type: "number", min: 50, max: 50000, value: values.max_size || 5000 })),
    el("td", { class: "primer" }, el("input", { value: values.region || "", placeholder: "例: chr12:109494123-109494444" })),
    el("td", {}, strand),
    el("td", {}, remove)]);
  remove.addEventListener("click", () => { row.remove(); refreshAmpliconOptions(); save(); });
  row.addEventListener("input", () => { refreshAmpliconOptions(); save(); });
  document.querySelector("#amplicons tbody").append(row);
  refreshAmpliconOptions();
}

const SAMPLE_KEYS = ["name", "group", "amplicon", "tag_f", "tag_r", "mother", "father", "fastq"];

function samples() {
  return [...document.querySelectorAll("#samples tbody tr")].map(row => {
    const fields = row.querySelectorAll("input, select");
    return Object.fromEntries(SAMPLE_KEYS.map((key, i) => [key, fields[i].value.trim()]));
  }).filter(s => ["name", "group", "tag_f", "tag_r", "mother", "father"].some(key => s[key]));
}

function tagPlaceholder(key) {
  if (document.getElementById("layout").value !== "nextera") return "ACGT...";
  return key === "tag_f" ? "例: S520" : "例: N712";
}

function applyLayout() {
  const layout = document.getElementById("layout").value;
  for (const node of document.querySelectorAll("[data-layout]")) node.classList.toggle("hidden", node.dataset.layout !== layout);
  document.getElementById("th-tag-f").textContent = layout === "nextera" ? "i5" : "タグF";
  document.getElementById("th-tag-r").textContent = layout === "nextera" ? "i7" : "タグR";
  for (const row of document.querySelectorAll("#samples tbody tr")) {
    const fields = row.querySelectorAll("input, select");
    fields[SAMPLE_KEYS.indexOf("tag_f")].placeholder = tagPlaceholder("tag_f");
    fields[SAMPLE_KEYS.indexOf("tag_r")].placeholder = tagPlaceholder("tag_r");
  }
}

function addSample(values = {}) {
  const remove = el("button", { class: "small danger", text: "削除" });
  const cells = SAMPLE_KEYS.map(key => {
    if (key === "amplicon" || key === "fastq") {
      const select = el("select", { "data-key": key, "data-value": values[key] || "" });
      return el("td", {}, select);
    }
    const input = el("input", { value: values[key] || "" });
    if (key === "mother" || key === "father") { input.setAttribute("list", "strain-list"); input.placeholder = key === "mother" ? "例: B6" : "例: JF1"; }
    if (key === "tag_f" || key === "tag_r") input.placeholder = tagPlaceholder(key);
    return el("td", { class: key.startsWith("tag") ? "tag" : "" }, input);
  });
  const row = el("tr", {}, [...cells, el("td", {}, remove)]);
  remove.addEventListener("click", () => { row.remove(); save(); });
  row.addEventListener("input", save);
  row.addEventListener("change", save);
  document.querySelector("#samples tbody").append(row);
  refreshAmpliconOptions();
  refreshFastqOptions();
}

function fillSelect(select, options) {
  const current = select.value || select.dataset.value || "";
  select.replaceChildren(...options.map(([value, label]) => el("option", { value, text: label })));
  if (options.some(([value]) => value === current)) select.value = current;
  select.dataset.value = select.value;
}

function refreshAmpliconOptions() {
  const names = amplicons().map(a => a.name).filter(Boolean);
  for (const select of document.querySelectorAll('#samples select[data-key="amplicon"]'))
    fillSelect(select, names.map(name => [name, name]));
}

function refreshFastqOptions() {
  const options = [["", "（すべて）"], ...inputs.folders.map(f => [f, f + "/"]), ...inputs.files.map(f => [f.path, f.path])];
  for (const select of document.querySelectorAll('#samples select[data-key="fastq"]')) {
    const value = select.value || select.dataset.value;
    if (value && !options.some(([v]) => v === value)) options.push([value, value + "（見つかりません）"]);
    fillSelect(select, options);
  }
}

function pasteSamples() {
  const text = document.getElementById("paste-text").value;
  const rows = text.split(/\r?\n/).map(line => line.split("\t")).filter(cells => cells.some(c => c.trim()));
  if (rows.length && /サンプル|name/i.test(rows[0][0])) rows.shift();
  if (!rows.length) return;
  document.querySelector("#samples tbody").replaceChildren();
  for (const cells of rows) addSample(Object.fromEntries(SAMPLE_KEYS.map((key, i) => [key, (cells[i] || "").trim()])));
  document.getElementById("paste-box").classList.add("hidden");
  save();
}

/* ---------- params & persistence ---------- */
function renderParams(values) {
  const grid = document.getElementById("params");
  grid.replaceChildren();
  for (const [key, label, step] of PARAMS) {
    const input = step === null ? el("input", { type: "checkbox", "data-param": key })
      : el("input", { type: "number", min: 0, step, "data-param": key });
    if (step === null) input.checked = values[key] !== false;
    else input.value = values[key] === null || values[key] === undefined ? "" : values[key];
    input.addEventListener("change", save);
    grid.append(el("label", {}, [label, input]));
  }
}

function params() {
  const result = {};
  for (const input of document.querySelectorAll("[data-param]")) {
    const key = input.dataset.param;
    result[key] = input.type === "checkbox" ? input.checked : (input.value === "" ? (key === "max_unconverted" ? null : undefined) : Number(input.value));
    if (result[key] === undefined) delete result[key];
  }
  return result;
}

function formState() {
  return { title: document.getElementById("title").value.trim(), threads: Number(document.getElementById("threads").value) || undefined,
           inputs: [...selected], amplicons: amplicons(), samples: samples(), params: params(),
           layout: document.getElementById("layout").value,
           auto_mother: document.getElementById("auto-mother").value.trim(),
           auto_father: document.getElementById("auto-father").value.trim() };
}

function save() {
  try { localStorage.setItem(STORE, JSON.stringify(formState())); } catch (error) { /* storage unavailable */ }
}

function restore(state) {
  document.getElementById("title").value = state.title || "";
  if (state.threads) document.getElementById("threads").value = state.threads;
  document.getElementById("layout").value = state.layout || "inline";
  document.getElementById("auto-mother").value = (state.auto && state.auto.mother) || state.auto_mother || "";
  document.getElementById("auto-father").value = (state.auto && state.auto.father) || state.auto_father || "";
  selected = new Set(state.inputs || []);
  document.querySelector("#amplicons tbody").replaceChildren();
  document.querySelector("#samples tbody").replaceChildren();
  for (const a of state.amplicons && state.amplicons.length ? state.amplicons : [{}]) addAmplicon(a);
  for (const s of state.samples && state.samples.length ? state.samples : [{}]) addSample(s);
  renderParams({ ...(health ? health.defaults : {}), ...(state.params || {}) });
  applyLayout();
  renderInputs();
  refreshFastqOptions();
}

/* ---------- jobs ---------- */
async function submit(mode) {
  const error = document.getElementById("form-error");
  error.textContent = "";
  try {
    const job = await api("/api/jobs", { body: { ...formState(), mode } });
    save();
    await loadJobs();
    if (mode !== "setup") showJob(job.id);
  } catch (e) { error.textContent = e.message; }
}

function progressText(job) {
  const p = job.progress || {};
  if (job.status === "failed") return job.error || "";
  if (p.stage === "primer_search" && p.total) return `ゲノム検索 ${p.chrom}番 (${p.done}/${p.total})`;
  if (p.stage === "download" && p.total) return `ダウンロード ${Math.round(100 * p.done / p.total)}%`;
  return job.message || "";
}

async function loadJobs() {
  let jobs = [];
  try { jobs = await api("/api/jobs"); } catch (error) { return; }
  const body = document.getElementById("jobs");
  body.replaceChildren();
  for (const job of jobs) {
    const actions = [];
    actions.push(el("button", { class: "small secondary", text: "表示", onclick: () => showJob(job.id) }));
    if (job.status === "queued" || job.status === "running")
      actions.push(el("button", { class: "small danger", text: "中止", onclick: async () => { await api(`/api/jobs/${job.id}/cancel`, { body: {} }); loadJobs(); } }));
    else
      actions.push(el("button", { class: "small danger", text: "削除", onclick: async () => {
        if (!confirm("この解析結果を削除しますか？")) return;
        await api(`/api/jobs/${job.id}/delete`, { body: {} });
        if (shownJob === job.id) document.getElementById("result").classList.add("hidden");
        loadJobs(); } }));
    if (job.mode !== "setup" && job.mode !== "explore")
      actions.push(el("button", { class: "small secondary", text: "設定を復元", onclick: async () => {
        const request = await api(`/api/jobs/${job.id}/request`);
        restore({ ...request, threads: request.threads }); save(); window.scrollTo(0, 0); } }));
    body.append(el("tr", {}, [
      el("td", { text: new Date(job.created).toLocaleString("ja-JP") }),
      el("td", { text: job.title || "-" }),
      el("td", { text: MODES[job.mode] || job.mode }),
      el("td", {}, el("span", { class: "badge " + job.status, text: STATUS[job.status] || job.status })),
      el("td", { class: job.status === "failed" ? "error" : "", text: progressText(job) }),
      el("td", {}, actions)]));
    if (shownJob === job.id && job.status !== (document.getElementById("result").dataset.status || "")) showJob(job.id, true);
  }
}

/* ---------- results ---------- */
function table(rows, columns) {
  if (!rows.length) return el("p", { class: "hint", text: "該当なし" });
  columns = columns || [...new Set(rows.flatMap(r => Object.keys(r)))];
  return el("div", { class: "scroll" }, el("table", {}, [
    el("thead", {}, el("tr", {}, columns.map(c => el("th", { text: c })))),
    el("tbody", {}, rows.map(r => el("tr", {}, columns.map(c => el("td", { text: r[c] === null || r[c] === undefined ? "" : String(r[c]) })))))]));
}

function fileUrl(job, name, download) { return `/api/jobs/${job}/files/${name.split("/").map(encodeURIComponent).join("/")}${download ? "?download=1" : ""}`; }

async function showJob(id, quiet) {
  shownJob = id;
  const section = document.getElementById("result");
  const body = document.getElementById("result-body");
  section.classList.remove("hidden");
  const job = await api(`/api/jobs/${id}`);
  section.dataset.status = job.status;
  document.getElementById("result-title").textContent = `結果: ${job.title || MODES[job.mode]}（${STATUS[job.status]}）`;
  const parts = [];
  if (job.status === "completed" && job.mode !== "setup") {
    try { parts.push(...renderResults(id, await api(`/api/jobs/${id}/results`))); }
    catch (error) { parts.push(el("p", { class: "error", text: error.message })); }
  } else if (job.error) parts.push(el("p", { class: "error", text: job.error }));
  const log = await api(`/api/jobs/${id}/log`);
  parts.push(el("details", job.status === "completed" ? {} : { open: true }, [el("summary", { text: "ログ" }), el("pre", { class: "log", text: log.text || "" })]));
  body.replaceChildren(...parts);
  if (!quiet) section.scrollIntoView({ behavior: "smooth" });
}

function sampleRowsAreEmpty() { return !samples().length; }

function importCandidate(c) {
  const names = amplicons().map(a => a.name);
  if (names.includes(c.name)) return alert(`アンプリコン「${c.name}」はすでに表にあります。`);
  const empty = [...document.querySelectorAll("#amplicons tbody tr")].find(row => [...row.querySelectorAll("input")].slice(0, 3).every(i => !i.value.trim()));
  if (empty) empty.remove();
  addAmplicon({ name: c.name, forward: c.forward, reverse: c.reverse, max_size: Math.max(5000, (c.size || 0) + 500),
                region: c.chrom ? `chr${c.chrom}:${c.start}-${c.end}` : "", strand: c.strand || "+" });
  save();
}

function importCombos(c) {
  const combos = c.combos.filter(r => r.i5 !== "なし" && r.i7 !== "なし" && /^S5/.test(r.i5) && /^N7/.test(r.i7));
  if (!combos.length) return alert("両側のインデックスがそろった組み合わせがありません。");
  if (sampleRowsAreEmpty()) document.querySelector("#samples tbody").replaceChildren();
  document.getElementById("layout").value = "nextera";
  if (!amplicons().some(a => a.name === c.name)) importCandidate(c);
  for (const r of combos) addSample({ name: `${r.i5}-${r.i7}`, amplicon: c.name, tag_f: r.i5, tag_r: r.i7 });
  applyLayout();
  save();
  alert(`${combos.length}サンプルを表に追加しました。サンプル名・群・母系統・父系統を入力してください。`);
}

function renderExplore(id, result) {
  const parts = [];
  if (result.warnings.length) parts.push(el("div", {}, result.warnings.map(w => el("p", { class: "warn", text: "注意: " + w }))));
  parts.push(el("div", { class: "row" }, [
    el("a", { href: fileUrl(id, "explore.xlsx", true) }, el("button", { class: "secondary", text: "調査結果(Excel)" })),
    el("span", { class: "hint", text: `調べたリード ${result.info.reads.toLocaleString()}${result.info.limited ? "（上限で打ち切り）" : ""}、両端アダプターのインサート ${result.info.inserts.toLocaleString()}` })]));
  parts.push(el("h3", { text: "リードの構造" }));
  parts.push(table(result.structure));
  parts.push(el("h3", { text: "アンプリコン候補" }));
  parts.push(el("p", { class: "hint", text: "推定プライマーはインサート両端の共通配列です（最長22塩基）。実際のプライマー配列が分かる場合はそちらを使ってください。" }));
  for (const c of result.candidates) {
    const strand = c.strand ? (c.strand === "+" ? "上鎖" : "下鎖") : "未変換のため鎖なし";
    const where = c.chrom ? `chr${c.chrom}:${c.start.toLocaleString()}-${c.end.toLocaleString()}（${strand}、${c.size} bp、完全長リード ${c.full_length}%）` : "";
    const buttons = c.chrom && c.converted ? [
      el("button", { class: "small secondary", text: "アンプリコン表に追加", onclick: () => importCandidate(c) }),
      el("button", { class: "small secondary", text: "インデックスの組み合わせをサンプル表に追加", onclick: () => importCombos(c) })] : [];
    parts.push(el("div", { class: "candidate" }, [
      el("h4", { text: `${c.name}  ${c.reads.toLocaleString()}リード（${c.fraction}%）  ${c.converted ? "バイサルファイト変換済み" : "未変換"}` }),
      el("p", { text: `F: ${c.forward}  R: ${c.reverse}  インサート長: ${c.peaks.map(([l, n]) => `${l}bp(${n})`).join(", ")}` }),
      where ? el("p", { text: "ゲノム位置: " + where }) : null,
      c.note ? el("p", { class: "warn", text: c.note }) : null,
      el("div", { class: "row" }, buttons),
      el("details", {}, [el("summary", { text: `インデックスの組み合わせ（${c.combos.length}）` }),
        table(c.combos.map(r => ({ i5: r.i5, i7: r.i7, リード数: r.reads })))])]));
  }
  parts.push(el("h3", { text: "インデックスの組み合わせ（全リード、10リード以上）" }));
  parts.push(el("details", {}, [el("summary", { text: `${result.combos.length}組み合わせ（「なし」は片側のインデックスが読めないリード）` }),
    table(result.combos.map(r => ({ i5: r.i5, i7: r.i7, リード数: r.reads })))]));
  return parts;
}

function renderResults(id, result) {
  if (result.mode === "explore") return renderExplore(id, result);
  const parts = [];
  if (result.warnings && result.warnings.length)
    parts.push(el("div", {}, result.warnings.map(w => el("p", { class: "warn", text: "注意: " + w }))));
  parts.push(el("h3", { text: "アンプリコン" }));
  parts.push(table(result.amplicons.map(a => ({
    名前: a.name, 位置: `chr${a.chrom}:${a.start.toLocaleString()}-${a.end.toLocaleString()}`, 鎖: a.strand === "+" ? "上鎖(+)" : "下鎖(-)",
    長さ: a.size, CpG数: a.site_positions.length,
    "SNP（母/父: 使用可/全体）": a.pairs.map(p => `${p.mother}/${p.father}: ${p.informative}/${p.snps.length}`).join("、"),
    indel: a.pairs.map(p => p.indels.length).join("、") }))));
  if (result.mode === "check") {
    parts.push(el("p", { class: "hint", text: "プライマーの位置とSNPを確認しました。問題がなければ「解析開始」を押してください。" }));
    return parts;
  }
  const downloads = [["summary.xlsx", "集計表(Excel)"], ["photo.xlsx", "写真(Excel, QUMA-san形式)"], ["cpg_methylation.csv", "CpG別CSV"],
                     ["samples.csv", "サンプル別CSV"], ["reads.tsv.gz", "リードごとの判定"]];
  parts.push(el("div", { class: "row" }, [
    ...downloads.map(([name, label]) => el("a", { href: fileUrl(id, name, true) }, el("button", { class: "secondary", text: label }))),
    el("a", { href: `/api/jobs/${id}/zip` }, el("button", { text: "すべてZIPで保存" }))]));
  parts.push(el("p", { class: "hint", text: `総リード ${Number(result.info.reads).toLocaleString()}、リード処理 ${result.info.seconds} 秒` }));
  if (result.info.igv && globalThis.mountIGV) {
    parts.push(el("h3", { text: "IGV（サンプル・タグ別、タグ不明）" }));
    const holder = el("div", {});
    const open = el("button", { text: "IGVを表示" });
    open.addEventListener("click", () => { open.remove(); globalThis.mountIGV(holder, id); });
    parts.push(el("div", {}, [open, holder]));
  }
  if (result.info.bam) {
    parts.push(el("h3", { text: "IGVでリードを見る" }));
    const buttons = [];
    if (result.info.sample_bams)
      buttons.push(el("a", { href: fileUrl(id, "bam_by_sample.zip", true) },
        el("button", { text: `サンプル別BAM（合格リードのみ、${result.info.sample_bams}サンプル＋タグ不明、ZIP）` })));
    buttons.push(el("a", { href: fileUrl(id, "reads.bam", true) }, el("button", { class: "secondary", text: "全リードのBAM（reads.bam）" })),
      el("a", { href: fileUrl(id, "reads.bam.bai", true) }, el("button", { class: "secondary", text: "reads.bam.bai（索引）" })));
    parts.push(el("div", { class: "row" }, buttons));
    parts.push(el("ol", { class: "hint" }, [
      "サンプル別BAMのZIPを展開し、IGVでゲノム「Mouse (GRCm39/mm39)」を選んで、File → Load from File で .bam をまとめて選びます。サンプルごとに1トラックになります（トラック名＝ファイル名＝サンプル名）。染色体名 12 と chr12 はIGVが自動で対応づけます。",
      "リードを右クリック →「Color alignments by」→「bisulfite mode」→「CG」で、CpGがメチル化は赤、非メチル化は青になります。",
      "「タグ不明_アンプリコン名.bam」は、そのアンプリコンに載ったのにサンプルに割り当てられなかったリード（タグが読めない／複数サンプルに一致）です。品質判定の前に外れたリードなので、不合格のリードも含みます。区別は「Group alignments by」→「tag」で YT（tag_unknown／tag_ambiguous）。",
      "各リードに読めたインデックスを記録しています（Nexteraは i5・i7、インラインは tf・tr タグ）。値はインデックス名（S520 など）か、読めなかった理由: no_index＝インデックスが付いていない（アダプターの外にP5/P7がない）、unreadable＝付いているが読めない、ambiguous＝複数に同じくらい一致、no_adapter＝プライマーの外にNexteraアダプターがない、no_primer＝プライマーが見つからない。「Group alignments by」→「tag」で i5 や i7 と入れると、インデックス別に並びます。リードをクリックすると表示される情報でも確認できます。",
      "アレル別に分けるには、右クリック →「Group alignments by」→「phase」（1: 母、2: 父）。",
      "サンプル別BAMはメチル化率の計算に使った合格リードだけです。不合格やアンプリコン外のリードも見たいときは全リードのBAMを開き、「Group alignments by」→「sample」や「tag」（YT＝判定結果）で分けます。",
      "IGVは参照（B6）の配列でCpGかどうかを決めます。SNPで片方の系統にしかないCpGは、IGVの色とツールの結果（その系統のリードだけで計算）が違って見えます。"]
      .map(text => el("li", { text }))));
  }
  parts.push(el("h3", { text: "サンプル別" }));
  parts.push(table(result.samples));
  const unassigned = Object.entries(result.unassigned || {}).map(([amplicon, counts]) => ({ アンプリコン: amplicon, ...counts }));
  if (unassigned.length) {
    parts.push(el("h3", { text: "サンプルに割り当てられなかったリード" }));
    parts.push(table(unassigned));
  }
  if (result.off_target && result.off_target.length) {
    parts.push(el("h3", { text: "アンプリコン外にマッピングされたリード（多い領域、1 kb単位）" }));
    parts.push(el("details", {}, [el("summary", { text: `${result.off_target.length}領域` }), table(result.off_target)]));
  }
  for (const image of result.images) {
    parts.push(el("h3", { text: `${image.amplicon}: メチル化率` }));
    parts.push(el("div", { class: "plots" }, [el("img", { src: fileUrl(id, image.means), alt: "" }), el("img", { src: fileUrl(id, image.profile), alt: "" })]));
  }
  parts.push(el("h3", { text: "写真（●メチル化 ○非メチル化、灰色は判定不能）" }));
  parts.push(el("div", { class: "photos" }, result.groups.filter(g => g.image).map(g =>
    el("figure", {}, [el("a", { href: fileUrl(id, g.image), target: "_blank" }, el("img", { src: fileUrl(id, g.image), alt: g.sample, loading: "lazy" })),
      el("figcaption", { text: `${g.sample}（${g.group || "-"}）${g.allele} ${g.strain}  n=${g.reads}  ${g.methylation ?? "-"}%` })]))));
  return parts;
}

/* ---------- start ---------- */
async function loadHealth() {
  const status = document.getElementById("status");
  try {
    health = await api("/api/health");
    status.textContent = `参照ゲノム: ${health.reference_ready ? "準備済み" : "未準備"} ／ CPU ${health.logical_cpus} ／ 空き容量 ${health.disk_free_gb} GB`;
    document.getElementById("setup-box").classList.toggle("hidden", health.reference_ready);
    if (!document.getElementById("threads").value) document.getElementById("threads").value = health.logical_cpus;
  } catch (error) { status.textContent = "サーバーに接続できません。"; }
}

async function loadStrains() {
  try {
    const { strains } = await api("/api/strains");
    document.getElementById("strain-list").replaceChildren(...["B6", "JF1", ...strains].map(s => el("option", { value: s })));
  } catch (error) { /* strains are optional hints */ }
}

document.addEventListener("DOMContentLoaded", async () => {
  await loadHealth();
  let state = {};
  try { state = JSON.parse(localStorage.getItem(STORE) || "{}"); } catch (error) { state = {}; }
  await loadInputs();
  restore(state);
  loadStrains();
  loadJobs();
  setInterval(loadJobs, 3000);
  setInterval(loadHealth, 30000);
  document.getElementById("inputs-refresh").addEventListener("click", loadInputs);
  document.getElementById("inputs-filter").addEventListener("input", renderInputs);
  document.getElementById("amplicon-add").addEventListener("click", () => addAmplicon());
  document.getElementById("sample-add").addEventListener("click", () => addSample());
  document.getElementById("paste-toggle").addEventListener("click", () => document.getElementById("paste-box").classList.toggle("hidden"));
  document.getElementById("paste-apply").addEventListener("click", pasteSamples);
  document.getElementById("bulk-apply").addEventListener("click", () => {
    const mother = document.getElementById("bulk-mother").value.trim(), father = document.getElementById("bulk-father").value.trim();
    for (const row of document.querySelectorAll("#samples tbody tr")) {
      const fields = row.querySelectorAll("input, select");
      if (mother) fields[5].value = mother;
      if (father) fields[6].value = father;
    }
    save();
  });
  document.getElementById("layout").addEventListener("change", () => { applyLayout(); save(); });
  document.getElementById("auto-mother").addEventListener("input", save);
  document.getElementById("auto-father").addEventListener("input", save);
  document.getElementById("title").addEventListener("input", save);
  document.getElementById("threads").addEventListener("input", save);
  document.getElementById("check-button").addEventListener("click", () => submit("check"));
  document.getElementById("explore-button").addEventListener("click", () => submit("explore"));
  document.getElementById("run-button").addEventListener("click", () => submit("analyze"));
  document.getElementById("setup-button").addEventListener("click", () => submit("setup"));
  document.getElementById("export-button").addEventListener("click", () => {
    const blob = new Blob([JSON.stringify(formState(), null, 1)], { type: "application/json" });
    const link = el("a", { href: URL.createObjectURL(blob), download: "ambs_settings.json" });
    link.click();
    URL.revokeObjectURL(link.href);
  });
  document.getElementById("import-file").addEventListener("change", async event => {
    const file = event.target.files[0];
    if (!file) return;
    try { restore(JSON.parse(await file.text())); save(); }
    catch (error) { document.getElementById("form-error").textContent = "設定ファイルを読み込めませんでした。"; }
    event.target.value = "";
  });
});
