"use strict";
/* Embedded IGV (igv.js 3.0.0, served locally) for a finished analysis: one track per
   sample and one "タグ不明" track per amplicon, switched with checkboxes. */
(() => {
  let library = null;
  let current = { browser: null, generation: 0 };

  function loadLibrary() {
    if (globalThis.igv) return Promise.resolve(globalThis.igv);
    if (!library) library = new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.src = "/vendor/igv-3.0.0.min.js";
      script.onload = () => resolve(globalThis.igv);
      script.onerror = () => { library = null; script.remove(); reject(new Error("IGVライブラリを読み込めませんでした。")); };
      document.head.append(script);
    });
    return library;
  }

  function node(tag, attrs = {}, children = []) {
    const element = document.createElement(tag);
    for (const [key, value] of Object.entries(attrs)) {
      if (key === "text") element.textContent = value;
      else if (key === "class") element.className = value;
      else element.setAttribute(key, value);
    }
    element.append(...[].concat(children).filter(Boolean));
    return element;
  }

  function trackConfig(track) {
    const { amplicon, unknown, ...config } = track;   // igv.js only gets its own options
    return config;
  }

  globalThis.mountIGV = (container, jobId) => {
    const generation = ++current.generation;
    if (current.browser && globalThis.igv) { globalThis.igv.removeBrowser(current.browser); current.browser = null; }
    const status = node("p", { class: "hint", role: "status", text: "IGVとデータを読み込んでいます…" });
    const controls = node("div", { class: "igv-controls" });
    const viewer = node("div", { class: "igv-viewer" });
    container.replaceChildren(
      node("p", { class: "hint", text: "CpGは赤＝メチル化、青＝非メチル化（写真と同じ判定）。リードの端の色付きの塩基は、ゲノムに並ばなかった部分（タグ・アダプターなど）です。リードをクリックすると、読めたインデックス（i5/i7、読めなければ理由）や判定（YT）が出ます。非メチル化CpGのTはこの表示に限りCとして描いています。" }),
      controls, status, viewer);

    (async () => {
      try {
        const response = await fetch(`/api/jobs/${jobId}/igv/config`);
        const data = await response.json();
        if (!response.ok) throw new Error(data.error || "表示データがありません。");
        const igv = await loadLibrary();
        if (generation !== current.generation) return;
        const browser = await igv.createBrowser(viewer, data.config);
        if (generation !== current.generation) { igv.removeBrowser(browser); return; }
        current.browser = browser;

        const jump = node("div", { class: "row" }, [node("span", { class: "hint", text: "アンプリコン:" }),
          ...data.amplicons.map(a => {
            const button = node("button", { class: "small secondary", type: "button", text: a.name });
            button.addEventListener("click", () => browser.search(a.locus));
            return button;
          })]);
        const groups = new Map();
        for (const track of data.alignments) {
          if (!groups.has(track.amplicon)) groups.set(track.amplicon, []);
          groups.get(track.amplicon).push(track);
        }
        const lists = [];
        let shown = 0;
        for (const [amplicon, tracks] of groups) {
          tracks.sort((a, b) => (b.unknown - a.unknown) || a.name.localeCompare(b.name, "ja"));
          const boxes = tracks.map(track => {
            const box = node("input", { type: "checkbox" });
            // start with the tag-unknown track and the first samples of the first amplicon
            box.checked = amplicon === data.alignments[0].amplicon && (track.unknown || shown++ < 3);
            if (box.checked) browser.loadTrack(trackConfig(track));
            box.addEventListener("change", () => box.checked
              ? browser.loadTrack(trackConfig(track)) : browser.removeTrackByName(track.name));
            return node("label", { class: track.unknown ? "igv-unknown" : "" }, [box, " " + track.name]);
          });
          lists.push(node("details", { open: "" }, [node("summary", { text: `${amplicon}（${tracks.length}トラック）` }),
            node("div", { class: "igv-tracks" }, boxes)]));
        }
        controls.replaceChildren(jump, ...lists);
        status.textContent = "ドラッグで移動、＋／−で拡大・縮小。表示するサンプル（タグ）はチェックで切り替えます。";
      } catch (error) {
        if (generation === current.generation) status.textContent = `IGVを表示できません：${error.message}`;
      }
    })();
  };
})();
