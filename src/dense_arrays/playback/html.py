"""Self-contained minimal playback canvas for dense-array realization plans."""

from __future__ import annotations

import html
import json
from collections.abc import Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from enum import Enum
from typing import Any

from .models import PlaybackPlan
from .positions import radial_path_positions
from .theme import PlaybackPresentation
from .timeline import complement_sequence, current_added_indices, revealed_indices


@dataclass(frozen=True, slots=True)
class PlaybackDocument:
    """One scene plus optional producer-rendered duplex frames."""

    plan: PlaybackPlan
    title: str
    subtitle: str = ""
    label_overrides: Mapping[str, str] = field(default_factory=dict)
    presentation: PlaybackPresentation = field(default_factory=PlaybackPresentation)
    duplex_svg_frames: tuple[str, ...] = ()


def _jsonable(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return {
            item.name: _jsonable(getattr(value, item.name)) for item in fields(value)
        }
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def _document_payload(document: PlaybackDocument) -> dict[str, Any]:
    payload = _jsonable(document)
    payload["graph_positions"] = radial_path_positions(len(document.plan.steps))
    payload["complement_sequence"] = complement_sequence(
        document.plan.realized_sequence
    )
    payload["timeline_frames"] = [
        {
            "revealed_indices": revealed_indices(document.plan.steps, step_index),
            "current_added_indices": current_added_indices(
                document.plan.steps, step_index
            ),
        }
        for step_index in range(len(document.plan.steps))
    ]
    return payload


def render_playback_html(
    plan: PlaybackPlan,
    *,
    title: str,
    subtitle: str = "",
    label_overrides: Mapping[str, str] | None = None,
) -> str:
    document = PlaybackDocument(
        plan=plan,
        title=title,
        subtitle=subtitle,
        label_overrides=dict(label_overrides or {}),
    )
    return render_playback_collection_html((document,), title=title)


def render_playback_collection_html(
    documents: tuple[PlaybackDocument, ...],
    *,
    title: str,
) -> str:
    if not documents:
        raise ValueError("at least one playback document is required")
    payload = json.dumps(
        [_document_payload(document) for document in documents],
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).replace("</", "<\\/")
    return _TEMPLATE.replace("__PAGE_TITLE__", html.escape(title, quote=True)).replace(
        "__PLAYBACK_PAYLOAD__",
        payload,
    )


_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=1200">
  <title>__PAGE_TITLE__</title>
  <style>
    :root { --ink:#4b5563; --line:#c3cad3; --active:#ef6c32; }
    * { box-sizing:border-box; }
    html,body { margin:0; min-height:100%; background:#ececec; color:var(--ink); font-family:"DejaVu Sans",sans-serif; }
    .review { min-width:1200px; padding:22px; }
    .canvas { width:min(1600px,calc(100vw - 44px)); aspect-ratio:16/9; margin:0 auto; display:grid; grid-template-columns:29% 71%; background:white; overflow:hidden; box-shadow:0 10px 34px rgba(0,0,0,.10); }
    .graph,.duplex { min-width:0; min-height:0; }
    .graph { padding:18px 2px 18px 18px; }
    .graph svg { width:100%; height:100%; display:block; overflow:visible; }
    .duplex { display:grid; place-items:center; padding:14px 18px 14px 4px; overflow:hidden; }
    .duplex > svg { display:block; width:100%; height:100%; max-width:100%; max-height:100%; }
    .edge { stroke:#c3cad3; stroke-width:2.2; fill:none; }
    .edge.done { stroke:#5f6b67; stroke-width:3.2; }
    .edge.hot { stroke:var(--active); stroke-width:5; }
    .node rect { fill:white; stroke:#c3cad3; stroke-width:2; rx:9; }
    .node.future { opacity:.24; }
    .node.done rect { stroke-width:2.5; }
    .node.hot rect { stroke:var(--active); stroke-width:5; }
    .node text { fill:white; font:700 15px "DejaVu Sans Mono",monospace; text-anchor:middle; dominant-baseline:middle; }
    .transport { width:min(1600px,calc(100vw - 44px)); margin:12px auto 0; display:grid; grid-template-columns:minmax(250px,.55fr) auto minmax(400px,1fr); gap:12px; align-items:center; }
    select,button { height:40px; border:1px solid #69716f; background:white; color:#343b39; font-weight:700; padding:0 12px; }
    button { min-width:82px; border-radius:999px; color:white; background:#343b39; }
    input { width:100%; accent-color:var(--active); }
  </style>
</head>
<body>
  <div class="review">
    <main class="canvas">
      <section class="graph" id="graph" aria-label="Realized path"></section>
      <section class="duplex" id="duplex" aria-label="Dense DNA array"></section>
    </main>
    <div class="transport">
      <select id="scene" aria-label="Select scene"></select>
      <button id="play" type="button">Play</button>
      <input id="timeline" type="range" min="0" value="0" aria-label="Playback step">
    </div>
  </div>
  <script id="playback-data" type="application/json">__PLAYBACK_PAYLOAD__</script>
  <script>
    const documents=JSON.parse(document.getElementById('playback-data').textContent);
    const colors=['#67BFA5','#D883A4','#7BA4D9','#C08A56','#5DA79F','#D1B06C','#74C0CB','#86A5D8'];
    const fixedUp='#7D86D1', fixedDown='#C886D1';
    const state={scene:0,step:0,timer:null};
    const el={graph:document.getElementById('graph'),duplex:document.getElementById('duplex'),scene:document.getElementById('scene'),play:document.getElementById('play'),timeline:document.getElementById('timeline')};
    const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
    function stepColor(step,index){ if(step.placement_kind==='fixed_element'){ const label=String(step.label??'').toLowerCase(); return label.includes('downstream')||label.includes('-10')?fixedDown:fixedUp; } return colors[index%colors.length]; }
    function stop(){ if(state.timer)clearInterval(state.timer); state.timer=null; el.play.textContent='Play'; }
    function nodeText(sequence,x,y){ const seq=String(sequence); if(seq.length<=12)return `<text x="${x}" y="${y}">${esc(seq)}</text>`; const cut=Math.ceil(seq.length/2); return `<text x="${x}" y="${y-10}">${esc(seq.slice(0,cut))}</text><text x="${x}" y="${y+11}">${esc(seq.slice(cut))}</text>`; }
    function renderGraph(doc){
      const steps=doc.plan.steps, positions=doc.graph_positions, width=500, height=900;
      let edges='',nodes='';
      for(let i=1;i<steps.length;i++){ const a=positions[i-1],b=positions[i]; const cls=i<state.step?'done':i===state.step?'hot':''; edges+=`<line class="edge ${cls}" x1="${a[0]*width}" y1="${a[1]*height}" x2="${b[0]*width}" y2="${b[1]*height}"/>`; }
      steps.forEach((step,index)=>{ const p=positions[index],x=p[0]*width,y=p[1]*height,color=stepColor(step,index),cls=index<state.step?'done':index===state.step?'hot':'future'; nodes+=`<g class="node ${cls}"><rect x="${x-94}" y="${y-39}" width="188" height="78" style="fill:${color};stroke:${index===state.step?'#ef6c32':color}"/>${nodeText(step.placement_sequence,x,y)}</g>`; });
      el.graph.innerHTML=`<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Realized k-mer path">${edges}${nodes}</svg>`;
    }
    function fallbackDuplex(doc){
      const plan=doc.plan,seq=plan.realized_sequence;
      const revealed=new Set(doc.timeline_frames[state.step].revealed_indices);
      const shown=[...seq].map((base,index)=>revealed.has(index)?base:'\u00b7').join('');
      const paired=[...doc.complement_sequence].map((base,index)=>revealed.has(index)?base:'\u00b7').join('');
      const bars=plan.steps.slice(0,state.step+1).map((step,index)=>{ const x=80+step.start/seq.length*1040,w=(step.end-step.start)/seq.length*1040,y=step.orientation==='rev'?435:150; return `<rect x="${x}" y="${y}" width="${w}" height="62" rx="7" fill="${stepColor(step,index)}"/><text x="${x+w/2}" y="${y+39}" text-anchor="middle" fill="white" font-family="DejaVu Sans Mono" font-size="22">${esc(step.placement_sequence)}</text>`; }).join('');
      return `<svg viewBox="0 0 1200 650"><text x="20" y="297" font-size="30" fill="#4b5563">5'</text><text x="20" y="362" font-size="30" fill="#4b5563">3'</text><text x="80" y="300" textLength="1040" lengthAdjust="spacing" font-family="DejaVu Sans Mono" font-size="22" fill="#4b5563">${shown}</text><text x="80" y="365" textLength="1040" lengthAdjust="spacing" font-family="DejaVu Sans Mono" font-size="22" fill="#4b5563">${paired}</text>${bars}</svg>`;
    }
    function render(){ const doc=documents[state.scene]; state.step=Math.max(0,Math.min(state.step,doc.plan.steps.length-1)); el.timeline.max=String(doc.plan.steps.length-1); el.timeline.value=String(state.step); renderGraph(doc); el.duplex.innerHTML=doc.duplex_svg_frames?.[state.step]||fallbackDuplex(doc); }
    documents.forEach((doc,index)=>{ const option=document.createElement('option'); option.value=String(index); option.textContent=doc.title; el.scene.appendChild(option); });
    el.scene.addEventListener('change',event=>{ stop();state.scene=Number(event.target.value);state.step=0;render(); });
    el.timeline.addEventListener('input',event=>{ stop();state.step=Number(event.target.value);render(); });
    el.play.addEventListener('click',()=>{ if(state.timer){stop();return;} el.play.textContent='Pause'; state.timer=setInterval(()=>{ const plan=documents[state.scene].plan; if(state.step<plan.steps.length-1){state.step+=1;render();return;} if(state.scene<documents.length-1){state.scene+=1;el.scene.value=String(state.scene);state.step=0;render();return;} stop(); },1050); });
    render();
  </script>
</body>
</html>
"""
