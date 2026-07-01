import os, sys
APP_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, APP_DIR)           # so `from predictor import ...` works regardless of cwd
def DATA(rel): return os.path.join(APP_DIR, rel)

import streamlit as st
import pandas as pd, numpy as np
import plotly.express as px, plotly.graph_objects as go
from rdkit import Chem
from rdkit.Chem import Draw, DataStructs
from predictor import HuKSARx, morgan_fp

try:
    from streamlit_ketcher import st_ketcher
    KETCHER_OK = True
except Exception:
    KETCHER_OK = False

_FAVICON = DATA("assets/huksa_logo.png")
st.set_page_config(page_title="HuKSA", layout="wide",
                   page_icon=_FAVICON if os.path.exists(_FAVICON) else "🎯")

st.markdown("""
<style>
  #MainMenu, footer {visibility: hidden;}
  .block-container {padding-top: 2.2rem; max-width: 1240px;}
  h1 {font-weight: 700; letter-spacing: -0.5px; margin-bottom: 0.2rem;}
  .huksa-tagline {color:#5B6670; font-size:0.95rem; line-height:1.45; margin-top:0.3rem;}
  .stButton > button {border-radius: 10px; font-weight: 600;}
  .rx-card {background:#fff; border:1px solid #E6EAEF; border-radius:14px; padding:18px 20px;}
  .rx-bar {height:14px; background:#EEF1F4; border-radius:3px;}
  .rx-fill {display:block; height:100%; border-radius:3px;}
</style>
""", unsafe_allow_html=True)

h_icon, h_text = st.columns([1, 11], vertical_alignment="center")
with h_icon:
    st.markdown(
        "<svg width='78' height='78' viewBox='0 0 72 72' style='display:block;margin-top:-6px;'>"
        "<g fill='none' stroke='#2c6fb0' stroke-width='2.2' stroke-linecap='round'>"
        "<line x1='36' y1='1' x2='36' y2='6'/><line x1='36' y1='66' x2='36' y2='71'/>"
        "<line x1='1' y1='36' x2='6' y2='36'/><line x1='66' y1='36' x2='71' y2='36'/></g>"
        "<circle cx='36' cy='36' r='30' fill='none' stroke='#2c6fb0' stroke-width='3.2'/>"
        "<circle cx='36' cy='36' r='21' fill='none' stroke='#27ae9e' stroke-width='3.2'/>"
        "<circle cx='36' cy='36' r='12' fill='none' stroke='#2c6fb0' stroke-width='2.8'/>"
        "<circle cx='59' cy='13' r='2.8' fill='#e8a32a'/>"
        "<circle cx='36' cy='36' r='5' fill='#e8a32a'/></svg>",
        unsafe_allow_html=True)
with h_text:
    st.title("HuKSA — A transparent kinome-wide selectivity profiling tool")
    st.markdown(
        "<p class='huksa-tagline'>A ligand-based method that predicts a molecule's likely kinase "
        "targets and selectivity directly from its structure. For any SMILES, HuKSA infers the probable kinase "
        "targets from the measured activity of the query's nearest structural analogues among 794 kinome-profiled "
        "compounds (464 kinases), and reports a coverage-robust selectivity score (S-score) together with an "
        "explicit confidence flag (High, Moderate, Low or Outlier), so that every prediction is accompanied by "
        "a measure of how far it can be trusted.</p>", unsafe_allow_html=True)

# ---------------------------------------------------------------- assets
@st.cache_resource(show_spinner="loading read-across engine…")
def load_assets():
    rx = HuKSARx(DATA("data"))
    atlas = pd.read_csv(DATA("data/atlas_background_data.csv"))
    # fingerprints of the map's reference points, so a query can be placed at its closest analog
    # (the map is visualization-only; no UMAP model is loaded, keeping the app lightweight/portable)
    atlas_fps = [morgan_fp(s) for s in atlas["SMILES"]]
    return rx, atlas, atlas_fps

rx, ATLAS, ATLAS_FPS = load_assets()

CONF = {"High":("#EAF3DE","#27500A"), "Moderate":("#FAEEDA","#633806"),
        "Low":("#FCEBEB","#791F1F"), "Outlier":("#FCEBEB","#791F1F")}
TIER = {"Highly Selective":("#E6F1FB","#0C447C"), "Selective":("#EAF3DE","#27500A"),
        "Promiscuous":("#FAEEDA","#633806"), "Unknown":("#F1EFE8","#444441")}

def query_map_xy(smiles):
    """Place the query on the UMAP map at its nearest reference analog's coordinates."""
    q = morgan_fp(smiles)
    if q is None:
        return None
    sims = [DataStructs.TanimotoSimilarity(q, f) if f is not None else 0.0 for f in ATLAS_FPS]
    i = int(np.argmax(sims))
    return (float(ATLAS.iloc[i]["UMAP_1"]) + 0.12, float(ATLAS.iloc[i]["UMAP_2"]) + 0.12)

def bar(pact):
    w = max(2, min(100, (pact - 5.5) / 3.7 * 100))
    return f"<span class='rx-bar'><span class='rx-fill' style='width:{w:.0f}%; background:#85B7EB;'></span></span>"

def result_card(res, top_n=8):
    ad, sel = res["applicability_domain"], res["selectivity"]
    cb, ct = CONF.get(ad["confidence"], CONF["Low"])
    tb, tt = TIER.get(sel["tier"], TIER["Unknown"])
    low = ad["confidence"] in ("Low", "Outlier")
    rows = ""
    for t in res["predicted_targets"][:top_n]:
        rows += (f"<div style='display:grid; grid-template-columns:96px 1fr 40px 56px; align-items:center; gap:10px; font-size:13px; margin:4px 0;'>"
                 f"<span style='font-weight:500;'>{t['gene']}</span>{bar(t['pred_pActivity'])}"
                 f"<span style='text-align:right; color:#5B6670;'>{t['pred_pActivity']:.1f}</span>"
                 f"<span style='text-align:right; color:#9aa6b4; font-size:11px;'>n={t['neighbours_measured']}</span></div>")
    warn = (f"<div style='font-size:12.5px; background:#FCEBEB; color:#791F1F; border-radius:8px; padding:9px 11px; margin:10px 0;'>"
            f"⚠ outside the applicability domain — the target ranking is unreliable; treat as a lead-only hint.</div>") if low else ""
    sval = f"{sel['s_score']:.3f}" if sel['s_score'] is not None else "—"
    st.markdown(
        f"<div class='rx-card'>"
        f"<div style='display:flex; justify-content:space-between; align-items:center; margin-bottom:10px;'>"
        f"<span style='font-size:15px; font-weight:600;'>prediction</span>"
        f"<span style='font-size:12px; padding:4px 12px; border-radius:8px; background:{cb}; color:{ct};'>{ad['confidence'].lower()} confidence</span></div>"
        f"<div style='display:flex; gap:20px; flex-wrap:wrap; font-size:13px; color:#5B6670; margin-bottom:12px;'>"
        f"<span>nearest analog · tanimoto <b style='color:#1A2027;'>{ad['nearest_neighbour_tanimoto']:.2f}</b></span>"
        f"<span>selectivity · S <b style='color:#1A2027;'>{sval}</b> "
        f"<span style='font-size:12px; padding:2px 9px; border-radius:8px; background:{tb}; color:{tt};'>{sel['tier'].lower()}</span></span></div>"
        f"{warn}"
        f"<div style='font-size:12px; color:#9aa6b4; margin-bottom:6px;'>predicted targets · weighted measured pActivity (k=5)</div>"
        f"{rows}</div>", unsafe_allow_html=True)

def atlas_map(query_xy=None):
    fig = px.scatter(ATLAS, x="UMAP_1", y="UMAP_2", color="S_Score",
                     color_continuous_scale="RdYlBu_r", opacity=0.35, template="plotly_white",
                     title="selectivity landscape  (query ✕ placed at its closest analog)")
    fig.update_traces(marker=dict(size=5, line=dict(width=0)))
    if query_xy is not None:
        fig.add_trace(go.Scatter(x=[query_xy[0]], y=[query_xy[1]], mode="markers", name="query",
                      marker=dict(color="white", size=16, symbol="x", line=dict(width=2, color="black"))))
    fig.update_layout(height=560, plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
                      margin=dict(l=20, r=20, t=50, b=20),
                      coloraxis_colorbar=dict(title="S-score<br>(low=selective)", thickness=14),
                      legend=dict(y=0.99, x=0.01, bgcolor="rgba(255,255,255,0.6)"))
    fig.update_xaxes(showgrid=True, gridcolor="rgba(0,0,0,0.05)", zeroline=False, title="UMAP 1")
    fig.update_yaxes(showgrid=True, gridcolor="rgba(0,0,0,0.05)", zeroline=False, title="UMAP 2")
    return fig

# ---------------------------------------------------------------- UI
st.divider()
EX = {"Tirabrutinib (BTK)":"CC#CC(=O)N1CC[C@H](C1)N2C3=NC=NC(=C3N(C2=O)C4=CC=C(C=C4)OC5=CC=CC=C5)N",
      "Sitravatinib (MET/AXL)":"COCCNCC1=CN=C(C=C1)C2=CC3=NC=CC(=C3S2)OC4=C(C=C(C=C4)NC(=O)C5(CC5)C(=O)NC6=CC=C(C=C6)F)F",
      "KU-0063794 (mTOR)":"C[C@@H]1CN(C[C@@H](O1)C)C2=NC3=C(C=CC(=N3)C4=CC(=C(C=C4)OC)CO)C(=N2)N5CCOCC5"}

tab1, tab2 = st.tabs(["🔬  Single molecule", "📚  Virtual screening"])

# ============================ SINGLE MOLECULE ============================
with tab1:
    c1, c2 = st.columns([3, 2])
    with c1:
        st.markdown("##### Query a molecule")
        ex = st.selectbox("load an example (optional)", ["—"] + list(EX), index=0)
        default = EX[ex] if ex != "—" else ""
        if KETCHER_OK:
            mode = st.radio("input mode", ["Paste SMILES", "Draw molecule"], horizontal=True)
            if mode == "Draw molecule":
                smiles = (st_ketcher(default, height=380, key="ketcher") or "").strip()
                st.caption(f"SMILES: `{smiles}`" if smiles else "draw a structure, then click **Apply** in the editor")
            else:
                smiles = st.text_input("SMILES", value=default, placeholder="paste a SMILES string…")
        else:
            smiles = st.text_input("SMILES", value=default, placeholder="paste a SMILES string…")
        go_btn = st.button("Predict", type="primary", use_container_width=True)
    if go_btn and smiles.strip():
        res = rx.predict(smiles.strip())
        if not res.get("valid"):
            st.error("Could not parse that SMILES.")
        else:
            with c2:
                m = Chem.MolFromSmiles(smiles.strip())
                if m is not None:
                    st.image(Draw.MolToImage(m, size=(300, 230)), caption="query molecule")
            st.divider()
            left, right = st.columns([1, 1])
            with left:
                result_card(res)
            with right:
                xy = query_map_xy(smiles.strip())
                st.plotly_chart(atlas_map(xy), use_container_width=True)
            with st.expander("full ranked target list (top 10) + how this works"):
                st.dataframe(pd.DataFrame(res["predicted_targets"]), hide_index=True, use_container_width=True)
                st.caption("Targets are predicted by Tanimoto-weighted read-across over the query's nearest "
                           "neighbours' measured kinase activity (no clustering). Confidence = nearest-neighbour "
                           "Tanimoto. Selectivity = read-across S-score. Validation (out-of-sample): the known "
                           "target is in the top-10 ~50% of the time and the confidence flag is calibrated — "
                           "trust High/Moderate, treat Low/Outlier as a hint only.")
    elif go_btn:
        st.warning("Enter a SMILES string first.")

# ============================ VIRTUAL SCREENING ============================
with tab2:
    st.markdown("##### Screen a library")
    st.caption("Paste SMILES (one per line) or upload a CSV with a SMILES column. Each molecule is scored "
               "by read-across; rank by selectivity, filter by confidence, or rank by predicted activity "
               "against a chosen kinase. Up to 100,000 molecules (large libraries take ~1 min per 50k).")
    col_in, col_opt = st.columns([3, 2])
    with col_in:
        up = st.file_uploader("upload CSV (must have a SMILES column)", type=["csv"])
        txt = st.text_area("…or paste SMILES (one per line)", height=150,
                           placeholder="CCOc1cc2ncc(C#N)...\nCC(C)(C)c1cnc(...)o1\n...")
    with col_opt:
        target = st.selectbox("rank by predicted activity for a kinase (optional)",
                              ["— overall (most selective first) —"] + rx.genes)
        only_domain = st.checkbox("keep only in-domain hits (High/Moderate confidence)", value=False)
        run = st.button("Screen library", type="primary", use_container_width=True)

    if run:
        smis = []
        if up is not None:
            try:
                dfu = pd.read_csv(up)
                scol = next((c for c in dfu.columns if c.strip().lower() in
                             ("smiles", "smi", "canonical_smiles", "isomeric_smiles", "structure")), None)
                if scol is None:
                    for c in dfu.select_dtypes(include="object").columns:
                        s0 = dfu[c].dropna()
                        if len(s0) and Chem.MolFromSmiles(str(s0.iloc[0])) is not None:
                            scol = c; break
                if scol:
                    smis += [str(s) for s in dfu[scol].dropna()]
                else:
                    st.error("No SMILES column found in the CSV.")
            except Exception as e:
                st.error(f"Could not read CSV: {e}")
        smis += [s.strip() for s in txt.splitlines() if s.strip()]
        smis = smis[:100000]

        if not smis:
            st.warning("Provide some SMILES (paste or upload a CSV).")
        else:
            by_target = target != "— overall (most selective first) —"
            tgt = target if by_target else None
            rows, n_bad = [], 0
            N = len(smis); step = max(1, N // 100)
            prog = st.progress(0.0, text=f"screening {N:,} molecules…")
            for i, smi in enumerate(smis):
                row = rx.screen_row(smi, target_gene=tgt)   # fast vectorized path
                if row is None:
                    n_bad += 1
                else:
                    rows.append(row)
                if i % step == 0 or i == N - 1:
                    prog.progress((i + 1) / N, text=f"screening {N:,} molecules…  {i+1:,}/{N:,}")
            prog.empty()

            df = pd.DataFrame(rows)
            if only_domain and len(df):
                df = df[df.confidence.isin(["High", "Moderate"])]
            if by_target and f"pAct[{target}]" in df:
                df = df.sort_values(f"pAct[{target}]", ascending=False, na_position="last")
            elif len(df):
                df = df.sort_values("S_score", ascending=True, na_position="last")
            df = df.reset_index(drop=True); df.index += 1

            msg = f"screened {len(rows)} valid molecules"
            if n_bad: msg += f" ({n_bad} unparseable skipped)"
            if by_target: msg += f" · ranked by predicted activity at {target}"
            st.success(msg)
            st.dataframe(df, use_container_width=True)
            st.download_button("⬇ download results (CSV)", df.to_csv(index=False).encode(),
                               "huksa_rx_screen.csv", "text/csv")
            st.caption("Predicted pActivity > 6 ≈ sub-µM engagement. Trust **High/Moderate** confidence rows; "
                       "**Low/Outlier** are out-of-domain guesses. A blank kinase column means that kinase was "
                       "not measured among the molecule's nearest neighbours.")
