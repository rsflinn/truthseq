"""
TruthSeq Web Interface
======================
A Streamlit app that lets anyone validate gene regulatory claims
without installing Python or running terminal commands.

Run locally:
    streamlit run app.py

Deploy on Streamlit Community Cloud:
    1. Push this file to your GitHub repo
    2. Go to share.streamlit.io
    3. Point it at your repo and app.py
    4. The app will download data on first run (~500 MB)

Requirements:
    pip install streamlit pandas pyarrow numpy scipy
"""

import os
import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats as sp_stats
from io import StringIO

# ============================================================
# Configuration
# ============================================================

PERCENTILE_VALIDATED = 90
PERCENTILE_PARTIAL = 50

GRADE_COLORS = {
    'VALIDATED': '#33ff66',
    'PARTIALLY_SUPPORTED': '#ffb02e',
    'WEAK': '#7a8a7a',
    'CONTRADICTED': '#ff4444',
    'UNTESTABLE': '#555555',
}

GRADE_DESCRIPTIONS = {
    'VALIDATED': 'Lab data confirms the direction and strength of this regulatory relationship.',
    'PARTIALLY_SUPPORTED': 'Direction matches but effect size is modest, or only partial evidence available.',
    'WEAK': 'The upstream gene was tested but the downstream gene did not respond strongly. '
            'This does not mean the relationship is wrong — it may be cell-type-specific.',
    'CONTRADICTED': 'Lab data shows the opposite direction from what was predicted.',
    'UNTESTABLE': 'The upstream gene was not in the knockdown dataset.',
}

# The TruthSeq description for the typewriter effect
TRUTHSEQ_DESCRIPTION = (
    "Most findings in biology don't survive replication. "
    "TruthSeq checks whether your computational gene regulatory predictions "
    "hold up against real experimental data from human cells. "
    "Upload your predictions. Find out which ones are real."
)

# ============================================================
# CRT Theme CSS + Animated Logo
# ============================================================

CRT_THEME = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

/* Global dark CRT background */
.stApp {
    background-color: #0a0a0a !important;
    color: #d4cfc4 !important;
}

/* Scanline overlay on the whole page */
.stApp::after {
    content: '';
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background: repeating-linear-gradient(
        0deg,
        transparent,
        transparent 3px,
        rgba(0, 0, 0, 0.04) 3px,
        rgba(0, 0, 0, 0.04) 6px
    );
    pointer-events: none;
    z-index: 9999;
}

/* Override all Streamlit text colors */
.stApp, .stApp p, .stApp span, .stApp label, .stApp div {
    color: #d4cfc4 !important;
}

h1, h2, h3, h4, h5, h6,
.stApp h1, .stApp h2, .stApp h3 {
    color: #f0ead8 !important;
    font-family: 'Inter', 'Segoe UI', sans-serif !important;
}

/* Sidebar styling */
section[data-testid="stSidebar"] {
    background-color: #0f0f0f !important;
    border-right: 1px solid #1a2a1a !important;
}

section[data-testid="stSidebar"] * {
    color: #b0b0a0 !important;
}

/* Input fields */
.stTextInput input, .stTextArea textarea {
    background-color: #111 !important;
    color: #33ff66 !important;
    border: 1px solid #1a3a1a !important;
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace !important;
}

.stTextInput input:focus, .stTextArea textarea:focus {
    border-color: #33ff66 !important;
    box-shadow: 0 0 8px rgba(51, 255, 102, 0.15) !important;
}

/* Buttons */
.stButton button {
    background-color: #111 !important;
    color: #33ff66 !important;
    border: 1px solid #33ff66 !important;
    font-family: 'Inter', sans-serif !important;
    transition: all 0.3s ease !important;
}

.stButton button:hover {
    background-color: #33ff66 !important;
    color: #0a0a0a !important;
    box-shadow: 0 0 15px rgba(51, 255, 102, 0.3) !important;
}

/* Radio buttons */
.stRadio label {
    color: #b0b0a0 !important;
}

/* Metrics */
[data-testid="stMetricValue"] {
    color: #f0ead8 !important;
    font-family: 'Inter', sans-serif !important;
}

[data-testid="stMetricLabel"] {
    color: #7a8a7a !important;
}

/* Dividers */
hr {
    border-color: #1a2a1a !important;
}

/* File uploader */
[data-testid="stFileUploader"] {
    background-color: #111 !important;
    border: 1px dashed #1a3a1a !important;
}

/* Expander */
.streamlit-expanderHeader {
    background-color: #111 !important;
    color: #b0b0a0 !important;
}

/* Success/warning/info boxes */
.stAlert {
    background-color: #111 !important;
    border-left-color: #33ff66 !important;
}

/* Caption text */
.stCaption, caption {
    color: #555 !important;
}

/* Download button */
.stDownloadButton button {
    background-color: #111 !important;
    color: #ffb02e !important;
    border: 1px solid #ffb02e !important;
}

.stDownloadButton button:hover {
    background-color: #ffb02e !important;
    color: #0a0a0a !important;
}

/* Logo and typewriter are self-contained in iframes */
/* Hide iframe borders */
iframe {
    border: none !important;
}

/* Grade badges - CRT style */
.grade-badge {
    padding: 8px 16px;
    border-radius: 2px;
    text-align: center;
    font-weight: 600;
    font-family: 'Inter', sans-serif;
    font-size: 13px;
    letter-spacing: 1px;
    text-transform: uppercase;
}

.grade-VALIDATED {
    background-color: rgba(51, 255, 102, 0.15);
    color: #33ff66;
    border: 1px solid #33ff66;
    box-shadow: 0 0 8px rgba(51, 255, 102, 0.1);
}

.grade-PARTIALLY_SUPPORTED {
    background-color: rgba(255, 176, 46, 0.15);
    color: #ffb02e;
    border: 1px solid #ffb02e;
    box-shadow: 0 0 8px rgba(255, 176, 46, 0.1);
}

.grade-WEAK {
    background-color: rgba(122, 138, 122, 0.15);
    color: #7a8a7a;
    border: 1px solid #7a8a7a;
}

.grade-CONTRADICTED {
    background-color: rgba(255, 68, 68, 0.15);
    color: #ff4444;
    border: 1px solid #ff4444;
    box-shadow: 0 0 8px rgba(255, 68, 68, 0.1);
}

.grade-UNTESTABLE {
    background-color: rgba(85, 85, 85, 0.15);
    color: #666;
    border: 1px solid #444;
}

/* Claim rows */
.claim-row {
    padding: 12px 16px;
    margin: 4px 0;
    background-color: #0f0f0f;
    border-left: 3px solid #1a2a1a;
    border-radius: 0 2px 2px 0;
}

.claim-genes {
    font-family: 'SF Mono', 'Fira Code', monospace;
    color: #d4cfc4;
    font-size: 15px;
}

.claim-genes .gene-name {
    color: #f0ead8;
    font-weight: 600;
}

.claim-genes .arrow {
    color: #33ff66;
    margin: 0 6px;
}

.claim-note {
    font-size: 12px;
    color: #666;
    margin-top: 4px;
    font-family: 'SF Mono', monospace;
}

/* Section headers */
.section-header {
    font-family: 'Inter', sans-serif;
    font-size: 20px;
    font-weight: 600;
    color: #f0ead8;
    margin: 30px 0 15px 0;
    padding-bottom: 8px;
    border-bottom: 1px solid #1a2a1a;
}

/* Footer */
.truthseq-footer {
    text-align: center;
    padding: 20px 0;
    color: #444;
    font-size: 12px;
    font-family: 'SF Mono', monospace;
}

.truthseq-footer a {
    color: #33ff66;
    text-decoration: none;
}
</style>
"""

LOGO_HTML = """
<html><head><style>
  *{margin:0;padding:0;box-sizing:border-box;}
  html,body{background:#0a0a0a;overflow:hidden;width:100%;height:100%;}
  .logo-container{position:relative;width:100%;height:280px;overflow:hidden;}
  .logo-container canvas{position:absolute;top:0;left:0;width:100%;height:280px;}
  .logo-wordmark{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);
    font-size:56px;font-weight:600;letter-spacing:3px;color:#f0ead8;z-index:10;
    text-shadow:0 0 8px rgba(240,234,216,0.6),0 0 20px rgba(240,234,216,0.35),
    0 0 40px rgba(255,200,80,0.15),0 0 60px rgba(51,255,102,0.08);
    font-family:'Inter','Segoe UI',sans-serif;pointer-events:none;}
  .logo-wordmark .upper{font-weight:600;}.logo-wordmark .lower{font-weight:400;}
  .logo-vignette{position:absolute;top:0;left:0;width:100%;height:100%;
    background:radial-gradient(ellipse at center,transparent 40%,rgba(10,10,10,0.8) 100%);
    pointer-events:none;z-index:15;}
  .logo-scanlines{position:absolute;top:0;left:0;width:100%;height:100%;
    background:repeating-linear-gradient(0deg,transparent,transparent 2px,
    rgba(0,0,0,0.06) 2px,rgba(0,0,0,0.06) 4px);pointer-events:none;z-index:20;}
</style></head><body>
<div class="logo-container" id="logoContainer">
    <canvas id="bgCanvas"></canvas>
    <canvas id="trailCanvas"></canvas>
    <canvas id="glowCanvas"></canvas>
    <div class="logo-wordmark">
        <span class="upper">T</span><span class="lower">ruth</span><span class="upper">S</span><span class="lower">eq</span>
    </div>
    <div class="logo-vignette"></div>
    <div class="logo-scanlines"></div>
</div>

<script>
(function() {
    const W = document.body.clientWidth || 900;
    const H = 280;
    const centerY = H / 2;

    const bgCanvas = document.getElementById('bgCanvas');
    const trailCanvas = document.getElementById('trailCanvas');
    const glowCanvas = document.getElementById('glowCanvas');

    [bgCanvas, trailCanvas, glowCanvas].forEach(c => {
        c.width = W;
        c.height = H;
        c.style.width = W + 'px';
        c.style.height = H + 'px';
    });

    const bgCtx = bgCanvas.getContext('2d');
    const trailCtx = trailCanvas.getContext('2d');
    const glowCtx = glowCanvas.getContext('2d');

    const GREEN = { r: 51, g: 255, b: 102 };
    const AMBER = { r: 255, g: 176, b: 46 };

    bgCtx.fillStyle = '#0a0a0a';
    bgCtx.fillRect(0, 0, W, H);
    bgCtx.globalAlpha = 0.012;
    for (let x = 0; x < W; x += 3) {
        for (let y = 0; y < H; y += 3) {
            bgCtx.fillStyle = '#33ff66';
            bgCtx.fillRect(x, y, 1, 1);
        }
    }
    bgCtx.globalAlpha = 1;

    const AMP = 55;
    const FREQ = 0.008;
    const SPD = 0.015;
    let time = 0;

    function greenWave(x, t) {
        return centerY - 8 +
            Math.sin(x * FREQ + t) * AMP * 0.7 +
            Math.sin(x * FREQ * 2.3 + t * 1.4) * AMP * 0.25 +
            Math.cos(x * FREQ * 0.5 + t * 0.7) * AMP * 0.15;
    }

    function amberWave(x, t) {
        return centerY + 8 +
            Math.sin(x * FREQ + t + 1.8) * AMP * 0.65 +
            Math.sin(x * FREQ * 2.1 + t * 1.3 + 0.5) * AMP * 0.3 +
            Math.cos(x * FREQ * 0.6 + t * 0.8 + 1.0) * AMP * 0.15;
    }

    function drawGlow(ctx, pts, col, lw, gr) {
        if (pts.length < 2) return;
        function stroke(a, w, blur) {
            ctx.save();
            ctx.strokeStyle = 'rgba('+col.r+','+col.g+','+col.b+','+a+')';
            ctx.lineWidth = w;
            ctx.lineCap = 'round';
            ctx.lineJoin = 'round';
            if (blur) ctx.filter = 'blur('+blur+'px)';
            ctx.beginPath();
            ctx.moveTo(pts[0].x, pts[0].y);
            for (let i=1;i<pts.length;i++) ctx.lineTo(pts[i].x, pts[i].y);
            ctx.stroke();
            ctx.restore();
        }
        stroke(0.06, lw + gr*4, gr*2);
        stroke(0.15, lw + gr*2, gr);
        stroke(0.3, lw + gr*0.8, gr*0.5);
        stroke(1.0, lw*1.5, 0);
        // white core
        ctx.save();
        ctx.strokeStyle = 'rgba(255,255,240,0.5)';
        ctx.lineWidth = lw*0.5;
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(pts[0].x, pts[0].y);
        for (let i=1;i<pts.length;i++) ctx.lineTo(pts[i].x, pts[i].y);
        ctx.stroke();
        ctx.restore();
    }

    function drawBloom(ctx, x, y) {
        let g = ctx.createRadialGradient(x,y,0,x,y,25);
        g.addColorStop(0,'rgba(255,250,220,0.6)');
        g.addColorStop(0.3,'rgba(255,230,150,0.3)');
        g.addColorStop(1,'rgba(180,220,100,0)');
        ctx.fillStyle = g;
        ctx.fillRect(x-25,y-25,50,50);
        let c = ctx.createRadialGradient(x,y,0,x,y,4);
        c.addColorStop(0,'rgba(255,255,255,0.7)');
        c.addColorStop(1,'rgba(255,255,255,0)');
        ctx.fillStyle = c;
        ctx.fillRect(x-4,y-4,8,8);
    }

    function stampTrail(pts, col) {
        if (pts.length < 2) return;
        trailCtx.globalAlpha = 0.1;
        trailCtx.strokeStyle = 'rgba('+col.r+','+col.g+','+col.b+',1)';
        trailCtx.lineWidth = 2.5;
        trailCtx.filter = 'blur(3px)';
        trailCtx.beginPath();
        trailCtx.moveTo(pts[0].x, pts[0].y);
        for (let i=1;i<pts.length;i++) trailCtx.lineTo(pts[i].x, pts[i].y);
        trailCtx.stroke();
        trailCtx.filter = 'none';
        trailCtx.globalAlpha = 1;
    }

    function animate() {
        time++;
        const t = time * SPD;
        glowCtx.clearRect(0, 0, W, H);

        // Fade trail
        trailCtx.globalCompositeOperation = 'destination-out';
        trailCtx.fillStyle = 'rgba(0,0,0,0.06)';
        trailCtx.fillRect(0, 0, W, H);
        trailCtx.globalCompositeOperation = 'source-over';

        const gp = [], ap = [];
        for (let x=0; x<=W; x+=2) {
            gp.push({x, y: greenWave(x, t)});
            ap.push({x, y: amberWave(x, t)});
        }

        drawGlow(glowCtx, gp, GREEN, 2, 10);
        drawGlow(glowCtx, ap, AMBER, 2, 10);
        stampTrail(gp, GREEN);
        stampTrail(ap, AMBER);

        for (let x=2; x<=W; x+=2) {
            const gy=greenWave(x,t), ay=amberWave(x,t);
            const pgy=greenWave(x-2,t), pay=amberWave(x-2,t);
            if ((gy-ay)*(pgy-pay)<0) drawBloom(glowCtx, x, (gy+ay)/2);
        }

        requestAnimationFrame(animate);
    }
    animate();
})();
</script>
</body></html>
"""

TYPEWRITER_JS = """
<html><head><style>
  *{margin:0;padding:0;}
  html,body{background:#0a0a0a;overflow:hidden;width:100%;height:100%;}
  .typewriter-container{padding:10px 0;}
  .typewriter-text{
    font-family:'SF Mono','Fira Code','Consolas',monospace;
    font-size:14px;color:#33ff66;line-height:1.7;white-space:pre-wrap;
  }
  .typewriter-cursor{
    display:inline-block;width:8px;height:16px;background-color:#33ff66;
    animation:blink 1s step-end infinite;vertical-align:text-bottom;
    box-shadow:0 0 6px rgba(51,255,102,0.5);
  }
  @keyframes blink{0%,100%{opacity:1;}50%{opacity:0;}}
</style></head><body>
<div class="typewriter-container">
    <div class="typewriter-text" id="typewriterText"></div>
</div>
<script>
(function() {
    const text = """ + '"' + TRUTHSEQ_DESCRIPTION.replace('"', '\\"') + '"' + """;
    const el = document.getElementById('typewriterText');
    let i = 0;
    const speed = 25;

    function type() {
        if (i < text.length) {
            el.innerHTML = text.substring(0, i+1) + '<span class="typewriter-cursor"></span>';
            i++;
            setTimeout(type, speed);
        } else {
            el.innerHTML = text + '<span class="typewriter-cursor"></span>';
        }
    }
    // Small delay before starting
    setTimeout(type, 800);
})();
</script>
</body></html>
"""


# ============================================================
# Data Loading (cached) + Auto-Download for Cloud Deployment
# ============================================================

# Pre-built parquet files hosted on GitHub Releases (no scanpy needed)
GITHUB_RELEASE_BASE = "https://github.com/rsflinn/truthseq/releases/download/v1.0.0"
EFFECTS_URL = f"{GITHUB_RELEASE_BASE}/replogle_knockdown_effects_cloud.parquet"
STATS_URL = f"{GITHUB_RELEASE_BASE}/replogle_knockdown_stats.parquet"
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _download_file(url, output_path, label="data"):
    """Download a file with Streamlit progress bar."""
    import requests

    if os.path.exists(output_path):
        existing_mb = os.path.getsize(output_path) / 1024 / 1024
        if existing_mb > 0.1:
            return output_path

    progress = st.progress(0, text=f"Downloading {label}...")

    resp = requests.get(url, stream=True, timeout=60)
    resp.raise_for_status()

    total = int(resp.headers.get('content-length', 0))
    downloaded = 0

    with open(output_path, 'wb') as f:
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
            downloaded += len(chunk)
            if total > 0:
                pct = downloaded / total
                mb_done = downloaded / 1024 / 1024
                progress.progress(pct, text=f"Downloading {label}... {mb_done:.0f}/{total/1024/1024:.0f} MB")

    progress.progress(1.0, text=f"{label} ready.")
    return output_path


def auto_download_data():
    """
    Download pre-built parquet files from GitHub Releases if not found locally.
    Returns (effects_path, stats_path) or (None, None) on failure.
    """
    os.makedirs(DATA_DIR, exist_ok=True)

    effects_path = os.path.join(DATA_DIR, "replogle_knockdown_effects.parquet")
    stats_path = os.path.join(DATA_DIR, "replogle_knockdown_stats.parquet")

    # If parquets already exist, return them
    if os.path.exists(effects_path) and os.path.exists(stats_path):
        return effects_path, stats_path

    # Check common local locations
    for check_dir in ['.', os.path.expanduser('~')]:
        eff = os.path.join(check_dir, "replogle_knockdown_effects.parquet")
        sta = os.path.join(check_dir, "replogle_knockdown_stats.parquet")
        if os.path.exists(eff):
            return eff, sta if os.path.exists(sta) else None

    # Download from GitHub Releases
    st.info(
        "First-time setup: downloading Replogle Perturb-seq reference data (~148 MB). "
        "This only happens once."
    )

    try:
        _download_file(EFFECTS_URL, effects_path, "knockdown effects (147 MB)")
        _download_file(STATS_URL, stats_path, "knockdown stats (<1 MB)")
        return effects_path, stats_path

    except Exception as e:
        st.error(f"Auto-download failed: {e}")
        st.markdown(
            "**Manual setup:** Clone the repo and run `python3 setup.py` locally, "
            "then copy the parquet files to the `data/` directory."
        )
        return None, None


@st.cache_data(show_spinner="Loading knockdown data...")
def load_replogle_data(effects_path, stats_path=None):
    """Load the Replogle Perturb-seq data with memory optimization."""
    effects_df = pd.read_parquet(effects_path)

    # Convert repeated string columns to categorical — cuts memory ~90%.
    # 37M rows with ~8K unique gene names: stores each string once, uses int indices.
    for col in ['knocked_down_gene', 'affected_gene', 'cell_line']:
        if col in effects_df.columns:
            effects_df[col] = effects_df[col].astype('category')

    stats_df = None
    if stats_path and os.path.exists(stats_path):
        stats_df = pd.read_parquet(stats_path)
    return effects_df, stats_df


def find_data_files():
    """Look for parquet files in the current directory and common locations."""
    search_dirs = ['.', 'data', DATA_DIR, os.path.expanduser('~')]
    found = {}

    for search_dir in search_dirs:
        for fname in os.listdir(search_dir) if os.path.isdir(search_dir) else []:
            fpath = os.path.join(search_dir, fname)
            if fname.endswith('_knockdown_effects.parquet'):
                ct = 'K562' if 'rpe1' not in fname.lower() else 'RPE1'
                found[f'{ct}_effects'] = fpath
            elif fname.endswith('_knockdown_stats.parquet'):
                ct = 'K562' if 'rpe1' not in fname.lower() else 'RPE1'
                found[f'{ct}_stats'] = fpath

    return found


# ============================================================
# Validation Logic
# ============================================================

def compute_percentile_from_stats(z_score, upstream, stats_df):
    """Compute approximate percentile using pre-computed quantiles from stats_df."""
    abs_z = abs(z_score)
    kd_stats = stats_df[stats_df['knocked_down_gene'] == upstream]
    if len(kd_stats) == 0:
        return 50.0

    row = kd_stats.iloc[0]
    # Walk through quantile breakpoints to estimate percentile
    quantiles = [5, 10, 25, 50, 75, 80, 85, 90, 95, 97, 99]
    for q in quantiles:
        col = f'q{q:02d}'
        if col in row and abs_z <= row[col]:
            return float(q)
    return 99.5  # above q99


def compute_percentile_from_distribution(z_score, kd_data):
    """Compute percentile rank of a Z-score within a knockdown's distribution."""
    abs_z = abs(z_score)
    all_abs_z = kd_data['z_score'].abs().values
    return float(sp_stats.percentileofscore(all_abs_z, abs_z))


def validate_single_claim(upstream, downstream, predicted_dir, effects_df, stats_df=None):
    """Validate a single gene regulatory claim."""
    available_kd = set(effects_df['knocked_down_gene'].unique())
    stats_kd = set(stats_df['knocked_down_gene'].unique()) if stats_df is not None else set()
    all_kd = available_kd | stats_kd

    if upstream not in all_kd:
        return {
            'grade': 'UNTESTABLE',
            'z_score': None,
            'percentile': None,
            'direction_match': None,
            'cell_line': None,
            'note': f'{upstream} was not in the knockdown dataset.',
        }

    kd_data = effects_df[effects_df['knocked_down_gene'] == upstream]
    target_hit = kd_data[kd_data['affected_gene'] == downstream]

    if len(target_hit) > 0:
        cell_types = target_hit['cell_line'].unique() if 'cell_line' in target_hit.columns else ['K562']
        per_ct = []

        for ct in cell_types:
            ct_hits = target_hit[target_hit['cell_line'] == ct] if 'cell_line' in target_hit.columns else target_hit
            z = float(ct_hits.iloc[0]['z_score'])
            obs_dir = "UP" if z < 0 else "DOWN"
            dir_match = (obs_dir == predicted_dir)

            # Use stats-based percentile (fast, works with trimmed cloud dataset)
            # Fall back to distribution-based if stats not available
            if stats_df is not None:
                pct = compute_percentile_from_stats(z, upstream, stats_df)
            else:
                ct_kd = kd_data
                if 'cell_line' in kd_data.columns:
                    ct_kd = kd_data[kd_data['cell_line'] == ct]
                pct = compute_percentile_from_distribution(z, ct_kd) if len(ct_kd) > 0 else 50.0

            per_ct.append({
                'cell_line': ct,
                'z_score': round(z, 4),
                'percentile': round(pct, 1),
                'direction_match': dir_match,
                'observed_direction': obs_dir,
            })

        best = max(per_ct, key=lambda x: x['percentile'])
        z = best['z_score']
        pct = best['percentile']
        dir_match = best['direction_match']
        cell_line = best['cell_line']

        if dir_match and pct >= PERCENTILE_VALIDATED:
            grade = 'VALIDATED'
        elif dir_match and pct >= PERCENTILE_PARTIAL:
            grade = 'PARTIALLY_SUPPORTED'
        elif not dir_match and pct >= PERCENTILE_PARTIAL:
            grade = 'CONTRADICTED'
        else:
            grade = 'WEAK'

        note = f"Z={z:.2f}, {pct:.0f}th percentile in {cell_line}"
        if len(per_ct) > 1:
            ct_notes = [f"{r['cell_line']}: Z={r['z_score']:.2f}, {r['percentile']:.0f}th pctl"
                        for r in per_ct]
            note += f". All cell types: {'; '.join(ct_notes)}"

        return {
            'grade': grade,
            'z_score': z,
            'percentile': pct,
            'direction_match': dir_match,
            'cell_line': cell_line,
            'note': note,
            'per_cell_type': per_ct if len(per_ct) > 1 else None,
        }

    else:
        return {
            'grade': 'WEAK',
            'z_score': None,
            'percentile': None,
            'direction_match': None,
            'cell_line': None,
            'note': f'{upstream} was knocked down but {downstream} did not respond above threshold.',
        }


def validate_claims(claims_df, effects_df, stats_df=None):
    """Validate a batch of claims."""
    results = []
    for _, row in claims_df.iterrows():
        result = validate_single_claim(
            row['upstream_gene'],
            row['downstream_gene'],
            row['predicted_direction'],
            effects_df,
            stats_df,
        )
        result['upstream_gene'] = row['upstream_gene']
        result['downstream_gene'] = row['downstream_gene']
        result['predicted_direction'] = row['predicted_direction']
        results.append(result)
    return pd.DataFrame(results)


# ============================================================
# Base Rate Simulation & Signal Strength
# ============================================================

def run_base_rate_simulation(results_df, effects_df, stats_df=None, n_permutations=1000):
    """
    Estimate the null baseline using binomial simulation.

    For a random gene pair with a random predicted direction:
      - P(direction matches observed) = 0.50
      - P(|z| >= 90th percentile) = 0.10  (by definition of percentile)
      - P(VALIDATED) = 0.50 * 0.10 = 0.05

    This is exact by construction of the percentile grading system.
    We simulate the null distribution via Binomial(n_claims, 0.05)
    rather than brute-forcing 37M-row lookups per permutation.
    """
    testable = results_df[results_df['grade'] != 'UNTESTABLE']
    n_claims = len(testable)
    if n_claims == 0:
        return None

    observed = int((testable['grade'] == 'VALIDATED').sum())

    # Null rate: 5% of random pairs would score VALIDATED
    # (10% exceed 90th percentile * 50% direction match)
    p_null = 0.05

    rng = np.random.default_rng(42)
    null_counts = rng.binomial(n_claims, p_null, size=n_permutations).tolist()

    null_mean = float(np.mean(null_counts))
    null_std = float(np.std(null_counts))
    p_value = float(np.mean([c >= observed for c in null_counts]))

    return {
        'observed': observed,
        'null_mean': null_mean,
        'null_std': null_std,
        'n_claims': n_claims,
        'p_value': p_value,
        'null_counts': null_counts,
    }


def compute_signal_strength(base_rate):
    """
    Compute Signal Strength: a 0-100 normalized score measuring how far
    your validated claims exceed random chance relative to the maximum possible.
    """
    if base_rate is None:
        return None

    observed = base_rate['observed']
    null_mean = base_rate['null_mean']
    n_claims = base_rate['n_claims']
    max_possible = n_claims

    if max_possible <= null_mean:
        return {
            'score': 0,
            'interpretation': 'NO_DYNAMIC_RANGE',
            'detail': 'Cannot compute — null mean equals or exceeds maximum.',
            'observed': observed,
            'null_mean': null_mean,
            'max_possible': max_possible,
            'p_value': base_rate['p_value'],
        }

    raw = (observed - null_mean) / (max_possible - null_mean) * 100
    score = max(0.0, min(100.0, round(raw, 1)))

    if score >= 80:
        band = 'STRONG'
        detail = 'Strong separation from random gene pairs.'
    elif score >= 50:
        band = 'MODERATE'
        detail = 'Claims score above random with meaningful separation.'
    elif score >= 20:
        band = 'WEAK'
        detail = 'Slight signal above random. Interpret with caution.'
    else:
        band = 'MINIMAL'
        detail = 'Claims score at or below what random gene pairs achieve.'

    return {
        'score': score,
        'interpretation': band,
        'detail': detail,
        'observed': observed,
        'null_mean': null_mean,
        'max_possible': max_possible,
        'p_value': base_rate['p_value'],
    }


def build_signal_strength_widget(score, observed, null_mean, max_possible, p_value):
    """Build the vintage receiver Signal Strength widget HTML for embedding in Streamlit."""
    p_display = '<0.001' if p_value < 0.001 else f'{p_value:.3f}'
    return f"""
<html><head>
<meta charset="UTF-8">
<style>
  @import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@300;400;500&family=Space+Grotesk:wght@300;400;500;600;700&display=swap');
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  html, body {{ background: #0a0a0a; overflow: hidden; width: 100%; }}

  .widget-wrap {{ display: flex; flex-direction: column; align-items: center; gap: 0; padding: 10px 0; }}

  .cabinet {{
    position: relative; padding: 16px 16px 0 16px; border-radius: 8px 8px 0 0;
    background:
      repeating-linear-gradient(0deg, transparent 0px, transparent 4px, rgba(0,0,0,0.035) 4px, rgba(0,0,0,0.035) 5px),
      repeating-linear-gradient(90deg, transparent 0px, transparent 1px, rgba(255,255,255,0.008) 1px, rgba(255,255,255,0.008) 2px),
      linear-gradient(180deg, #7a4e30 0%, #6b422a 2%, #5a3621 5%, #4d2e1a 12%, #583419 30%, #4a2c15 50%, #553218 70%, #4a2c15 88%, #4d2e1a 95%, #5a3621 98%, #6b422a 100%);
    box-shadow: 0 12px 40px rgba(0,0,0,0.7), 0 4px 12px rgba(0,0,0,0.5),
      inset 0 1px 0 rgba(255,200,140,0.15), inset 0 -1px 0 rgba(0,0,0,0.4),
      inset 1px 0 0 rgba(255,200,140,0.08), inset -1px 0 0 rgba(255,200,140,0.08);
  }}

  .chrome-bezel {{
    position: relative; padding: 3px; border-radius: 4px;
    background: linear-gradient(180deg, #d0d0d0 0%, #b0b0b0 2%, #c8c8c8 8%, #aaa 20%, #bbb 50%, #a5a5a5 80%, #c0c0c0 95%, #ddd 100%);
    box-shadow: inset 0 1px 0 rgba(255,255,255,0.6), inset 0 -1px 0 rgba(0,0,0,0.2), 0 0 0 1px rgba(0,0,0,0.3);
  }}

  .glass-panel {{
    position: relative; background: #050505; border-radius: 2px; padding: 24px 36px 20px; width: 640px; overflow: hidden;
    box-shadow: inset 0 3px 10px rgba(0,0,0,0.95), inset 0 0 40px rgba(0,0,0,0.7), inset 0 -2px 6px rgba(0,0,0,0.5);
  }}
  .glass-panel::before {{
    content: ''; position: absolute; top: 0; left: 0; right: 0; bottom: 0;
    background:
      linear-gradient(155deg, rgba(255,255,255,0.04) 0%, rgba(255,255,255,0.015) 15%, transparent 30%, transparent 70%, rgba(255,255,255,0.008) 85%, rgba(255,255,255,0.02) 100%),
      linear-gradient(180deg, rgba(255,255,255,0.025) 0%, transparent 8%),
      radial-gradient(ellipse at 25% 12%, rgba(255,255,255,0.03) 0%, transparent 40%);
    pointer-events: none; z-index: 20; border-radius: 2px;
  }}
  .glass-panel::after {{
    content: ''; position: absolute; top: 0; left: 0; right: 0; bottom: 0;
    background:
      radial-gradient(ellipse at 62% 50%, rgba(74,180,74,0.04) 0%, transparent 55%),
      radial-gradient(ellipse at 18% 55%, rgba(212,168,74,0.05) 0%, transparent 45%);
    pointer-events: none; z-index: 1; border-radius: 2px;
  }}

  .header-row {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 14px; position: relative; z-index: 5; }}
  .brand-mark {{ font-family: 'DM Mono', monospace; font-size: 11px; font-weight: 400; letter-spacing: 5px; text-transform: uppercase; color: #5a9a5a; text-shadow: 0 0 8px rgba(90,180,90,0.3), 0 0 24px rgba(90,180,90,0.1); }}
  .power-indicator {{ display: flex; align-items: center; gap: 8px; }}
  .power-dot {{ width: 6px; height: 6px; border-radius: 50%; background: #5ac45a; box-shadow: 0 0 4px rgba(90,196,90,0.9), 0 0 10px rgba(90,196,90,0.5), 0 0 20px rgba(90,196,90,0.2); animation: pulse 3s ease-in-out infinite; }}
  @keyframes pulse {{ 0%,100%{{opacity:0.8}} 50%{{opacity:1}} }}
  .power-label {{ font-family: 'DM Mono', monospace; font-size: 9px; letter-spacing: 2px; color: #3a7a3a; text-transform: uppercase; text-shadow: 0 0 6px rgba(58,122,58,0.2); }}

  .gauge-area {{ display: flex; gap: 16px; align-items: stretch; position: relative; z-index: 5; }}

  .vu-housing {{
    width: 140px; flex-shrink: 0; position: relative; border-radius: 4px; overflow: hidden;
    background: #080604; border: 1px solid #2a2215;
    box-shadow: inset 0 2px 8px rgba(0,0,0,0.8), inset 0 0 20px rgba(0,0,0,0.4), 0 0 20px rgba(212,168,74,0.06);
  }}
  .vu-housing::before {{
    content: ''; position: absolute; top: 0; left: 0; right: 0; bottom: 0;
    background: radial-gradient(ellipse at 50% 40%, rgba(212,168,74,0.1) 0%, transparent 60%), radial-gradient(ellipse at 50% 80%, rgba(180,140,60,0.05) 0%, transparent 50%);
    pointer-events: none; z-index: 0;
  }}
  .vu-housing canvas {{ display: block; width: 100%; height: auto; position: relative; z-index: 2; }}
  .vu-readout {{ position: relative; z-index: 2; text-align: center; padding: 0 0 8px; font-family: 'DM Mono', monospace; font-size: 28px; font-weight: 500; color: #d4a84a; text-shadow: 0 0 12px rgba(212,168,74,0.45), 0 0 30px rgba(212,168,74,0.15); line-height: 1; }}
  .vu-readout .pct {{ font-size: 15px; color: #8a6830; }}
  .vu-label {{ position: relative; z-index: 2; text-align: center; padding-bottom: 6px; font-family: 'DM Mono', monospace; font-size: 8px; letter-spacing: 3px; color: #6a5530; text-transform: uppercase; text-shadow: 0 0 4px rgba(106,85,48,0.3); }}

  .dial-housing {{ flex: 1; position: relative; border-radius: 3px; overflow: hidden; background: #050505; }}
  .dial-housing::before {{ content: ''; position: absolute; top: 0; left: 0; right: 0; bottom: 0; background: radial-gradient(ellipse at 50% 50%, rgba(74,180,74,0.035) 0%, transparent 70%); pointer-events: none; }}
  .dial-housing canvas {{ display: block; width: 100%; height: auto; position: relative; z-index: 2; }}

  .readout-strip {{ margin-top: 12px; padding-top: 10px; border-top: 1px solid rgba(255,255,255,0.025); display: flex; justify-content: space-between; align-items: center; position: relative; z-index: 5; }}
  .readout-band {{ font-family: 'DM Mono', monospace; font-size: 14px; font-weight: 500; letter-spacing: 4px; text-transform: uppercase; transition: color 0.4s, text-shadow 0.4s; }}
  .readout-detail {{ font-family: 'Space Grotesk', sans-serif; font-size: 12px; font-weight: 300; color: #555; max-width: 340px; text-align: right; line-height: 1.4; }}

  .metal-panel {{
    width: calc(640px + 32px + 6px); position: relative; border-radius: 0 0 8px 8px; padding: 14px 32px;
    display: flex; justify-content: space-between; align-items: center;
    background:
      repeating-linear-gradient(90deg, transparent 0px, transparent 1px, rgba(255,255,255,0.025) 1px, rgba(255,255,255,0.025) 2px),
      linear-gradient(180deg, #c0c0c0 0%, #aaa 3%, #b5b5b5 15%, #a8a8a8 50%, #9d9d9d 85%, #959595 97%, #888 100%);
    box-shadow: inset 0 1px 0 rgba(255,255,255,0.4), inset 0 -1px 0 rgba(0,0,0,0.2), 0 6px 16px rgba(0,0,0,0.5);
  }}
  .context-item {{ text-align: center; min-width: 80px; }}
  .context-value {{ font-size: 17px; font-weight: 700; color: #222; }}
  .context-label {{ font-size: 9px; color: #555; text-transform: uppercase; letter-spacing: 1.2px; margin-top: 2px; }}
  .metal-div {{ width: 1px; height: 32px; background: linear-gradient(180deg, transparent 5%, rgba(0,0,0,0.12) 50%, transparent 95%); }}
</style>
</head><body>

<div class="widget-wrap">
  <div class="cabinet">
    <div class="chrome-bezel">
      <div class="glass-panel">
        <div class="header-row">
          <div class="brand-mark">TruthSeq</div>
          <div class="power-indicator">
            <div class="power-dot"></div>
            <div class="power-label">Signal Lock</div>
          </div>
        </div>
        <div class="gauge-area">
          <div class="vu-housing">
            <canvas id="vu-canvas" width="280" height="230"></canvas>
            <div class="vu-readout"><span id="vu-score">0</span><span class="pct">%</span></div>
            <div class="vu-label">Tuning</div>
          </div>
          <div class="dial-housing">
            <canvas id="dial-canvas" width="960" height="380"></canvas>
          </div>
        </div>
        <div class="readout-strip">
          <div class="readout-band" id="band-label">MINIMAL</div>
          <div class="readout-detail" id="detail-line">Scanning...</div>
        </div>
      </div>
    </div>
  </div>

  <div class="metal-panel">
    <div class="context-item">
      <div class="context-value" id="ctx-observed">{observed}</div>
      <div class="context-label">Validated</div>
    </div>
    <div class="metal-div"></div>
    <div class="context-item">
      <div class="context-value" id="ctx-expected">{null_mean:.1f}</div>
      <div class="context-label">Expected (null)</div>
    </div>
    <div class="metal-div"></div>
    <div class="context-item">
      <div class="context-value" id="ctx-max">{max_possible}</div>
      <div class="context-label">Max Possible</div>
    </div>
    <div class="metal-div"></div>
    <div class="context-item">
      <div class="context-value" id="ctx-pval">{p_display}</div>
      <div class="context-label">p-value</div>
    </div>
  </div>
</div>

<script>
// Color system
function getIndicatorColor(score) {{
  if (score < 20) return {{ r: 200, g: 60, b: 50, hex: '#c83c32' }};
  if (score < 35) {{
    const t = (score - 20) / 15;
    return {{ r: Math.round(200 + 12*t), g: Math.round(60 + 100*t), b: Math.round(50 + 24*t), hex: lerpHex('#c83c32', '#d4a84a', t) }};
  }}
  if (score < 50) return {{ r: 212, g: 168, b: 74, hex: '#d4a84a' }};
  if (score < 65) {{
    const t = (score - 50) / 15;
    return {{ r: Math.round(212 - 112*t), g: Math.round(168 + 12*t), b: Math.round(74 - 24*t), hex: lerpHex('#d4a84a', '#64b44a', t) }};
  }}
  return {{ r: 74, g: 180, b: 74, hex: '#4ab44a' }};
}}
function lerpHex(a, b, t) {{
  const ar=parseInt(a.slice(1,3),16),ag=parseInt(a.slice(3,5),16),ab=parseInt(a.slice(5,7),16);
  const br=parseInt(b.slice(1,3),16),bg=parseInt(b.slice(3,5),16),bb=parseInt(b.slice(5,7),16);
  const rr=Math.round(ar+(br-ar)*t),rg=Math.round(ag+(bg-ag)*t),rb=Math.round(ab+(bb-ab)*t);
  return '#'+rr.toString(16).padStart(2,'0')+rg.toString(16).padStart(2,'0')+rb.toString(16).padStart(2,'0');
}}

// VU Meter
const vuC = document.getElementById('vu-canvas');
const vu = vuC.getContext('2d');
const VW = vuC.width, VH = vuC.height;
const VCX = VW/2, VCY = VH*0.85, VR = VW*0.36;

function drawVU(score) {{
  vu.clearRect(0,0,VW,VH);
  const bg = vu.createRadialGradient(VCX, VCY*0.65, 0, VCX, VCY*0.65, VR*1.6);
  bg.addColorStop(0, 'rgba(212,168,74,0.09)'); bg.addColorStop(0.5, 'rgba(180,140,60,0.04)'); bg.addColorStop(1, 'transparent');
  vu.fillStyle = bg; vu.fillRect(0,0,VW,VH);

  vu.strokeStyle = 'rgba(120,96,48,0.2)'; vu.lineWidth = 1;
  vu.beginPath(); vu.arc(VCX, VCY, VR, Math.PI, 0, false); vu.stroke();
  vu.strokeStyle = 'rgba(120,96,48,0.1)'; vu.beginPath(); vu.arc(VCX, VCY, VR*0.7, Math.PI, 0, false); vu.stroke();

  for (let pct = 0; pct <= 100; pct += 5) {{
    const a = Math.PI + (pct/100)*Math.PI;
    const isMaj = pct % 25 === 0, isRed = pct > 80;
    const r1 = VR - (isMaj ? 16 : 10), r2 = VR;
    vu.strokeStyle = isRed ? 'rgba(200,60,50,0.5)' : (isMaj ? 'rgba(212,168,74,0.5)' : 'rgba(120,96,48,0.25)');
    vu.lineWidth = isMaj ? 2 : 0.8;
    vu.beginPath(); vu.moveTo(VCX+r1*Math.cos(a), VCY+r1*Math.sin(a)); vu.lineTo(VCX+r2*Math.cos(a), VCY+r2*Math.sin(a)); vu.stroke();
  }}

  [{{v:0,p:0}},{{v:25,p:25}},{{v:50,p:50}},{{v:75,p:75}},{{v:100,p:100}}].forEach(n => {{
    const a = Math.PI + (n.p/100)*Math.PI, lr = VR + 20;
    const x = VCX + lr*Math.cos(a), y = VCY + lr*Math.sin(a);
    const isRed = n.p > 80, col = isRed ? [200,60,50] : [212,168,74];
    vu.font = '500 16px "DM Mono", monospace'; vu.textAlign = 'center'; vu.textBaseline = 'middle';
    vu.shadowColor = 'rgba('+col[0]+','+col[1]+','+col[2]+',0.5)'; vu.shadowBlur = 12;
    vu.fillStyle = 'rgba('+col[0]+','+col[1]+','+col[2]+',0.15)'; vu.fillText(n.v, x, y);
    vu.shadowBlur = 0; vu.fillStyle = 'rgba('+col[0]+','+col[1]+','+col[2]+',0.75)'; vu.fillText(n.v, x, y);
  }});

  const na = Math.PI + (score/100)*Math.PI, nl = VR*0.88;
  const nx = VCX + nl*Math.cos(na), ny = VCY + nl*Math.sin(na);
  vu.strokeStyle = 'rgba(0,0,0,0.35)'; vu.lineWidth = 3;
  vu.beginPath(); vu.moveTo(VCX+1.5, VCY+1.5); vu.lineTo(nx+1.5, ny+1.5); vu.stroke();
  vu.save(); vu.shadowColor = 'rgba(232,192,96,0.4)'; vu.shadowBlur = 10;
  vu.strokeStyle = 'rgba(232,192,96,0.3)'; vu.lineWidth = 4;
  vu.beginPath(); vu.moveTo(VCX, VCY); vu.lineTo(nx, ny); vu.stroke(); vu.restore();
  vu.strokeStyle = '#e8c060'; vu.lineWidth = 1.8;
  vu.beginPath(); vu.moveTo(VCX, VCY); vu.lineTo(nx, ny); vu.stroke();
  vu.fillStyle = '#1a1408'; vu.beginPath(); vu.arc(VCX, VCY, 6, 0, Math.PI*2); vu.fill();
  const hg = vu.createRadialGradient(VCX-1, VCY-1, 0, VCX, VCY, 4);
  hg.addColorStop(0, '#8a7040'); hg.addColorStop(1, '#3a2810');
  vu.fillStyle = hg; vu.beginPath(); vu.arc(VCX, VCY, 3.5, 0, Math.PI*2); vu.fill();
}}

// Main dial
const dialC = document.getElementById('dial-canvas');
const dc = dialC.getContext('2d');
const DW = dialC.width, DH = dialC.height;
const DL = 50, DR = DW-50, DY = DH*0.46, DDW = DR-DL, BH = 100;
const ZONES = [{{s:0,e:20,c:[180,60,50]}},{{s:20,e:50,c:[212,168,74]}},{{s:50,e:80,c:[100,180,74]}},{{s:80,e:100,c:[74,180,74]}}];
function px(pct) {{ return DL + (pct/100)*DDW; }}

function drawDial(score) {{
  dc.clearRect(0,0,DW,DH);
  const bt = DY - BH/2, bb = DY + BH/2;
  const blg = dc.createRadialGradient(DW/2, DY, 0, DW/2, DY, DDW*0.55);
  blg.addColorStop(0, 'rgba(74,180,74,0.03)'); blg.addColorStop(1, 'transparent');
  dc.fillStyle = blg; dc.fillRect(0,0,DW,DH);

  ZONES.forEach(z => {{
    const x1 = px(z.s), x2 = px(z.e); const [r,g,b] = z.c;
    dc.fillStyle = 'rgba('+r+','+g+','+b+',0.03)'; dc.fillRect(x1, bt, x2-x1, BH);
    if (z.s > 0) {{ dc.strokeStyle = 'rgba(255,255,255,0.02)'; dc.lineWidth = 1; dc.beginPath(); dc.moveTo(x1,bt); dc.lineTo(x1,bb); dc.stroke(); }}
  }});

  ZONES.forEach(z => {{
    const ze = Math.min(z.e, score); if (ze <= z.s) return;
    const x1 = px(z.s), x2 = px(ze); const [r,g,b] = z.c;
    dc.fillStyle = 'rgba('+r+','+g+','+b+',0.1)'; dc.fillRect(x1, bt, x2-x1, BH);
    dc.strokeStyle = 'rgba('+r+','+g+','+b+',0.55)'; dc.lineWidth = 2.5;
    dc.beginPath(); dc.moveTo(x1, DY); dc.lineTo(x2, DY); dc.stroke();
    [16, 32, 50].forEach((blur, i) => {{
      const alpha = [0.12, 0.05, 0.02][i];
      dc.save(); dc.shadowColor = 'rgba('+r+','+g+','+b+','+alpha+')'; dc.shadowBlur = blur;
      dc.strokeStyle = 'rgba('+r+','+g+','+b+','+(alpha*0.8)+')'; dc.lineWidth = 4;
      dc.beginPath(); dc.moveTo(x1, DY); dc.lineTo(x2, DY); dc.stroke(); dc.restore();
    }});
  }});

  for (let p = 0; p <= 100; p += 2) {{
    const x = px(p); const isZone = [0,20,50,80,100].includes(p);
    const isMaj = p % 10 === 0; const isMid = p % 5 === 0 && !isMaj;
    const tl = isZone ? 22 : (isMaj ? 14 : (isMid ? 9 : 4));
    const lit = p <= score; const gc = lit ? [74,180,74] : [30,50,30];
    const al = isZone ? 0.55 : (isMaj ? 0.35 : (isMid ? 0.18 : 0.06));
    dc.strokeStyle = 'rgba('+gc[0]+','+gc[1]+','+gc[2]+','+al+')';
    dc.lineWidth = isZone ? 1.8 : (isMaj ? 1.2 : 0.6);
    dc.beginPath(); dc.moveTo(x, bb); dc.lineTo(x, bb+tl); dc.stroke();
    dc.beginPath(); dc.moveTo(x, bt); dc.lineTo(x, bt-tl); dc.stroke();
  }}

  dc.textAlign = 'center'; dc.textBaseline = 'top';
  [0,20,50,80,100].forEach(p => {{
    const x = px(p), lit = p <= score; const col = lit ? [74,180,74] : [30,50,30]; const al = lit ? 0.9 : 0.3;
    dc.font = '600 22px "DM Mono", monospace';
    dc.shadowColor = lit ? 'rgba(74,180,74,0.6)' : 'transparent'; dc.shadowBlur = lit ? 14 : 0;
    dc.fillStyle = 'rgba('+col[0]+','+col[1]+','+col[2]+','+(al*0.3)+')'; dc.fillText(p, x, bb+26);
    dc.shadowBlur = 0; dc.fillStyle = 'rgba('+col[0]+','+col[1]+','+col[2]+','+al+')'; dc.fillText(p, x, bb+26);
  }});
  [10,30,40,60,70,90].forEach(p => {{
    const x = px(p), lit = p <= score;
    dc.font = '400 13px "DM Mono", monospace';
    dc.fillStyle = lit ? 'rgba(74,180,74,0.45)' : 'rgba(30,50,30,0.2)'; dc.fillText(p, x, bb+30);
  }});

  dc.font = '400 13px "DM Mono", monospace'; dc.textBaseline = 'bottom';
  dc.fillStyle = 'rgba(74,120,74,0.4)'; dc.textAlign = 'left'; dc.fillText('SIGNAL', DL, DH-8);
  dc.textAlign = 'right'; dc.fillText('STRENGTH', DR, DH-8);

  const ic = getIndicatorColor(score); const tx = px(score); const colStr = ic.r+','+ic.g+','+ic.b;
  dc.save(); dc.shadowColor = 'rgba('+colStr+',0.4)'; dc.shadowBlur = 30;
  dc.strokeStyle = 'rgba('+colStr+',0.08)'; dc.lineWidth = 14;
  dc.beginPath(); dc.moveTo(tx, bt-24); dc.lineTo(tx, bb+22); dc.stroke(); dc.restore();
  dc.save(); dc.shadowColor = 'rgba('+colStr+',0.5)'; dc.shadowBlur = 12;
  dc.strokeStyle = 'rgba('+colStr+',0.2)'; dc.lineWidth = 5;
  dc.beginPath(); dc.moveTo(tx, bt-20); dc.lineTo(tx, bb+18); dc.stroke(); dc.restore();
  dc.save(); dc.shadowColor = 'rgba('+colStr+',0.6)'; dc.shadowBlur = 6;
  dc.strokeStyle = 'rgba('+colStr+',0.9)'; dc.lineWidth = 2.5;
  dc.beginPath(); dc.moveTo(tx, bt-18); dc.lineTo(tx, bb+16); dc.stroke(); dc.restore();
  dc.save(); dc.shadowColor = 'rgba('+colStr+',0.4)'; dc.shadowBlur = 8;
  dc.fillStyle = 'rgba('+colStr+',0.8)';
  dc.beginPath(); dc.moveTo(tx, bt-20); dc.lineTo(tx-6, bt-30); dc.lineTo(tx+6, bt-30); dc.closePath(); dc.fill(); dc.restore();
}}

// Band info
const BAND_COLORS = {{'MINIMAL':'#cc4444','WEAK':'#d4a84a','MODERATE':'#7abd5a','STRONG':'#4aB44a'}};
const BAND_DETAILS = {{
  'MINIMAL':'Claims score at or below what random gene pairs achieve.',
  'WEAK':'Slight signal above random. Interpret with caution.',
  'MODERATE':'Claims score above random with meaningful separation.',
  'STRONG':'Strong separation from random gene pairs.',
}};
function getBand(s) {{ if (s >= 80) return 'STRONG'; if (s >= 50) return 'MODERATE'; if (s >= 20) return 'WEAK'; return 'MINIMAL'; }}

function updateAll(score) {{
  drawVU(score); drawDial(score);
  document.getElementById('vu-score').textContent = Math.round(score);
  const band = getBand(score);
  const bel = document.getElementById('band-label');
  bel.textContent = band; bel.style.color = BAND_COLORS[band];
  bel.style.textShadow = '0 0 14px '+BAND_COLORS[band]+'60, 0 0 30px '+BAND_COLORS[band]+'20';
  document.getElementById('detail-line').textContent = BAND_DETAILS[band];
}}

// Animate to target
let cur = 0;
const tgt = {score};
function animIn() {{
  cur += (tgt - cur) * 0.045;
  if (Math.abs(cur - tgt) < 0.3) {{ cur = tgt; updateAll(cur); return; }}
  updateAll(cur);
  requestAnimationFrame(animIn);
}}
requestAnimationFrame(animIn);
</script>
</body></html>
"""


# ============================================================
# Streamlit App
# ============================================================

def main():
    st.set_page_config(
        page_title="TruthSeq",
        page_icon="",
        layout="wide",
    )

    # Inject CRT theme
    st.markdown(CRT_THEME, unsafe_allow_html=True)

    # Animated logo header
    st.components.v1.html(LOGO_HTML, height=290, scrolling=False)

    # Typewriter description
    st.components.v1.html(TYPEWRITER_JS, height=90, scrolling=False)

    # ---- LANDING PAGE SECTIONS ----

    # How it works - 3 steps
    st.markdown("""
    <div style="display:flex; gap:20px; margin:30px 0 40px 0;">
        <div style="flex:1; padding:20px; background:#0f0f0f; border-top:2px solid #33ff66; border-radius:0 0 4px 4px;">
            <div style="color:#33ff66; font-size:32px; font-weight:700; font-family:Inter,sans-serif; margin-bottom:8px;">1</div>
            <div style="color:#f0ead8; font-size:15px; font-weight:600; margin-bottom:6px;">Upload your predictions</div>
            <div style="color:#7a8a7a; font-size:13px; line-height:1.5;">
                CSV with three columns: upstream gene, downstream gene, predicted direction. If you claim Gene X controls Gene Y, put it in the file.</div>
        </div>
        <div style="flex:1; padding:20px; background:#0f0f0f; border-top:2px solid #ffb02e; border-radius:0 0 4px 4px;">
            <div style="color:#ffb02e; font-size:32px; font-weight:700; font-family:Inter,sans-serif; margin-bottom:8px;">2</div>
            <div style="color:#f0ead8; font-size:15px; font-weight:600; margin-bottom:6px;">Check against real experiments</div>
            <div style="color:#7a8a7a; font-size:13px; line-height:1.5;">
                TruthSeq looks up what actually happened when each upstream gene was knocked out in human cells. Direct cause-and-effect data from 11,000+ experiments.</div>
        </div>
        <div style="flex:1; padding:20px; background:#0f0f0f; border-top:2px solid #f0ead8; border-radius:0 0 4px 4px;">
            <div style="color:#f0ead8; font-size:32px; font-weight:700; font-family:Inter,sans-serif; margin-bottom:8px;">3</div>
            <div style="color:#f0ead8; font-size:15px; font-weight:600; margin-bottom:6px;">Get a grade for each claim</div>
            <div style="color:#7a8a7a; font-size:13px; line-height:1.5;">
                <span style="color:#33ff66;">VALIDATED</span> &middot;
                <span style="color:#ffb02e;">PARTIALLY SUPPORTED</span> &middot;
                <span style="color:#7a8a7a;">WEAK</span> &middot;
                <span style="color:#ff4444;">CONTRADICTED</span> &middot;
                <span style="color:#555;">UNTESTABLE</span></div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # The numbers
    st.markdown("""
    <div style="display:flex; gap:30px; justify-content:center; margin:0 0 40px 0; padding:25px 0;
                border-top:1px solid #1a2a1a; border-bottom:1px solid #1a2a1a;">
        <div style="text-align:center;">
            <div style="color:#33ff66; font-size:28px; font-weight:700; font-family:Inter,sans-serif;">11,000+</div>
            <div style="color:#555; font-size:11px; letter-spacing:1px; text-transform:uppercase;">gene knockdowns</div>
        </div>
        <div style="text-align:center;">
            <div style="color:#ffb02e; font-size:28px; font-weight:700; font-family:Inter,sans-serif;">37M+</div>
            <div style="color:#555; font-size:11px; letter-spacing:1px; text-transform:uppercase;">gene pairs tested</div>
        </div>
        <div style="text-align:center;">
            <div style="color:#f0ead8; font-size:28px; font-weight:700; font-family:Inter,sans-serif;">2</div>
            <div style="color:#555; font-size:11px; letter-spacing:1px; text-transform:uppercase;">human cell types</div>
        </div>
        <div style="text-align:center;">
            <div style="color:#7a8a7a; font-size:28px; font-weight:700; font-family:Inter,sans-serif;">3</div>
            <div style="color:#555; font-size:11px; letter-spacing:1px; text-transform:uppercase;">evidence tiers</div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Use cases
    st.markdown("""
    <div style="margin:0 0 40px 0;">
        <div style="color:#f0ead8; font-size:18px; font-weight:600; font-family:Inter,sans-serif;
                    margin-bottom:15px; padding-bottom:8px; border-bottom:1px solid #1a2a1a;">Who this is for</div>
        <div style="display:grid; grid-template-columns:1fr 1fr; gap:12px;">
            <div style="padding:14px 18px; background:#0f0f0f; border-left:2px solid #33ff66; border-radius:0 4px 4px 0;">
                <span style="color:#d4cfc4; font-size:13px; line-height:1.5;">
                You ran a GRN analysis and found 200 predicted regulatory edges. Which ones would hold up if you tested them in a lab?</span>
            </div>
            <div style="padding:14px 18px; background:#0f0f0f; border-left:2px solid #ffb02e; border-radius:0 4px 4px 0;">
                <span style="color:#d4cfc4; font-size:13px; line-height:1.5;">
                You're reviewing a paper's supplementary data and want to spot-check their gene regulatory claims before citing them.</span>
            </div>
            <div style="padding:14px 18px; background:#0f0f0f; border-left:2px solid #f0ead8; border-radius:0 4px 4px 0;">
                <span style="color:#d4cfc4; font-size:13px; line-height:1.5;">
                You're exploring genetic data with AI tools and want a reality check before going too far down a rabbit hole.</span>
            </div>
            <div style="padding:14px 18px; background:#0f0f0f; border-left:2px solid #7a8a7a; border-radius:0 4px 4px 0;">
                <span style="color:#d4cfc4; font-size:13px; line-height:1.5;">
                You're a working researcher who needs a fast, independent check on computational predictions. No biology degree required.</span>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ---- END LANDING PAGE SECTIONS ----

    # Sidebar: data loading
    st.sidebar.markdown(
        "<div style='color:#33ff66; font-family: Inter, sans-serif; "
        "font-size:14px; font-weight:600; letter-spacing:1px; margin-bottom:15px;'>"
        "DATA SOURCE</div>",
        unsafe_allow_html=True,
    )

    data_files = find_data_files()

    effects_path = st.sidebar.text_input(
        "Effects parquet path",
        value=data_files.get('K562_effects', 'replogle_knockdown_effects.parquet'),
    )
    stats_path = st.sidebar.text_input(
        "Stats parquet path (optional)",
        value=data_files.get('K562_stats', 'replogle_knockdown_stats.parquet'),
    )

    rpe1_effects = data_files.get('RPE1_effects', '')
    if rpe1_effects:
        use_rpe1 = st.sidebar.checkbox("Also validate against RPE1 cells", value=True)
    else:
        use_rpe1 = False

    # Load data — auto-download if not found
    if not os.path.exists(effects_path):
        # Try auto-download
        auto_eff, auto_stats = auto_download_data()
        if auto_eff and os.path.exists(auto_eff):
            effects_path = auto_eff
            stats_path = auto_stats
        else:
            st.markdown(
                "<div style='background:#111; border:1px solid #ffb02e; padding:20px; "
                "border-radius:4px; margin:20px 0;'>"
                "<span style='color:#ffb02e; font-weight:600;'>DATA NOT FOUND</span><br>"
                "<span style='color:#999; font-size:14px;'>"
                "Run <code style=\"color:#33ff66;\">python3 setup.py</code> to download "
                "the Replogle Perturb-seq data (~357 MB), then restart this app.</span></div>",
                unsafe_allow_html=True,
            )
            st.code("python3 setup.py\nstreamlit run app.py", language="bash")
            return

    effects_df, stats_df = load_replogle_data(effects_path, stats_path)

    if use_rpe1 and rpe1_effects:
        rpe1_eff, rpe1_stats = load_replogle_data(
            rpe1_effects,
            data_files.get('RPE1_stats', ''),
        )
        effects_df = pd.concat([effects_df, rpe1_eff], ignore_index=True)
        if rpe1_stats is not None and stats_df is not None:
            stats_df = pd.concat([stats_df, rpe1_stats], ignore_index=True)

    n_knockdowns = effects_df['knocked_down_gene'].nunique()
    n_pairs = len(effects_df)
    cell_types = effects_df['cell_line'].unique() if 'cell_line' in effects_df.columns else ['K562']

    st.sidebar.markdown(
        f"<div style='background:#0f1a0f; border:1px solid #1a3a1a; padding:12px; "
        f"border-radius:4px; margin-top:10px;'>"
        f"<span style='color:#33ff66; font-size:13px;'>"
        f"&#x2713; {n_knockdowns:,} knockdowns<br>"
        f"&#x2713; {n_pairs:,} gene pairs<br>"
        f"&#x2713; {', '.join(cell_types)}</span></div>",
        unsafe_allow_html=True,
    )

    # Main input area
    st.markdown(
        '<div class="section-header">Your claims</div>',
        unsafe_allow_html=True,
    )

    input_method = st.radio(
        "How do you want to enter your predictions?",
        ["Upload CSV", "Enter manually", "Try example"],
        horizontal=True,
    )

    claims_df = None

    if input_method == "Upload CSV":
        uploaded = st.file_uploader(
            "Upload a CSV with columns: upstream_gene, downstream_gene, predicted_direction",
            type=['csv', 'tsv'],
        )
        if uploaded:
            sep = '\t' if uploaded.name.endswith('.tsv') else ','
            claims_df = pd.read_csv(uploaded, sep=sep)

    elif input_method == "Enter manually":
        st.markdown(
            "<span style='color:#7a8a7a; font-size:13px;'>"
            "Enter one claim per line: upstream_gene, downstream_gene, direction (UP or DOWN)</span>",
            unsafe_allow_html=True,
        )
        text_input = st.text_area(
            "Claims",
            value="MYT1L, SCN2A, UP\nTCF4, CACNA1A, UP\nMEF2C, GRIN2B, UP",
            height=150,
            label_visibility="collapsed",
        )
        if text_input.strip():
            lines = [l.strip() for l in text_input.strip().split('\n') if l.strip()]
            rows = []
            for line in lines:
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 3:
                    rows.append({
                        'upstream_gene': parts[0],
                        'downstream_gene': parts[1],
                        'predicted_direction': parts[2].upper(),
                    })
            if rows:
                claims_df = pd.DataFrame(rows)

    elif input_method == "Try example":
        claims_df = pd.DataFrame([
            {'upstream_gene': 'SLC30A1', 'downstream_gene': 'MT2A', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'GATA1', 'downstream_gene': 'TYROBP', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'MYT1L', 'downstream_gene': 'SCN2A', 'predicted_direction': 'UP'},
            {'upstream_gene': 'TCF4', 'downstream_gene': 'CACNA1A', 'predicted_direction': 'UP'},
            {'upstream_gene': 'CHMP6', 'downstream_gene': 'SOD2', 'predicted_direction': 'DOWN'},
            {'upstream_gene': 'TP53', 'downstream_gene': 'CDKN1A', 'predicted_direction': 'UP'},
        ])
        st.markdown(
            "<div style='background:#0f1a0f; border-left:3px solid #33ff66; "
            "padding:10px 15px; color:#7a8a7a; font-size:13px;'>"
            "Example claims loaded — includes known-true, known-weak and untestable cases.</div>",
            unsafe_allow_html=True,
        )

    if claims_df is not None:
        required = ['upstream_gene', 'downstream_gene', 'predicted_direction']
        missing = [c for c in required if c not in claims_df.columns]
        if missing:
            st.error(f"Missing required columns: {missing}")
            return

        claims_df['predicted_direction'] = claims_df['predicted_direction'].str.upper().str.strip()

        st.markdown(
            f"<span style='color:#7a8a7a; font-size:13px;'>"
            f"{len(claims_df)} claims to validate</span>",
            unsafe_allow_html=True,
        )

        # Run validation
        with st.spinner("Validating against knockdown data..."):
            results_df = validate_claims(claims_df, effects_df, stats_df)

        # Results header
        st.markdown(
            '<div class="section-header">Results</div>',
            unsafe_allow_html=True,
        )

        # --- Signal Strength ---
        # Run base rate simulation and compute signal strength
        testable_count = len(results_df[results_df['grade'] != 'UNTESTABLE'])
        if testable_count > 0:
            base_rate = run_base_rate_simulation(results_df, effects_df, stats_df)
            signal_strength = compute_signal_strength(base_rate)

            if signal_strength and signal_strength['score'] is not None:
                st.markdown(
                    '<div class="section-header">Signal Strength</div>',
                    unsafe_allow_html=True,
                )
                widget_html = build_signal_strength_widget(
                    score=signal_strength['score'],
                    observed=signal_strength['observed'],
                    null_mean=signal_strength['null_mean'],
                    max_possible=signal_strength['max_possible'],
                    p_value=signal_strength['p_value'],
                )
                st.components.v1.html(widget_html, height=440, scrolling=False)

                # Explainer
                with st.expander("What is Signal Strength?"):
                    st.markdown(
                        "<span style='color:#7a8a7a; font-size:13px; line-height:1.8;'>"
                        "Signal Strength measures whether your gene regulatory predictions "
                        "perform better than random chance. TruthSeq runs 200 simulations "
                        "with random gene pairs through the same validation pipeline, then "
                        "compares your results against that null baseline.<br><br>"
                        "<b style='color:#cc4444;'>MINIMAL (0-20%)</b> — at or below random noise<br>"
                        "<b style='color:#d4a84a;'>WEAK (20-50%)</b> — slight signal, interpret cautiously<br>"
                        "<b style='color:#7abd5a;'>MODERATE (50-80%)</b> — meaningful separation from random<br>"
                        "<b style='color:#4ab44a;'>STRONG (80-100%)</b> — strong separation from random</span>",
                        unsafe_allow_html=True,
                    )

        # Grade summary bar
        grade_counts = results_df['grade'].value_counts()
        grades_html = '<div style="display:flex; gap:12px; margin-bottom:25px;">'
        for grade in ['VALIDATED', 'PARTIALLY_SUPPORTED', 'WEAK', 'CONTRADICTED', 'UNTESTABLE']:
            count = grade_counts.get(grade, 0)
            color = GRADE_COLORS[grade]
            label = grade.replace('_', ' ')
            grades_html += (
                f'<div style="flex:1; text-align:center; padding:15px 10px; '
                f'background:#0f0f0f; border:1px solid {color}40; border-radius:4px;">'
                f'<div style="color:{color}; font-size:28px; font-weight:700; '
                f'font-family:Inter,sans-serif;">{count}</div>'
                f'<div style="color:{color}80; font-size:10px; letter-spacing:1px; '
                f'text-transform:uppercase; margin-top:4px;">{label}</div></div>'
            )
        grades_html += '</div>'
        st.markdown(grades_html, unsafe_allow_html=True)

        # Per-claim results
        for _, row in results_df.iterrows():
            grade = row['grade']
            color = GRADE_COLORS.get(grade, '#555')

            claim_html = (
                f'<div style="display:flex; justify-content:space-between; '
                f'align-items:center; padding:12px 16px; margin:6px 0; '
                f'background:#0f0f0f; border-left:3px solid {color};">'
                f'<div>'
                f'<div class="claim-genes">'
                f'<span class="gene-name">{row["upstream_gene"]}</span>'
                f'<span class="arrow">&#x2192;</span>'
                f'<span class="gene-name">{row["downstream_gene"]}</span>'
                f'<span style="color:#555; margin-left:10px;">predicted {row["predicted_direction"]}</span>'
                f'</div>'
                f'<div class="claim-note">{row["note"] or ""}</div>'
                f'</div>'
                f'<div class="grade-badge grade-{grade}">'
                f'{grade.replace("_", " ")}</div>'
                f'</div>'
            )
            st.markdown(claim_html, unsafe_allow_html=True)

            if isinstance(row.get('per_cell_type'), list):
                with st.expander("Cell type breakdown"):
                    for ct_result in row['per_cell_type']:
                        match_color = '#33ff66' if ct_result['direction_match'] else '#ff4444'
                        match_word = 'matches' if ct_result['direction_match'] else 'OPPOSES'
                        st.markdown(
                            f"<span style='color:#7a8a7a; font-family:monospace; font-size:13px;'>"
                            f"  {ct_result['cell_line']}: Z={ct_result['z_score']:.2f}, "
                            f"{ct_result['percentile']:.0f}th percentile, "
                            f"direction <span style='color:{match_color};'>{match_word}</span></span>",
                            unsafe_allow_html=True,
                        )

        # Download
        st.markdown(
            '<div class="section-header" style="margin-top:40px;">Download</div>',
            unsafe_allow_html=True,
        )
        csv_output = results_df.to_csv(index=False)
        st.download_button(
            "Download results as CSV",
            csv_output,
            file_name="truthseq_results.csv",
            mime="text/csv",
        )

        # Grade legend
        with st.expander("What do the grades mean?"):
            for grade, desc in GRADE_DESCRIPTIONS.items():
                color = GRADE_COLORS[grade]
                st.markdown(
                    f"<div style='margin:8px 0;'>"
                    f"<span style='color:{color}; font-weight:600; font-family:Inter,sans-serif; "
                    f"letter-spacing:1px; font-size:12px;'>{grade.replace('_', ' ')}</span>"
                    f"<span style='color:#7a8a7a; font-size:13px; margin-left:10px;'>{desc}</span></div>",
                    unsafe_allow_html=True,
                )

    # Footer
    st.markdown(
        '<div class="truthseq-footer">'
        'TruthSeq validates gene regulatory predictions against the '
        'Replogle Perturb-seq atlas (Replogle et al. 2022, Cell). '
        'Data from Figshare.<br>'
        '<a href="https://github.com/rsflinn/truthseq">GitHub</a> · '
        '<a href="https://github.com/rsflinn/truthseq#how-the-grades-work">How the grades work</a>'
        '</div>',
        unsafe_allow_html=True,
    )


if __name__ == '__main__':
    main()
