import streamlit as st
import matplotlib.pyplot as plt
import io
import re
import urllib.request
import json
import pandas as pd
from stmol import showmol
import py3Dmol

# CONFIG & ULTRA-SLEEK MODERN PRESTIGE THEME
st.set_page_config(page_title="Bioinformatics Suite", page_icon="🧬", layout="wide")

st.markdown("""
    <style>
    /* 1. Efek Tekstur Kertas Jeruk Eksklusif pada Background Utama */
    .stApp {
        background-color: #F8FAFC !important;
        background-image: url("data:image/svg+xml,%3Csvg viewBox='0 0 200 200' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='noiseFilter'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.85' numOctaves='3' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23noiseFilter)' opacity='0.045'/%3E%3C/svg%3E") !important;
    }
    
    /* 2. Reset Gap Atas & Padding Layout */
    .block-container {
        padding-top: 2rem !important;
        padding-bottom: 2rem !important;
        padding-left: 3.5rem !important;
        padding-right: 3.5rem !important;
    }
    
    /* 3. Desain Typografi Header Gradasi Premium */
    .dashboard-title {
        font-family: 'Inter', -apple-system, sans-serif;
        background: linear-gradient(135deg, #4F46E5 0%, #0EA5E9 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 2.6rem !important;
        font-weight: 800 !important;
        letter-spacing: -0.5px;
        margin-bottom: 6px !important;
    }
    
    .dashboard-subtitle {
        font-family: 'Inter', sans-serif;
        color: #475569 !important;
        font-size: 1.1rem;
        line-height: 1.5;
        margin-bottom: 18px !important;
    }
    .dashboard-subtitle b { color: #0EA5E9 !important; font-weight: 600; }
    
    /* Description Wrapper Elegan Menyatu dengan Tekstur */
    .dashboard-desc {
        font-family: 'Inter', sans-serif;
        color: #1E293B !important;
        font-size: 1rem;
        line-height: 1.6;
        text-align: justify;
        background: rgba(255, 255, 255, 0.6);
        backdrop-filter: blur(4px);
        padding: 16px 20px;
        border-radius: 12px;
        border: 1px solid rgba(0, 0, 0, 0.06);
        margin-bottom: 25px;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.02);
    }

    /* 4. Penyelarasan Komponen Sidebar & Logo IPB */
    [data-testid="stSidebar"] {
        background-color: #0F172A !important; /* Sidebar tetap dark premium untuk kontras tinggi */
        border-right: 1px solid rgba(255, 255, 255, 0.05);
    }
    
    [data-testid="stSidebar"] [data-testid="stImage"] {
        background-color: #FFFFFF !important;
        padding: 12px !important;
        border-radius: 12px !important;
        box-shadow: 0px 10px 25px rgba(99, 102, 241, 0.15) !important;
        margin-bottom: 25px !important;
    }
    
    .sidebar-menu-title {
        color: #818CF8 !important;
        font-size: 0.9rem !important;
        font-weight: 700 !important;
        letter-spacing: 1px;
        margin-bottom: 15px !important;
        text-transform: uppercase;
    }
    

    /* 5. Desain Modern Elemen Tombol (Button Hover Effect) */
    div.stButton > button {
        background: linear-gradient(135deg, #4F46E5 0%, #3730A3 100%) !important;
        color: #FFFFFF !important;
        font-weight: 600 !important;
        letter-spacing: 0.3px !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 12px 24px !important;
        transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1) !important;
        box-shadow: 0px 4px 12px rgba(79, 70, 229, 0.2) !important;
    }
    div.stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0px 8px 20px rgba(79, 70, 229, 0.35) !important;
        filter: brightness(1.1);
    }
    
    /* 6. FIX TOTAL: MEWARNAI PUTIH TEKS LABEL UPLOADER & VISIBILITAS INDIKATOR */
    [data-testid="stSidebar"] [data-testid="stWidgetLabel"] p,
    [data-testid="stSidebar"] label p,
    .stFileUploader label p {
        color: #FFFFFF !important;
        font-weight: 600 !important;
        font-size: 0.95rem !important;
    }
    
    [data-testid="stFileUploader"] section {
        border: 2px dashed rgba(79, 70, 229, 0.4) !important;
        border-radius: 12px !important;
        background: rgba(255, 255, 255, 0.85) !important;
        padding: 15px !important;
    }
    [data-testid="stFileUploader"] section div {
        color: #0F172A !important;
        font-weight: 500 !important;
    }
    [data-testid="stFileUploader"] section small {
        color: #334155 !important;
        font-weight: 600 !important;
    }
    [data-testid="stFileUploaderDropzone"] {
        color: #0F172A !important;
    }
    [data-testid="stFileUploader"] button {
        background-color: #4F46E5 !important;
        color: white !important;
    }
    
    hr {
        border: 0;
        height: 1px;
        background: linear-gradient(to right, rgba(0,0,0,0), rgba(0,0,0,0.1), rgba(0,0,0,0));
        margin: 25px 0 !important;
    }
    </style>
""", unsafe_allow_html=True)

# CORE BIOINFORMATICS FUNCTIONS
def parse_fasta(file_name, file_content):
    records = []
    current_header = None
    current_seq = []
    for line in file_content.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if current_header:
                records.append({"file": file_name, "header": current_header, "sequence": "".join(current_seq)})
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line.upper())
    if current_header:
        records.append({"file": file_name, "header": current_header, "sequence": "".join(current_seq)})
    return records

def parse_fastq(file_name, file_content):
    records = []
    lines = file_content.splitlines()
    for i in range(0, len(lines), 4):
        if i + 1 < len(lines):
            header = lines[i].strip()[1:]
            seq = lines[i+1].strip().upper()
            if header and seq:
                records.append({"file": file_name, "header": header, "sequence": seq})
    return records

def calculate_nucleotide_frequency(sequence):
    freq_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    for base in sequence:
        if base in freq_dict: freq_dict[base] += 1
        else: freq_dict[base] = freq_dict.get(base, 0) + 1
    return freq_dict

def calculate_gc_content(sequence):
    if not sequence: return 0.0
    return ((sequence.count('G') + sequence.count('C')) / len(sequence)) * 100

def fetch_pdb_by_nucleotide_blast(sequence):
    try:
        clean_seq = re.sub(r'[^ATCGU]', '', sequence.upper()).replace('U', 'T')
        if len(clean_seq) < 15: return None
        if len(clean_seq) > 600: clean_seq = clean_seq[:600]
            
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        query_payload = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 100.0,
                    "identity_cutoff": 0.05,
                    "target": "pdb_dna_sequence",
                    "value": clean_seq
                }
            },
            "return_type": "entry"
        }
        req_data = json.dumps(query_payload).encode('utf-8')
        req = urllib.request.Request(url, data=req_data, headers={'Content-Type': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as response:
            res_data = json.loads(response.read().decode('utf-8'))
            if "result_set" in res_data and len(res_data["result_set"]) > 0:
                return str(res_data["result_set"][0]["identifier"]).upper()
    except Exception:
        pass
    return None

# SIDEBAR CONTROLS
st.sidebar.image("logo ipb horizontal.png", use_container_width=True)
st.sidebar.markdown("<p class='sidebar-menu-title'> Control Panel</p>", unsafe_allow_html=True)
uploaded_files = st.sidebar.file_uploader(
    "Unggah file sekuens (FASTA/FASTQ)", 
    type=["fasta", "fastq", "fa", "fq", "txt", "fna"],
    accept_multiple_files=True,
    key="bio_pipeline_uploader"
)

# HEADER MAIN CONTENT
st.markdown("""
    <div>
        <h1 class="dashboard-title">Integrated Bioinformatics Analytics Suite</h1>
        <p class="dashboard-subtitle">
            <b>Bioinformatics Study Program </b> | Mata Kuliah Struktur Data &nbsp;&bull;&nbsp; <b>IPB University</b>
        </p>
        <p class="dashboard-desc">
            Otomasi pipeline analisis GC content, manajemen komputasi list, pencatatan frekuensi kamus, ekspor berkas integrasi CSV, serta pemodelan spasial visualisasi interaktif 3D makromolekul secara komprehensif.
        </p>
    </div>
""", unsafe_allow_html=True)

st.markdown("<hr>", unsafe_allow_html=True)

raw_sequence_list = []

if uploaded_files:
    for f in uploaded_files:
        content = io.StringIO(f.getvalue().decode("utf-8")).read()
        filename = f.name.lower()
        is_fastq = filename.endswith("fastq") or filename.endswith("fq") or content.startswith("@")
        
        file_records = parse_fastq(f.name, content) if is_fastq else parse_fasta(f.name, content)
        raw_sequence_list.extend(file_records)

if raw_sequence_list:
    processed_pipeline_list = []
    
    for item in raw_sequence_list:
        seq = item["sequence"]
        freq = calculate_nucleotide_frequency(seq)
        gc_val = calculate_gc_content(seq)
        
        processed_pipeline_list.append({
            "File Sumber": item["file"],
            "Header ID": item["header"],
            "Panjang Rantai (bp)": len(seq),
            "GC Content (%)": round(gc_val, 2),
            "Frekuensi A": freq.get('A', 0),
            "Frekuensi T": freq.get('T', 0),
            "Frekuensi C": freq.get('C', 0),
            "Frekuensi G": freq.get('G', 0),
            "raw_sequence": seq
        })
        
    processed_pipeline_list.sort(key=lambda x: x["GC Content (%)"], reverse=True)
    df_all = pd.DataFrame(processed_pipeline_list)
    
    # MULTI-COLUMN DESIGN
    col_left, col_right = st.columns([1.1, 0.9], gap="large")
    
    with col_left:
        st.markdown("<h3 style='color: #4F46E5; font-size: 1.4rem; font-weight:700;'> Pipeline Analysis Output</h3>", unsafe_allow_html=True)
        st.markdown("Seluruh data sekuens dari berkas otomatis diurutkan berdasarkan GC Content secara menurun.")
        
        df_plot = df_all.head(10)
        
        fig, ax = plt.subplots(figsize=(10, 4.2))
        fig.patch.set_facecolor('none')  
        ax.set_facecolor('none')
        
        colors = ['#4F46E5' if i < 3 else '#94A3B8' for i in range(len(df_plot))]
        bars = ax.bar(range(len(df_plot)), df_plot["GC Content (%)"], color=colors, edgecolor='none', width=0.5, alpha=0.95)
        
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height}%', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 4), textcoords="offset points", ha='center', va='bottom', 
                        fontsize=8, fontweight='bold', color='#475569')
                        
        ax.set_title("Grafik Kelimpahan Distribusi Nilai GC Content", fontsize=11, fontweight='bold', color='#4F46E5', pad=15)
        ax.set_ylabel("Persentase GC (%)", fontsize=9, color='#475569')
        ax.set_xlabel("Peringkat Sekuens", fontsize=9, color='#475569')
        ax.set_xticks(range(len(df_plot)))
        ax.set_xticklabels(df_plot.index + 1, color='#475569')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color((0, 0, 0, 0.1))
        ax.spines['bottom'].set_color((0, 0, 0, 0.1))
        ax.tick_params(colors='#475569')
        ax.set_ylim(0, 110)
        ax.grid(True, linestyle="--", alpha=0.15, color='#000000', axis='y')
        st.pyplot(fig)
        
        st.markdown("<br>", unsafe_allow_html=True)
        csv_data = df_all.drop(columns=["raw_sequence"]).to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 EXPORT DATA PIPELINE (.CSV)",
            data=csv_data,
            file_name="rekap_pipeline_mini_project.csv",
            mime="text/csv",
            use_container_width=True
        )
        
        st.markdown("<br><h3 style='color: #4F46E5; font-size: 1.3rem; font-weight:700;'> Top 3 High GC Sequences</h3>", unsafe_allow_html=True)
        df_top3_display = df_all.head(3)
        st.dataframe(df_top3_display.drop(columns=["raw_sequence"]), use_container_width=True)

    with col_right:
        st.markdown("<h3 style='color: #4F46E5; font-size: 1.4rem; font-weight:700;'> 3D Spatial Conformation</h3>", unsafe_allow_html=True)
        st.markdown("Visualisasi interaktif struktur 3D protein homolog hasil konversi ID.")
        
        df_top3 = df_all.head(3)
        top3_options = {}
        for idx, row in df_top3.iterrows():
            label = f"Peringkat {idx+1} | {row['Header ID'][:25]}... ({row['GC Content (%)']}% GC)"
            top3_options[label] = row
            
        selected_3d_label = st.selectbox("Pilih Target Sekuens untuk Pemodelan 3D:", list(top3_options.keys()))
        selected_target_data = top3_options[selected_3d_label]
        
        clean_header = selected_target_data["Header ID"].lower()
        clean_file = selected_target_data["File Sumber"].lower()
        target_seq = selected_target_data["raw_sequence"]
        
        match_explicit = re.search(r'\b([0-9][a-z0-9]{3})\b', clean_header + " " + clean_file)
        if match_explicit:
            detected_pdb = match_explicit.group(1).upper()
        else:
            mapping_pdb = {"igf1": "3GD3", "igf-1": "3GD3", "1bna": "1BNA", "1imx": "1IMX"}
            matched_local = False
            for key, value in mapping_pdb.items():
                if key in clean_header or key in clean_file:
                    detected_pdb = value
                    matched_local = True
                    break
            if not matched_local:
                with st.spinner("Mencari model homolog melalui API RCSB..."):
                    api_result = fetch_pdb_by_nucleotide_blast(target_seq)
                    detected_pdb = api_result if api_result else "1BNA"

        gene_to_pdb_mapping = {"3479": "3GD3", "3480": "1IGR", "7157": "1TUP", "1956": "1IVO", "351": "1A82"}
        user_input = st.text_input("Masukkan Kode PDB ID atau Gene ID:", value=detected_pdb).strip()
        
        pdb_id = None
        if user_input.isdigit():
            if user_input in gene_to_pdb_mapping:
                pdb_id = gene_to_pdb_mapping[user_input]
                st.success(f" Gene ID '{user_input}' dikonversi ke PDB ID: {pdb_id}")
            else:
                st.error(f"❌ Gene ID '{user_input}' tidak ditemukan di database lokal.")
        else:
            if user_input: pdb_id = user_input

        style_type = st.selectbox("Gaya Pemodelan (Style Mode):", ["cartoon", "stick", "sphere", "line"])
        color_scheme = st.selectbox("Skema Warna Struktur Rantai:", ["spectrum", "chain", "ss"])

        if pdb_id and len(pdb_id) == 4:
            try:
                with st.spinner(f"Menarik koordinat model '{pdb_id.upper()}'..."):
                    view = py3Dmol.view(query=f'pdb:{pdb_id.lower()}', width=450, height=400)
                    if style_type == "cartoon": view.setStyle({'cartoon': {'color': color_scheme}})
                    elif style_type == "stick": view.setStyle({'stick': {}})
                    elif style_type == "sphere": view.setStyle({'sphere': {}})
                    else: view.setStyle({'line': {}})
                    view.setBackgroundColor('#1E293B') # Kontras gelap kokoh khusus canvas 3D mol agar molekul menyala terang
                    view.center()
                    view.zoomTo()
                    st.info(f"Visualisasi aktif referensi: PDB ID {pdb_id.upper()}")
                    showmol(view, height=400, width=450)
            except Exception:
                st.error("Gagal memuat model PDB.")
        else:
            st.warning("Sistem membutuhkan masukan kode registrasi PDB ID atau Gene ID yang valid.")
else:
    st.info("Sistem Pipeline Siap. Silakan unggah berkas sekuens (FASTA/FASTQ) di panel sidebar untuk memulai analisis otomatis.")