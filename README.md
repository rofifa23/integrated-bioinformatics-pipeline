# Integrated Bioinformatics Analytics & Pipeline Suite 

link website: https://pipeline-bioinformatika-terintegrasi.streamlit.app/

Repositori ini berisi kode program dan dokumentasi teknis untuk *Mini Project* mata kuliah Struktur Data Bioinformatika, Departemen Biokimia, IPB University. Aplikasi ini mengotomasi *pipeline* analisis sekuens asam nukleat mentah hingga visualisasi konformasi spasial makromolekul secara interaktif.

##  Fitur Utama Aplikasi
1. **Otomasi File Parser (FASTA/FASTQ)**: Membaca berkas sekuens mentah secara asinkronus dan memisahkan baris *header identifier* dengan data sekuens biologis.
2. **Kalkulasi GC Content Kuantitatif**: Menghitung frekuensi kemunculan nukleotida (basa A, T, C, G) menggunakan struktur data *dictionary* untuk merekam skor kelimpahan, serta menghitung rasio persentase nilai GC.
3. **Pipeline Data Sorting & Integration**: Menggabungkan seluruh hasil parsing ke dalam format *List of Dictionaries*, mengurutkannya berdasarkan nilai GC Content tertinggi secara menurun (*Max Heap Priority Principle*), dan menyediakan fitur ekspor data terintegrasi ke berkas `.CSV`.
4. **Visualisasi Interaktif 3D Spasial**: Integrasi langsung dengan API RCSB Protein Data Bank (PDB) untuk memetakan ID gene homolog atau PDB ID menjadi model struktural protein 3D interaktif yang dapat dirotasi secara *real-time*.

##  Struktur Berkas Projek
```text
├── IGF1_datasets/                  # Folder dataset sekuens uji coba
├── Aplication project1.py          # Script utama aplikasi Streamlit
├── logo ipb horizontal.png         # Aset logo IPB University untuk UI
├── nutrition_jurnal pendukung...   # Dokumen literatur pendukung analisis
└── README.md                       # Dokumentasi repositori (Berkas ini)

Prasyarat & Cara Menjalankan Aplikasi
Bagian ini menjelaskan langkah penyiapan lingkungan lokal untuk mengeksekusi program.

Unduh/Clone Repositori ini

Bash
git clone [https://github.com/rofifa23/integrated-bioinformatics-pipeline.git](https://github.com/rofifa23/integrated-bioinformatics-pipeline.git)
cd integrated-bioinformatics-pipeline
Instalasi Pustaka (Dependencies)
Jalankan perintah berikut pada terminal untuk menginstal pustaka yang diperlukan:

Bash
pip install streamlit matplotlib pandas stmol py3Dmol
Eksekusi Aplikasi
Jalankan server lokal Streamlit dengan perintah:

Bash
streamlit run "Aplication project1.py"

Teknologi & Library yang Digunakan
Bahasa Pemrograman: Python 3.13
Framework Dashboard: Streamlit (Kustomisasi Premium CSS UI)
Visualisasi Statistik: Matplotlib (Light theme minimalis)
Manajemen Data: Pandas DataFrame
Rendering 3D Makromolekul: py3Dmol & stmol
