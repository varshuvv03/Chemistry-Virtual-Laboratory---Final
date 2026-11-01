# Cheminformatics Virtual Laboratory

**Developed by:** Varshini A (RA2611030050040) & Vismitha Parasa (RA2611030050047)  
**Institution:** SRM Institute of Science and Technology, Tiruchirappalli Campus  
**Department:** Department of Cybersecurity • Faculty of Engineering and Technology  

---

## Overview

This application is an exact, high-performance conversion of the Cheminformatics Virtual Laboratory created **without using Streamlit**. It uses **Flask** on the backend and modern **HTML5, CSS3, JavaScript, WebGL (3Dmol.js), and Chart.js** on the frontend, delivering instantaneous page transitions, smooth interactive 3D molecular rotations, and complete offline capability.

Original Streamlit file is safely backed up as `app_streamlit.py`.

---

## Features & Architecture

1. **Navigation (8 Dedicated Pages)**:
   - 🏠 **Home**: Welcome dashboard, responsive institutional header with SRM logo, 10 example molecule preview cards with 2D structures and SMILES codes, developer credentials table.
   - 📖 **Theory & Concepts**: 4 tabs covering Cheminformatics foundations, SMILES syntax and rules table with 6 2D rendered example cards, Molecular Descriptors & Lipinski's Rule of 5, and 3D Molecular Models with 4 interactive 3Dmol.js models of Aspirin side-by-side.
   - 🔬 **Molecular Visualization**: Real-time 2D molecular diagram and molecular formula/bond statistics alongside a full WebGL interactive 3D viewer (Ball & Stick, Stick, Space Filling, Wireframe, with optional van der Waals surface).
   - 📊 **Molecular Properties Report**: Complete calculation of all 21 molecular descriptors, 6 key drug-likeness metrics, CSV export, single-molecule PDF report download, and 3D viewer with static Matplotlib image generator.
   - 📈 **Structure–Property Analysis**: Multi-compound comparative matrix (13 analytical dimensions), CSV export, dynamic interactive property bar chart, Lipinski chemical space scatter plot (MW vs LogP), and Lipinski Rule of 5 compliance cards.
   - 🚀 **Advanced & Innovative Tools (Cybersecurity & Advanced Cheminformatics)**:
     - 🔐 Cryptographic Chemical Hashes (SHA-256, SHA3-256, MD5, InChIKey, InChI for molecular data integrity)
     - 🎯 Quantitative Estimate of Drug-Likeness (QED Score with continuous desirability parameter profiling)
     - 🧬 Bemis-Murcko Scaffold Deconstruction (Side-by-side original compound vs extracted cyclic core)
     - ⚡ Molecular Similarity Engine (Topological Daylight/Morgan fingerprint similarity: Tanimoto, Dice, Cosine metrics)
   - 📝 **Assessment Quiz**: All 24 questions with automated scoring, percentage grading, performance badges, and detailed review of missed questions with correct answers.
   - 📄 **Final Student Molecular Analysis Report**: Automated compilation and instant download of the complete 10-page institutional laboratory report PDF with official SRMIST headers, dual-pass NumberedCanvas pagination, tables, and faculty evaluation signature blocks.

---

## How to Run

1. **Install requirements** (if not already installed):
   ```bash
   pip install -r requirements.txt
   ```

2. **Start the application**:
   ```bash
   python app.py
   ```

3. **Open in browser**:
   Navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000)

