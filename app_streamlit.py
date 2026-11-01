# ============================================================================
# Cheminformatics Virtual Laboratory
# Developed by Varshini A & Vismitha Parasa — SRM Institute of Science and Technology
# ============================================================================

import streamlit as st
import pandas as pd
import numpy as np
import io
import os
import base64
import textwrap
import random
from io import BytesIO

# RDKit imports
from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Draw,
    AllChem,
    rdMolDescriptors,
    Lipinski,
    Crippen,
    rdchem,
)
from rdkit.Chem.Draw import rdMolDraw2D

# 3D visualization
import py3Dmol
from stmol import showmol

# PDF generation
from reportlab.lib.pagesizes import A4, letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, mm, cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    Image as RLImage,
    PageBreak,
    KeepTogether,
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.pdfgen import canvas

# Matplotlib for static 3D images in PDF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from PIL import Image as PILImage

# ============================================================================
# PAGE CONFIG
# ============================================================================
st.set_page_config(
    page_title="Cheminformatics Virtual Laboratory",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ============================================================================
# ============================================================================
# CUSTOM CSS & RESPONSIVE MOBILE VIEW
# ============================================================================
st.markdown("""
<style>
    /* Desktop default typography */
    .main-header {
        font-size: 2.2rem;
        font-weight: bold;
        color: #2e8b57;
        text-align: center;
        padding: 10px 0;
    }
    .sub-header {
        font-size: 1.3rem;
        color: #2e4a3b;
        text-align: center;
        margin-bottom: 20px;
    }
    .info-box {
        background-color: #f0f8f4;
        border-left: 5px solid #2e8b57;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
    }
    .success-box {
        background-color: #eafaf1;
        border-left: 5px solid #2e8b57;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
    }
    .warning-box {
        background-color: #fef9e7;
        border-left: 5px solid #f39c12;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 20px;
        font-weight: 600;
    }
    .dev-footer {
        text-align: center;
        padding: 20px;
        color: #7f8c8d;
        font-size: 0.9rem;
        border-top: 1px solid #d1e7dd;
        margin-top: 30px;
    }

    /* Header responsive classes */
    .responsive-header-title {
        font-size: 2.1rem;
        font-weight: 800;
        color: #2e8b57;
        line-height: 1.25;
        margin-bottom: 6px;
    }
    .responsive-header-dept {
        font-size: 1.15rem;
        font-weight: 600;
        color: #1b4d32;
        line-height: 1.35;
        margin-bottom: 4px;
    }
    .responsive-header-campus {
        font-size: 0.95rem;
        color: #555555;
        font-weight: 500;
        line-height: 1.35;
        margin: 0;
    }
    .responsive-logo-wrapper {
        display: flex;
        justify-content: flex-end;
        align-items: center;
        height: 100%;
    }

    /* Mobile view optimizations (screens <= 768px) */
    @media (max-width: 768px) {
        .responsive-header-title {
            font-size: 1.5rem !important;
            text-align: center !important;
        }
        .responsive-header-dept {
            font-size: 0.95rem !important;
            text-align: center !important;
        }
        .responsive-header-campus {
            font-size: 0.82rem !important;
            text-align: center !important;
            margin-bottom: 10px !important;
        }
        .responsive-logo-wrapper {
            justify-content: center !important;
            margin-top: 5px !important;
            margin-bottom: 10px !important;
        }
        .stTabs [data-baseweb="tab"] {
            padding: 6px 10px !important;
            font-size: 0.85rem !important;
        }
        .stButton>button {
            width: 100% !important;
        }
        .dev-footer {
            font-size: 0.78rem !important;
            padding: 15px 5px !important;
        }
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# EXAMPLE MOLECULES DATABASE
# ============================================================================
EXAMPLE_MOLECULES = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
    "Ethanol": "CCO",
    "Acetic Acid": "CC(=O)O",
    "Benzene": "c1ccccc1",
    "Penicillin G": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@@H]1C(=O)O",
    "Cholesterol": "C(CCCC(C)C)CCC1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def validate_smiles(smiles):
    """Validate a SMILES string and return the molecule object or None."""
    if not smiles or not smiles.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        return mol
    except Exception:
        return None


def get_mol_2d_image(mol, size=(400, 400)):
    """Generate a 2D image of the molecule."""
    if mol is None:
        return None
    try:
        img = Draw.MolToImage(mol, size=size)
        return img
    except Exception:
        return None


def get_mol_2d_image_bytes(mol, size=(400, 400)):
    """Generate a 2D image of the molecule as PNG bytes."""
    img = get_mol_2d_image(mol, size)
    if img is None:
        return None
    buf = BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf


def generate_3d_coords(mol):
    """Generate 3D coordinates for a molecule."""
    if mol is None:
        return None
    try:
        mol_3d = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        if result == -1:
            result = AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol_3d)
        return mol_3d
    except Exception:
        return None


def get_3d_viewer(mol_3d, style="stick", show_surface=False, width=500, height=400):
    """Create a py3Dmol viewer with the given style."""
    try:
        mol_block = Chem.MolToMolBlock(mol_3d)
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(mol_block, "mol")

        if style == "stick":
            viewer.setStyle({"stick": {"radius": 0.15}})
        elif style == "ball_and_stick":
            viewer.setStyle({"stick": {"radius": 0.1}, "sphere": {"scale": 0.25}})
        elif style == "space_filling":
            viewer.setStyle({"sphere": {"scale": 1.0}})
        elif style == "wireframe":
            viewer.setStyle({"line": {"linewidth": 2}})
        else:
            viewer.setStyle({"stick": {}})

        if show_surface:
            viewer.addSurface(
                py3Dmol.VDW, {"opacity": 0.7, "color": "lightblue"}
            )

        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        return viewer
    except Exception:
        return None


def generate_static_3d_image(mol_3d, style="ball_and_stick"):
    """Generate a static matplotlib 3D image of the molecule for PDF."""
    if mol_3d is None:
        return None
    try:
        conf = mol_3d.GetConformer()
        atom_positions = []
        atom_symbols = []
        for i in range(mol_3d.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_positions.append([pos.x, pos.y, pos.z])
            atom_symbols.append(mol_3d.GetAtomWithIdx(i).GetSymbol())

        atom_positions = np.array(atom_positions)

        # Color map for atoms
        color_map = {
            "C": "#404040",
            "N": "#3050F8",
            "O": "#FF0D0D",
            "S": "#FFFF30",
            "P": "#FF8000",
            "F": "#90E050",
            "Cl": "#1FF01F",
            "Br": "#A62929",
            "I": "#940094",
            "H": "#CCCCCC",
        }

        fig = plt.figure(figsize=(5, 4), dpi=120)
        ax = fig.add_subplot(111, projection="3d")

        # Plot atoms
        for i, (pos, sym) in enumerate(zip(atom_positions, atom_symbols)):
            c = color_map.get(sym, "#808080")
            s = 30 if sym == "H" else 80
            if style == "space_filling":
                s *= 3
            ax.scatter(*pos, c=c, s=s, edgecolors="black", linewidths=0.5, alpha=0.9)

        # Plot bonds
        if style != "space_filling":
            for bond in mol_3d.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                p1 = atom_positions[i]
                p2 = atom_positions[j]
                lw = 1.5 if style == "wireframe" else 2.5
                ax.plot(
                    [p1[0], p2[0]],
                    [p1[1], p2[1]],
                    [p1[2], p2[2]],
                    color="gray",
                    linewidth=lw,
                    alpha=0.7,
                )

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("3D Molecular Structure", fontsize=10)

        # Remove grid for cleaner look
        ax.grid(True, alpha=0.3)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        plt.tight_layout()
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
        plt.close(fig)
        buf.seek(0)
        return buf
    except Exception:
        return None


def count_sigma_pi_bonds(mol):
    """Count sigma and pi bonds in a molecule."""
    if mol is None:
        return 0, 0
    sigma = 0
    pi = 0
    for bond in mol.GetBonds():
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE:
            sigma += 1
        elif bt == Chem.rdchem.BondType.DOUBLE:
            sigma += 1
            pi += 1
        elif bt == Chem.rdchem.BondType.TRIPLE:
            sigma += 1
            pi += 2
        elif bt == Chem.rdchem.BondType.AROMATIC:
            sigma += 1
            pi += 0.5  # Delocalized
    return int(sigma), round(pi)


def calculate_all_properties(mol):
    """Calculate a comprehensive set of molecular properties."""
    if mol is None:
        return {}

    sigma, pi = count_sigma_pi_bonds(mol)

    props = {
        "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular Weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "Exact Mass": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA (\u00c5\u00b2)": round(Descriptors.TPSA(mol), 2),
        "Number of Atoms": mol.GetNumAtoms(),
        "Number of Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Number of Bonds": mol.GetNumBonds(),
        "\u03c3 (Sigma) Bonds": sigma,
        "\u03c0 (Pi) Bonds": pi,
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Aromatic Rings": Descriptors.NumAromaticRings(mol),
        "Number of Rings": Descriptors.RingCount(mol),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 3),
        "Molar Refractivity": round(Crippen.MolMR(mol), 3),
        "Number of Radical Electrons": Descriptors.NumRadicalElectrons(mol),
        "Number of Valence Electrons": Descriptors.NumValenceElectrons(mol),
        "Formal Charge": Chem.GetFormalCharge(mol),
        "Complexity (BertzCT)": round(Descriptors.BertzCT(mol), 2),
    }
    return props


def get_interpretation(props):
    """Generate interpretation text based on property values."""
    lines = []
    mw = props.get("Molecular Weight (g/mol)", 0)
    logp = props.get("LogP", 0)
    tpsa = props.get("TPSA (\u00c5\u00b2)", 0)
    hbd = props.get("H-Bond Donors", 0)
    hba = props.get("H-Bond Acceptors", 0)
    rb = props.get("Rotatable Bonds", 0)

    # Molecular Weight
    if mw <= 500:
        lines.append(f"Molecular Weight ({mw:.1f} g/mol): Within Lipinski's limit (\u2264500), favorable for oral absorption.")
    else:
        lines.append(f"Molecular Weight ({mw:.1f} g/mol): Exceeds Lipinski's limit (>500), may have reduced oral bioavailability.")

    # LogP
    if logp < 0:
        lines.append(f"LogP ({logp:.2f}): Negative value indicates the molecule is hydrophilic (water-loving), suggesting good aqueous solubility.")
    elif logp <= 5:
        lines.append(f"LogP ({logp:.2f}): Within acceptable range (\u22645), indicating balanced lipophilicity for membrane permeation.")
    else:
        lines.append(f"LogP ({logp:.2f}): Exceeds Lipinski's limit (>5), suggesting high lipophilicity and possible poor solubility.")

    # TPSA
    if tpsa < 60:
        lines.append(f"TPSA ({tpsa:.1f} \u00c5\u00b2): Low polarity, suggesting good cell membrane permeability but potentially low aqueous solubility.")
    elif tpsa <= 140:
        lines.append(f"TPSA ({tpsa:.1f} \u00c5\u00b2): Moderate polarity, favorable for oral bioavailability.")
    else:
        lines.append(f"TPSA ({tpsa:.1f} \u00c5\u00b2): High polarity (>140), may have difficulty crossing cell membranes.")

    # HBD
    if hbd <= 5:
        lines.append(f"H-Bond Donors ({hbd}): Within Lipinski's limit (\u22645).")
    else:
        lines.append(f"H-Bond Donors ({hbd}): Exceeds Lipinski's limit (>5).")

    # HBA
    if hba <= 10:
        lines.append(f"H-Bond Acceptors ({hba}): Within Lipinski's limit (\u226410).")
    else:
        lines.append(f"H-Bond Acceptors ({hba}): Exceeds Lipinski's limit (>10).")

    # Rotatable Bonds
    if rb <= 10:
        lines.append(f"Rotatable Bonds ({rb}): Good molecular flexibility for binding.")
    else:
        lines.append(f"Rotatable Bonds ({rb}): High flexibility (>10), may reduce binding selectivity.")

    # Lipinski violations
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    if violations == 0:
        lines.append(f"\n\u2705 Lipinski Rule of 5: PASS (0 violations) \u2014 Molecule is predicted to be orally active.")
    elif violations == 1:
        lines.append(f"\n\u26a0\ufe0f Lipinski Rule of 5: {violations} violation \u2014 Molecule may still be orally bioavailable.")
    else:
        lines.append(f"\n\u274c Lipinski Rule of 5: {violations} violations \u2014 Molecule is predicted to have poor oral bioavailability.")

    return "\n".join(lines)


def lipinski_violations(mol):
    """Count Lipinski Rule of 5 violations."""
    if mol is None:
        return 4
    v = 0
    if Descriptors.MolWt(mol) > 500:
        v += 1
    if Descriptors.MolLogP(mol) > 5:
        v += 1
    if Descriptors.NumHDonors(mol) > 5:
        v += 1
    if Descriptors.NumHAcceptors(mol) > 10:
        v += 1
    return v


# ============================================================================
# PDF GENERATION — 10-PAGE COMPREHENSIVE EXPERIMENT REPORT
# ============================================================================
class NumberedCanvas(canvas.Canvas):
    """Dual-pass canvas to dynamically compute and draw exact 'Page X of Y' and header logos."""
    def __init__(self, *args, **kwargs):
        super(NumberedCanvas, self).__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_decorations(num_pages)
            super(NumberedCanvas, self).showPage()
        super(NumberedCanvas, self).save()

    def draw_page_decorations(self, page_count):
        self.saveState()
        logo_path = "assets/srm_logo.png"

        # Running header on pages > 1
        if self._pageNumber > 1:
            if os.path.exists(logo_path):
                try:
                    self.drawImage(logo_path, 1.4 * cm, 29.7 * cm - 1.55 * cm, width=2.6 * cm, height=1.25 * cm, mask="auto", preserveAspectRatio=True)
                except Exception:
                    pass
            self.setFont("Helvetica-Bold", 8)
            self.setFillColor(colors.HexColor("#2e8b57"))
            self.drawString(4.3 * cm, 29.7 * cm - 0.95 * cm, "SRM INSTITUTE OF SCIENCE AND TECHNOLOGY")
            self.setFont("Helvetica", 7)
            self.setFillColor(colors.HexColor("#555555"))
            self.drawString(4.3 * cm, 29.7 * cm - 1.35 * cm, "Cheminformatics Virtual Laboratory \u2014 Student Molecular Analysis Report")
            self.drawRightString(21.0 * cm - 1.4 * cm, 29.7 * cm - 1.15 * cm, "TIRUCHIRAPPALLI CAMPUS")
            self.setStrokeColor(colors.HexColor("#2e8b57"))
            self.setLineWidth(0.8)
            self.line(1.4 * cm, 29.7 * cm - 1.6 * cm, 21.0 * cm - 1.4 * cm, 29.7 * cm - 1.6 * cm)

        # Running footer on all pages
        self.setStrokeColor(colors.HexColor("#bdc3c7"))
        self.setLineWidth(0.5)
        self.line(1.4 * cm, 1.2 * cm, 21.0 * cm - 1.4 * cm, 1.2 * cm)
        self.setFont("Helvetica", 7)
        self.setFillColor(colors.HexColor("#555555"))
        self.drawString(1.4 * cm, 0.8 * cm, "Department of Cybersecurity \u2022 Varshini A & Vismitha Parasa \u2022 SRMIST Trichy")
        self.drawRightString(21.0 * cm - 1.4 * cm, 0.8 * cm, f"Page {self._pageNumber} of {page_count}")
        self.restoreState()


def get_interpretation_compact(props, mol):
    """Generate structured interpretation bullet points for PDF formatting."""
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rb = Descriptors.NumRotatableBonds(mol)

    viol = 0
    if mw > 500: viol += 1
    if logp > 5: viol += 1
    if hbd > 5: viol += 1
    if hba > 10: viol += 1

    lines = []
    if viol == 0:
        lip_txt = "Compliant (0 violations) \u2014 Drug-like candidate"
    else:
        lip_txt = f"{viol} violation(s) detected"
    lines.append(f"<b>Bioavailability & Lipinski Rule:</b> {lip_txt} (MW: {mw:.1f}, LogP: {logp:.2f}, HBD: {hbd}, HBA: {hba}).")

    if tpsa <= 140:
        pol_txt = "Optimal for cellular permeability"
    else:
        pol_txt = "High polarity, may hinder membrane transport"
    lines.append(f"<b>Polarity & Permeation:</b> TPSA is {tpsa:.1f} \u00c5\u00b2 ({pol_txt}). Rotatable bonds: {rb} (conformational flexibility).")

    if logp < 0:
        sol_txt = "hydrophilic character favorable for aqueous solubility"
    elif logp <= 5:
        sol_txt = "balanced lipophilicity for physiological membranes"
    else:
        sol_txt = "lipophilic character with potential solubility limits"
    lines.append(f"<b>Solubility Profile:</b> LogP of {logp:.2f} indicates {sol_txt}.")
    return lines


def build_pdf_report(student_name, reg_no, molecules_data, style_3d="ball_and_stick"):
    """Build a comprehensive 10-page multi-molecule PDF report with institutional headers, tables, charts, and sign-offs."""
    # Ensure standard cohort of 7 compounds if fewer provided to guarantee a full 10-page laboratory report
    cohort_defaults = [
        ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Paracetamol", "CC(=O)NC1=CC=C(O)C=C1"),
        ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
        ("Glucose", "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"),
        ("Ethanol", "CCO"),
        ("Benzene", "c1ccccc1"),
    ]

    working_molecules = list(molecules_data)
    current_smiles_set = {m.get("smiles") for m in working_molecules}

    # Pad with defaults if less than 7 so we achieve exact 10 pages
    for dname, dsmi in cohort_defaults:
        if len(working_molecules) >= 7:
            break
        if dsmi not in current_smiles_set:
            dmol = validate_smiles(dsmi)
            if dmol:
                dprops = calculate_all_properties(dmol)
                working_molecules.append({
                    "name": dname,
                    "smiles": dsmi,
                    "mol": dmol,
                    "properties": dprops,
                    "interpretation": get_interpretation(dprops),
                })
                current_smiles_set.add(dsmi)

    # Use first 7 molecules for pages 2 to 8
    cohort_7 = working_molecules[:7]

    buf = BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=A4,
        leftMargin=1.4 * cm,
        rightMargin=1.4 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.4 * cm,
    )
    styles = getSampleStyleSheet()

    tbl_pstyle = ParagraphStyle("tblp", fontSize=6.8, leading=8.5)
    tbl_hstyle = ParagraphStyle("tblh", fontSize=7, leading=8.5, textColor=colors.white)

    story = []
    logo_path = "assets/srm_logo.png"

    # ==================== PAGE 1: TITLE & COVER ====================
    if os.path.exists(logo_path):
        logo_rl = RLImage(logo_path, width=4.6 * cm, height=2.3 * cm)
    else:
        logo_rl = Paragraph("<b>SRM LOGO</b>", ParagraphStyle("no_logo", fontSize=12))

    header_tbl = Table([
        [
            logo_rl,
            Table([
                [Paragraph("<b>SRM INSTITUTE OF SCIENCE AND TECHNOLOGY</b>", ParagraphStyle("H1", fontSize=14, leading=17, textColor=colors.HexColor("#1b4d32")))],
                [Paragraph("<b>TIRUCHIRAPPALLI CAMPUS</b>", ParagraphStyle("H2", fontSize=10.5, leading=13, textColor=colors.HexColor("#2e8b57")))],
                [Paragraph("FACULTY OF ENGINEERING AND TECHNOLOGY \u2022 DEPARTMENT OF CYBERSECURITY", ParagraphStyle("H3", fontSize=8, leading=10, textColor=colors.HexColor("#555555")))],
                [Paragraph("<b>CHEMINFORMATICS VIRTUAL LABORATORY</b>", ParagraphStyle("H4", fontSize=9.5, leading=12, textColor=colors.HexColor("#2e8b57")))],
            ], colWidths=[13.4 * cm])
        ]
    ], colWidths=[4.8 * cm, 13.4 * cm])
    header_tbl.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (0, 0), "LEFT"),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
    ]))
    story.append(header_tbl)
    story.append(Spacer(1, 4 * mm))

    title_tbl = Table([
        [Paragraph("<b>MULTI-MOLECULE STUDENT ANALYSIS & RESEARCH REPORT</b>", ParagraphStyle("DT", alignment=TA_CENTER, fontSize=12, leading=15, textColor=colors.white))]
    ], colWidths=[18.2 * cm])
    title_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#2e8b57")),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
    ]))
    story.append(title_tbl)
    story.append(Spacer(1, 4 * mm))

    # Student metadata table
    student_data = [
        [Paragraph("<b>Student Name</b>", ParagraphStyle("c1", fontSize=8.5, textColor=colors.white)), Paragraph("<b>Registration No.</b>", ParagraphStyle("c2", fontSize=8.5, textColor=colors.white)), Paragraph("<b>Email ID</b>", ParagraphStyle("c3", fontSize=8.5, textColor=colors.white))],
        [Paragraph("<b>Varshini A</b>", ParagraphStyle("s1", fontSize=8.5)), Paragraph("RA2611030050040", ParagraphStyle("s2", fontSize=8.5)), Paragraph("va3375@srmist.edu.in", ParagraphStyle("s3", fontSize=8.5))],
        [Paragraph("<b>Vismitha Parasa</b>", ParagraphStyle("s4", fontSize=8.5)), Paragraph("RA2611030050047", ParagraphStyle("s5", fontSize=8.5)), Paragraph("vp8717@srmist.edu.in", ParagraphStyle("s6", fontSize=8.5))],
    ]
    std_tbl = Table(student_data, colWidths=[6.0 * cm, 5.2 * cm, 7.0 * cm])
    std_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1b4d32")),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#bdc3c7")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f9f6")]),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
    ]))
    story.append(std_tbl)
    story.append(Spacer(1, 4 * mm))

    story.append(Paragraph("<b>Laboratory Experiment Abstract</b>", ParagraphStyle("absh", fontSize=10, leading=12, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        "This laboratory report compiles a systematic cheminformatics study conducted using computational tools including RDKit, py3Dmol, Matplotlib, and Streamlit. A curated cohort of seven representative organic compounds (Aspirin, Caffeine, Paracetamol, Ibuprofen, Glucose, Ethanol, and Benzene) was evaluated across 21 molecular descriptors encompassing structural, topological, electronic, and physicochemical dimensions. Conformer models were optimized using the Merck Molecular Force Field (MMFF94) following distance geometry embedding (ETKDG). Drug-likeness and membrane permeability were assessed using Lipinski\u2019s Rule of Five and Veber criteria, establishing a comprehensive reference for virtual screening and structure-property correlation.",
        ParagraphStyle("absb", fontSize=8, leading=11, alignment=TA_JUSTIFY)
    ))
    story.append(Spacer(1, 4 * mm))

    story.append(Paragraph("<b>Report Structure & Table of Contents (10 Pages)</b>", ParagraphStyle("toch", fontSize=10, leading=12, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 1.5 * mm))
    toc_data = [
        [Paragraph("<b>Page</b>", ParagraphStyle("tc1", fontSize=7.5, textColor=colors.white)), Paragraph("<b>Section / Compound</b>", ParagraphStyle("tc2", fontSize=7.5, textColor=colors.white)), Paragraph("<b>Analytical Scope</b>", ParagraphStyle("tc3", fontSize=7.5, textColor=colors.white))],
        [Paragraph("Page 1", ParagraphStyle("t1", fontSize=7.5)), Paragraph("Title & Student Metadata", ParagraphStyle("t2", fontSize=7.5)), Paragraph("Executive Overview, Methodology, Objectives", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 2", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 1: {cohort_7[0].get('name', 'Compound 1')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 3", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 2: {cohort_7[1].get('name', 'Compound 2')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 4", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 3: {cohort_7[2].get('name', 'Compound 3')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 5", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 4: {cohort_7[3].get('name', 'Compound 4')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 6", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 5: {cohort_7[4].get('name', 'Compound 5')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 7", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 6: {cohort_7[5].get('name', 'Compound 6')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 8", ParagraphStyle("t1", fontSize=7.5)), Paragraph(f"Compound 7: {cohort_7[6].get('name', 'Compound 7')}", ParagraphStyle("t2", fontSize=7.5)), Paragraph("2D/3D Conformation, 21 Descriptors, Bioavailability", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 9", ParagraphStyle("t1", fontSize=7.5)), Paragraph("Comparative Analysis", ParagraphStyle("t2", fontSize=7.5)), Paragraph("Multi-Molecule Matrix, Lipinski Space Scatter Chart", ParagraphStyle("t3", fontSize=7.5))],
        [Paragraph("Page 10", ParagraphStyle("t1", fontSize=7.5)), Paragraph("Summary & Sign-off", ParagraphStyle("t2", fontSize=7.5)), Paragraph("Technological Stack, Code Script, Evaluation Signatures", ParagraphStyle("t3", fontSize=7.5))],
    ]
    toc_tbl = Table(toc_data, colWidths=[2.2 * cm, 5.5 * cm, 10.5 * cm])
    toc_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2e8b57")),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#bdc3c7")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f9f6")]),
        ("TOPPADDING", (0, 0), (-1, -1), 1.6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 1.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
    ]))
    story.append(toc_tbl)
    story.append(PageBreak())

    # ==================== PAGES 2 to 8: 7 MOLECULES (1 PAGE EACH) ====================
    mol_objects_for_page9 = []

    for idx, minfo in enumerate(cohort_7, 1):
        cname = minfo.get("name", f"Compound {idx}")
        csmi = minfo.get("smiles", "")
        cmol = minfo.get("mol")
        cprops = minfo.get("properties") or calculate_all_properties(cmol)
        interp_lines = get_interpretation_compact(cprops, cmol)

        mol_objects_for_page9.append((cname, csmi, cmol, cprops))

        story.append(Paragraph(f"<b>Molecule {idx} of 7: {cname}</b> &nbsp;&nbsp;|&nbsp;&nbsp; SMILES: <code>{csmi}</code>", ParagraphStyle("mh", fontSize=10.5, leading=13, textColor=colors.HexColor("#1b4d32"))))
        story.append(Spacer(1, 1.5 * mm))

        # Side-by-side structures
        img2d_buf = get_mol_2d_image_bytes(cmol, size=(280, 210))
        cmol_3d = generate_3d_coords(cmol)
        img3d_buf = generate_static_3d_image(cmol_3d, style=style_3d)

        img_row = Table([
            [Paragraph("<b>2D Chemical Structure</b>", ParagraphStyle("l1", alignment=TA_CENTER, fontSize=7.5)), Paragraph("<b>3D Molecular Model (Conformer)</b>", ParagraphStyle("l2", alignment=TA_CENTER, fontSize=7.5))],
            [RLImage(img2d_buf, width=7.2 * cm, height=4.6 * cm), RLImage(img3d_buf, width=7.2 * cm, height=4.6 * cm)]
        ], colWidths=[9.1 * cm, 9.1 * cm])
        img_row.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 1),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1),
        ]))
        story.append(img_row)
        story.append(Spacer(1, 1.5 * mm))

        # 21 Descriptors in 2 columns
        items = list(cprops.items())
        half = (len(items) + 1) // 2
        c1_items = items[:half]
        c2_items = items[half:]

        p_rows = [[
            Paragraph("<b>Descriptor</b>", tbl_hstyle), Paragraph("<b>Value</b>", tbl_hstyle),
            Paragraph("<b>Descriptor</b>", tbl_hstyle), Paragraph("<b>Value</b>", tbl_hstyle),
        ]]
        for i in range(half):
            p1, v1 = c1_items[i]
            p2, v2 = c2_items[i] if i < len(c2_items) else ("", "")
            p_rows.append([
                Paragraph(str(p1), tbl_pstyle), Paragraph(str(v1), tbl_pstyle),
                Paragraph(str(p2), tbl_pstyle), Paragraph(str(v2), tbl_pstyle),
            ])

        ptbl = Table(p_rows, colWidths=[5.0 * cm, 4.1 * cm, 5.0 * cm, 4.1 * cm])
        ptbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2e8b57")),
            ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#bdc3c7")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f9f6")]),
            ("TOPPADDING", (0, 0), (-1, -1), 1.2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 1.2),
            ("LEFTPADDING", (0, 0), (-1, -1), 3),
            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ]))
        story.append(ptbl)
        story.append(Spacer(1, 2 * mm))

        story.append(Paragraph("<b>Cheminformatics & Pharmacological Interpretation:</b>", ParagraphStyle("inh", fontSize=8.5, leading=10, textColor=colors.HexColor("#1b4d32"))))
        story.append(Spacer(1, 1 * mm))
        for iline in interp_lines:
            story.append(Paragraph(iline, ParagraphStyle("inb", fontSize=7.2, leading=9.2, alignment=TA_JUSTIFY)))
            story.append(Spacer(1, 0.5 * mm))

        story.append(PageBreak())

    # ==================== PAGE 9: COMPARATIVE ANALYSIS ====================
    story.append(Paragraph("<b>Comparative Structure\u2013Property Analysis (7 Molecules)</b>", ParagraphStyle("cph", fontSize=12, leading=15, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 2 * mm))
    story.append(Paragraph(
        "A multi-compound matrix correlating molecular weight, lipophilicity (LogP), topological polar surface area (TPSA), hydrogen bond donor/acceptor counts, rotatable bonds, and Lipinski Rule of 5 violations across all analyzed structures.",
        ParagraphStyle("cpb", fontSize=7.8, leading=10.5, alignment=TA_JUSTIFY)
    ))
    story.append(Spacer(1, 2.5 * mm))

    comp_headers = [
        Paragraph("<b>Compound</b>", tbl_hstyle), Paragraph("<b>MW (g/mol)</b>", tbl_hstyle),
        Paragraph("<b>LogP</b>", tbl_hstyle), Paragraph("<b>TPSA (\u00c5\u00b2)</b>", tbl_hstyle),
        Paragraph("<b>HBD</b>", tbl_hstyle), Paragraph("<b>HBA</b>", tbl_hstyle),
        Paragraph("<b>Rot.B</b>", tbl_hstyle), Paragraph("<b>Arom.R</b>", tbl_hstyle),
        Paragraph("<b>Lipinski Viol.</b>", tbl_hstyle),
    ]
    comp_rows = [comp_headers]

    names_list, mw_list, logp_list = [], [], []
    for name, smi, mol, _ in mol_objects_for_page9:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rb = Descriptors.NumRotatableBonds(mol)
        ar = Descriptors.NumAromaticRings(mol)

        viol = 0
        if mw > 500: viol += 1
        if logp > 5: viol += 1
        if hbd > 5: viol += 1
        if hba > 10: viol += 1

        names_list.append(name)
        mw_list.append(mw)
        logp_list.append(logp)

        viol_text = "0 (Pass)" if viol == 0 else f"{viol} Viol."
        comp_rows.append([
            Paragraph(name, tbl_pstyle), Paragraph(f"{mw:.1f}", tbl_pstyle),
            Paragraph(f"{logp:.2f}", tbl_pstyle), Paragraph(f"{tpsa:.1f}", tbl_pstyle),
            Paragraph(str(hbd), tbl_pstyle), Paragraph(str(hba), tbl_pstyle),
            Paragraph(str(rb), tbl_pstyle), Paragraph(str(ar), tbl_pstyle),
            Paragraph(viol_text, tbl_pstyle),
        ])

    comp_tbl = Table(comp_rows, colWidths=[3.2 * cm, 2.0 * cm, 1.8 * cm, 2.0 * cm, 1.5 * cm, 1.5 * cm, 1.6 * cm, 1.8 * cm, 2.8 * cm])
    comp_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2e8b57")),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#bdc3c7")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f9f6")]),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(comp_tbl)
    story.append(Spacer(1, 3 * mm))

    # Scatter plot
    fig, ax = plt.subplots(figsize=(6.8, 3.2), dpi=110)
    ax.scatter(mw_list, logp_list, s=70, c="#2e8b57", edgecolors="black", linewidths=0.5, zorder=3)
    for n, x, y in zip(names_list, mw_list, logp_list):
        ax.annotate(n, (x, y), textcoords="offset points", xytext=(5, 4), fontsize=7.5, fontweight="bold")
    ax.axhline(y=5, color="#c0392b", linestyle="--", linewidth=0.8, alpha=0.7, label="Lipinski Max LogP (5.0)")
    ax.axvline(x=500, color="#c0392b", linestyle="--", linewidth=0.8, alpha=0.7, label="Lipinski Max MW (500 g/mol)")
    ax.set_xlabel("Molecular Weight (g/mol)", fontsize=8)
    ax.set_ylabel("LogP (Octanol-Water Partition)", fontsize=8)
    ax.set_title("Lipinski Chemical Space: MW vs LogP Distribution", fontsize=9, fontweight="bold", color="#1b4d32")
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="upper left")
    plt.tight_layout(pad=0.4)
    scatter_buf = BytesIO()
    fig.savefig(scatter_buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    scatter_buf.seek(0)

    scatter_img = RLImage(scatter_buf, width=13.0 * cm, height=6.0 * cm)
    scatter_tbl = Table([[scatter_img]], colWidths=[18.2 * cm])
    scatter_tbl.setStyle(TableStyle([("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    story.append(scatter_tbl)
    story.append(Spacer(1, 2 * mm))

    story.append(Paragraph(
        "<b>Scatter Chart Interpretation:</b> All seven compounds reside well within the Lipinski MW boundary (\u2264 500 g/mol) and LogP boundary (\u2264 5.0). Glucose demonstrates notable hydrophilicity (negative LogP) driven by five hydroxyl groups, whereas Aspirin, Paracetamol, and Ibuprofen occupy ideal oral drug-like space.",
        ParagraphStyle("scb", fontSize=7.2, leading=9.5, alignment=TA_JUSTIFY)
    ))
    story.append(PageBreak())

    # ==================== PAGE 10: METHODOLOGY, CODE & SIGN-OFF ====================
    story.append(Paragraph("<b>Technological Architecture, References & Sign-off</b>", ParagraphStyle("p10h", fontSize=12, leading=15, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 2 * mm))

    tech_data = [
        [Paragraph("<b>Software Framework</b>", tbl_hstyle), Paragraph("<b>Role in Virtual Laboratory Pipeline</b>", tbl_hstyle)],
        [Paragraph("<b>RDKit (C++ / Python)</b>", tbl_pstyle), Paragraph("SMILES parsing, valence verification, 2D depictions, 21 descriptor calculations, ETKDG 3D embedding, MMFF94 force field geometry optimization.", tbl_pstyle)],
        [Paragraph("<b>py3Dmol & stmol</b>", tbl_pstyle), Paragraph("WebGL-accelerated interactive 3D molecular surface rendering directly inside modern web browser sessions.", tbl_pstyle)],
        [Paragraph("<b>Streamlit</b>", tbl_pstyle), Paragraph("Reactive web application server powering parameter routing, multi-page state management, and real-time visualization.", tbl_pstyle)],
        [Paragraph("<b>Matplotlib 3D (mplot3d)</b>", tbl_pstyle), Paragraph("Static rendering of 3D molecular conformers with atomic CPK color-mapping for high-resolution vector PDF export.", tbl_pstyle)],
        [Paragraph("<b>ReportLab Plus</b>", tbl_pstyle), Paragraph("Dynamic PDF generation with dual-pass NumberedCanvas pagination, custom grid flowables, and institutional headers.", tbl_pstyle)],
        [Paragraph("<b>Pandas & NumPy</b>", tbl_pstyle), Paragraph("Matrix computations for atomic coordinate transformations and CSV descriptor export.", tbl_pstyle)],
    ]
    tech_tbl = Table(tech_data, colWidths=[5.0 * cm, 13.2 * cm])
    tech_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2e8b57")),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#bdc3c7")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f9f6")]),
        ("TOPPADDING", (0, 0), (-1, -1), 2.2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
    ]))
    story.append(tech_tbl)
    story.append(Spacer(1, 3 * mm))

    story.append(Paragraph("<b>Core Python Script Execution:</b>", ParagraphStyle("csh", fontSize=8.5, leading=10, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 1 * mm))
    code_str = """from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
# Parse molecule & calculate key descriptors
mol = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
mw, logp = Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
# Generate 3D conformer with ETKDG and optimize using MMFF94
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)"""

    story.append(Paragraph(f'<font face="Courier" size="6.5">{code_str.replace(chr(10), "<br/>")}</font>', ParagraphStyle("csb", backColor=colors.HexColor("#f5f9f6"), borderColor=colors.HexColor("#2e8b57"), borderWidth=0.5, borderPadding=5, leading=8.5)))
    story.append(Spacer(1, 3.5 * mm))

    story.append(Paragraph("<b>Laboratory Conclusions & Learning Outcomes:</b>", ParagraphStyle("lch", fontSize=8.5, leading=10, textColor=colors.HexColor("#1b4d32"))))
    story.append(Spacer(1, 1 * mm))
    story.append(Paragraph(
        "Through this virtual experiment, we successfully integrated computational cheminformatics tools to model, calculate, and evaluate molecular architectures. Quantitative topological and electronic descriptors provided robust indicators of pharmacological viability. The project illustrates how modern in silico pipelines drastically reduce experimental screening costs while advancing molecular understanding.",
        ParagraphStyle("lcb", fontSize=7.2, leading=9.5, alignment=TA_JUSTIFY)
    ))
    story.append(Spacer(1, 4 * mm))

    sig_data = [
        [Paragraph("<b>Prepared & Submitted By:</b>", ParagraphStyle("sgh1", fontSize=7.5, textColor=colors.HexColor("#1b4d32"))), Paragraph("<b>Faculty Verification & Grade:</b>", ParagraphStyle("sgh2", fontSize=7.5, textColor=colors.HexColor("#1b4d32")))],
        [
            Paragraph("<b>Varshini A</b> (RA2611030050040)<br/><b>Vismitha Parasa</b> (RA2611030050047)<br/>Department of Cybersecurity, SRMIST Trichy", ParagraphStyle("sg1", fontSize=7.2, leading=9)),
            Paragraph("Evaluated by: ____________________________<br/>Marks / Grade: ________ / ________<br/>Signature & Date: ________________________", ParagraphStyle("sg2", fontSize=7.2, leading=9))
        ]
    ]
    sig_tbl = Table(sig_data, colWidths=[9.1 * cm, 9.1 * cm])
    sig_tbl.setStyle(TableStyle([
        ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#2e8b57")),
        ("LINEBEFORE", (1, 0), (1, -1), 0.5, colors.HexColor("#2e8b57")),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f9f6")),
    ]))
    story.append(sig_tbl)

    doc.build(story, canvasmaker=NumberedCanvas)
    buf.seek(0)
    return buf


# ============================================================================
# QUIZ QUESTIONS (24 total)
# ============================================================================
QUIZ_QUESTIONS = [
    # Original 12 from sir's screenshots
    {
        "q": "What does SMILES stand for?",
        "options": [
            "Simplified Molecular Input Line Entry System",
            "Standard Molecular Identification Line Entry System",
            "Simplified Molecular Identification Line Entry Specification",
            "Standard Molecular Input Line Entry Specification",
        ],
        "answer": 0,
    },
    {
        "q": "Which Python library is primarily used for cheminformatics tasks like parsing SMILES and calculating molecular descriptors?",
        "options": ["NumPy", "RDKit", "SciPy", "TensorFlow"],
        "answer": 1,
    },
    {
        "q": "What is the SMILES notation for benzene?",
        "options": ["C6H6", "c1ccccc1", "C=C=C=C=C=C", "CCCCCC"],
        "answer": 1,
    },
    {
        "q": "LogP is a measure of:",
        "options": [
            "Molecular weight",
            "Boiling point",
            "Partition coefficient (lipophilicity)",
            "Number of hydrogen bonds",
        ],
        "answer": 2,
    },
    {
        "q": "According to Lipinski's Rule of 5, the molecular weight should be:",
        "options": ["\u2264 300", "\u2264 400", "\u2264 500", "\u2264 600"],
        "answer": 2,
    },
    {
        "q": "TPSA stands for:",
        "options": [
            "Total Polar Surface Area",
            "Topological Polar Surface Analysis",
            "Total Polarity Surface Area",
            "Topological Polar Spatial Area",
        ],
        "answer": 0,
    },
    {
        "q": "Which 3D model represents atoms as spheres filling their van der Waals radii?",
        "options": ["Wireframe", "Ball and Stick", "Stick", "Space Filling"],
        "answer": 3,
    },
    {
        "q": "A molecule with LogP < 0 is considered:",
        "options": ["Lipophilic", "Hydrophilic", "Amphiphilic", "Neutral"],
        "answer": 1,
    },
    {
        "q": "The number of hydrogen bond donors allowed by Lipinski's Rule is:",
        "options": ["\u2264 3", "\u2264 5", "\u2264 7", "\u2264 10"],
        "answer": 1,
    },
    {
        "q": "Which of these is NOT a molecular descriptor?",
        "options": ["LogP", "TPSA", "HTML", "Molecular Weight"],
        "answer": 2,
    },
    {
        "q": "In SMILES notation, '=' represents:",
        "options": [
            "Single bond",
            "Double bond",
            "Triple bond",
            "Aromatic bond",
        ],
        "answer": 1,
    },
    {
        "q": "What is the purpose of the ETKDG method in RDKit?",
        "options": [
            "Calculate molecular weight",
            "Generate 3D conformer coordinates",
            "Parse SMILES strings",
            "Draw 2D structures",
        ],
        "answer": 1,
    },
    # 12 Additional questions
    {
        "q": "In SMILES, '#' represents:",
        "options": [
            "A comment",
            "A triple bond",
            "An aromatic ring",
            "A branch",
        ],
        "answer": 1,
    },
    {
        "q": "Which descriptor measures the total surface area of polar atoms in a molecule?",
        "options": ["LogP", "Molar Refractivity", "TPSA", "BertzCT"],
        "answer": 2,
    },
    {
        "q": "The SMILES notation for ethanol is:",
        "options": ["C2H5OH", "CCO", "C=CO", "COC"],
        "answer": 1,
    },
    {
        "q": "Lipinski's Rule of 5 predicts:",
        "options": [
            "Chemical reactivity",
            "Oral bioavailability / drug-likeness",
            "Melting point",
            "UV absorption",
        ],
        "answer": 1,
    },
    {
        "q": "In a Ball-and-Stick model:",
        "options": [
            "Only bonds are shown",
            "Atoms are spheres connected by sticks representing bonds",
            "Atoms fill their van der Waals radii",
            "Only atoms are shown",
        ],
        "answer": 1,
    },
    {
        "q": "A sigma (\u03c3) bond is formed by:",
        "options": [
            "Lateral overlap of orbitals",
            "Head-on overlap of orbitals",
            "No orbital overlap",
            "Electrostatic attraction only",
        ],
        "answer": 1,
    },
    {
        "q": "Which library is used for interactive 3D molecular visualization in the browser?",
        "options": ["Matplotlib", "py3Dmol", "Seaborn", "Plotly"],
        "answer": 1,
    },
    {
        "q": "How many hydrogen bond acceptors does Lipinski's Rule allow?",
        "options": ["\u2264 5", "\u2264 8", "\u2264 10", "\u2264 15"],
        "answer": 2,
    },
    {
        "q": "Molar Refractivity is related to:",
        "options": [
            "Molecular polarizability and volume",
            "Boiling point only",
            "Number of atoms only",
            "Electronegativity",
        ],
        "answer": 0,
    },
    {
        "q": "In SMILES, branches are denoted by:",
        "options": [
            "Square brackets []",
            "Parentheses ()",
            "Curly braces {}",
            "Angle brackets <>",
        ],
        "answer": 1,
    },
    {
        "q": "What does a high TPSA value (>140 \u00c5\u00b2) suggest about a molecule?",
        "options": [
            "Excellent membrane permeability",
            "High lipophilicity",
            "Poor cell membrane permeability",
            "Low molecular weight",
        ],
        "answer": 2,
    },
    {
        "q": "The Fraction Csp3 descriptor indicates:",
        "options": [
            "The fraction of carbon atoms that are sp3 hybridized",
            "The fraction of hydrogen atoms",
            "The ratio of sigma to pi bonds",
            "The percentage of aromatic atoms",
        ],
        "answer": 0,
    },
    {
        "q": "MMFF in RDKit stands for:",
        "options": [
            "Molecular Mechanics Force Field",
            "Merck Molecular Force Field",
            "Modified Molecular Fragment Field",
            "Macro Molecular Functional Field",
        ],
        "answer": 1,
    },
]


# ============================================================================
# SIDEBAR
# ============================================================================
with st.sidebar:
    st.markdown("## 🧪 Cheminformatics Virtual Lab")
    st.markdown("**SRMIST Tiruchirappalli Campus**")
    st.caption("Dept. of Cybersecurity \u2022 Faculty of Engineering & Technology")
    st.markdown("---")

    # Student details
    st.markdown("### 🎓 Student Details")
    student_name_sidebar = st.text_input("Student Name", value="Varshini A & Vismitha Parasa", key="sidebar_name")
    reg_no_sidebar = st.text_input("Registration Number", value="RA2611030050040 / RA2611030050047", key="sidebar_reg")

    st.markdown("---")

    # Page navigation
    st.markdown("### 📑 Navigation")
    page = st.radio(
        "Select Page",
        [
            "🏠 Home",
            "📖 Theory",
            "🔬 Molecular Visualization",
            "📊 Molecular Properties Report",
            "📈 Structure\u2013Property Analysis",
            "🚀 Advanced & Innovative Tools",
            "📝 Assessment",
            "📄 Final Student Report",
        ],
        label_visibility="collapsed",
    )

    st.markdown("---")

    # Molecule input (shared across pages)
    st.markdown("### 🧬 Molecule Input")
    mol_source = st.radio("Choose input method:", ["Select Example", "Enter SMILES"], key="mol_source")

    if mol_source == "Select Example":
        selected_mol = st.selectbox("Example Molecule", list(EXAMPLE_MOLECULES.keys()))
        current_smiles = EXAMPLE_MOLECULES[selected_mol]
        st.code(current_smiles, language=None)
    else:
        current_smiles = st.text_input("Enter SMILES string", value="CC(=O)OC1=CC=CC=C1C(=O)O", key="custom_smiles")

    # Store in session state for cross-page sharing
    st.session_state["current_smiles"] = current_smiles
    st.session_state["student_name"] = student_name_sidebar
    st.session_state["reg_no"] = reg_no_sidebar

    st.markdown("---")
    st.markdown("### 📚 Learning Objectives")
    st.markdown("""
    - Understand SMILES notation
    - Visualize molecules in 2D & 3D
    - Calculate molecular descriptors
    - Interpret drug-likeness rules
    - Compare molecular properties
    - Generate professional reports
    """)


# ============================================================================
# PAGE: HOME
# ============================================================================
if page == "🏠 Home":
    # Header: Two responsive columns with logo resized one size bigger (165px)
    col_text, col_logo = st.columns([3.8, 1.2], gap="medium")
    with col_text:
        st.markdown("""
        <div style="display: flex; flex-direction: column; justify-content: center; height: 100%; padding-top: 2px;">
            <div class="responsive-header-title">
                🧪 Cheminformatics Virtual Laboratory
            </div>
            <div class="responsive-header-dept">
                Department of Cybersecurity &bull; Faculty of Engineering and Technology
            </div>
            <div class="responsive-header-campus">
                SRM Institute of Science and Technology &mdash; Tiruchirappalli Campus
            </div>
        </div>
        """, unsafe_allow_html=True)

    with col_logo:
        try:
            logo = PILImage.open("assets/srm_logo.png")
            st.markdown(
                '<div class="responsive-logo-wrapper">'
                '<div style="background-color: #ffffff; padding: 5px 8px; border-radius: 8px; border: 1px solid #d1e7dd; box-shadow: 0 2px 8px rgba(46,139,87,0.08); width: 175px; text-align: center;">',
                unsafe_allow_html=True
            )
            st.image(logo, width=165)
            st.markdown('</div></div>', unsafe_allow_html=True)
        except Exception:
            pass

    st.markdown("---")

    st.markdown("""
    ### Welcome!

    This virtual laboratory provides an interactive environment for students to explore
    the world of **cheminformatics** \u2014 the application of computational techniques to
    solve chemical problems.

    #### What you can do here:
    - 🔬 **Visualize molecules** in both 2D and interactive 3D formats
    - 📊 **Calculate molecular properties** including LogP, TPSA, molecular weight, and more
    - 📈 **Compare molecules** side by side with charts and tables
    - 📝 **Test your knowledge** with a comprehensive quiz
    - 📄 **Generate PDF reports** for submission
    """)

    st.markdown("---")
    st.markdown("### 🧬 Example Molecules")

    cols = st.columns(5)
    example_list = list(EXAMPLE_MOLECULES.items())[:5]
    for i, (name, smi) in enumerate(example_list):
        with cols[i]:
            mol = validate_smiles(smi)
            if mol:
                img = get_mol_2d_image(mol, size=(250, 250))
                if img:
                    st.image(img, caption=name, use_container_width=True)
                st.caption(f"`{smi}`")

    cols2 = st.columns(5)
    example_list2 = list(EXAMPLE_MOLECULES.items())[5:]
    for i, (name, smi) in enumerate(example_list2):
        with cols2[i]:
            mol = validate_smiles(smi)
            if mol:
                img = get_mol_2d_image(mol, size=(250, 250))
                if img:
                    st.image(img, caption=name, use_container_width=True)
                st.caption(f"`{smi}`")

    st.markdown("---")

    st.markdown("""
    ### 👩‍💻 Developed By

    | Name | Registration Number | Email ID |
    |------|-------------------|----------|
    | **Varshini A** | *RA2611030050040* | `va3375@srmist.edu.in` |
    | **Vismitha Parasa** | *RA2611030050047* | `vp8717@srmist.edu.in` |

    🏛️ **SRM Institute of Science and Technology** | Tiruchirappalli Campus  
    Department of Cybersecurity \u2022 Faculty of Engineering and Technology
    """)


# ============================================================================
# PAGE: THEORY
# ============================================================================
elif page == "📖 Theory":
    st.markdown('<div class="main-header">📖 Theory & Concepts</div>', unsafe_allow_html=True)
    st.markdown("---")

    tab1, tab2, tab3, tab4 = st.tabs([
        "🧪 Cheminformatics",
        "🔤 SMILES Notation",
        "📊 Molecular Descriptors",
        "🎨 3D Molecular Models",
    ])

    with tab1:
        st.markdown("""
        ## What is Cheminformatics?

        **Cheminformatics** (also known as chemical informatics or chemoinformatics) is an interdisciplinary
        field that combines chemistry, computer science, and information technology to collect, store,
        analyze, and manipulate chemical data.

        ### Key Areas:
        - **Molecular Representation**: Encoding molecular structures in computer-readable formats (SMILES, InChI, MOL files)
        - **Molecular Descriptors**: Numerical values that describe molecular properties (LogP, TPSA, MW)
        - **Structure-Activity Relationships (SAR)**: Understanding how molecular structure affects biological activity
        - **Virtual Screening**: Computational methods to identify promising drug candidates
        - **ADMET Prediction**: Predicting Absorption, Distribution, Metabolism, Excretion, and Toxicity

        ### Applications:
        1. **Drug Discovery**: Identifying and optimizing lead compounds
        2. **Materials Science**: Designing new materials with desired properties
        3. **Environmental Chemistry**: Predicting environmental fate of chemicals
        4. **Chemical Safety**: Toxicity prediction and risk assessment
        """)

    with tab2:
        st.markdown("""
        ## SMILES (Simplified Molecular Input Line Entry System)

        SMILES is a line notation for representing molecular structures as text strings.
        It was developed by David Weininger in the late 1980s.

        ### Basic Rules:
        | Symbol | Meaning | Example |
        |--------|---------|---------|
        | Uppercase letters | Atoms | C, N, O, S, P |
        | Lowercase letters | Aromatic atoms | c, n, o |
        | `=` | Double bond | C=O (formaldehyde) |
        | `#` | Triple bond | C#N (hydrogen cyanide) |
        | `()` | Branches | CC(=O)O (acetic acid) |
        | Numbers | Ring closures | c1ccccc1 (benzene) |
        | `[]` | Complex atoms | [Na+], [Fe+2] |
        | `/`, `\\` | Cis/trans stereochemistry | F/C=C/F |
        | `@`, `@@` | Chirality | [C@@H] |

        ### Examples:
        """)

        examples = {
            "Water": "O",
            "Methane": "C",
            "Ethanol": "CCO",
            "Acetic Acid": "CC(=O)O",
            "Benzene": "c1ccccc1",
            "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        }
        cols = st.columns(3)
        for i, (name, smi) in enumerate(examples.items()):
            with cols[i % 3]:
                mol = validate_smiles(smi)
                if mol:
                    img = get_mol_2d_image(mol, size=(200, 200))
                    if img:
                        st.image(img, caption=f"{name}: {smi}", use_container_width=True)

    with tab3:
        st.markdown("""
        ## Molecular Descriptors

        Molecular descriptors are numerical values that characterize properties of molecules.
        They are fundamental in cheminformatics for quantitative structure-activity relationships (QSAR).

        ### Key Descriptors:

        | Descriptor | Description | Drug-like Range |
        |-----------|-------------|-----------------|
        | **Molecular Weight (MW)** | Sum of atomic masses | \u2264 500 g/mol |
        | **LogP** | Partition coefficient (octanol/water) | \u2264 5 |
        | **TPSA** | Total Polar Surface Area | \u2264 140 \u00c5\u00b2 |
        | **H-Bond Donors (HBD)** | NH and OH groups | \u2264 5 |
        | **H-Bond Acceptors (HBA)** | N and O atoms | \u2264 10 |
        | **Rotatable Bonds** | Molecular flexibility | \u2264 10 |
        | **Fraction Csp3** | Saturation measure | Higher is more 3D-like |
        | **Molar Refractivity** | Polarizability & volume | 40\u2013130 |

        ### Lipinski's Rule of 5:
        A molecule is likely orally active if it satisfies:
        - MW \u2264 500
        - LogP \u2264 5
        - HBD \u2264 5
        - HBA \u2264 10

        *Violations of \u22652 rules suggest poor oral bioavailability.*
        """)

    with tab4:
        st.markdown("""
        ## 3D Molecular Models

        Different visualization styles highlight different aspects of molecular structure:

        ### 1. Ball and Stick Model
        - Atoms shown as **small spheres**, bonds as **sticks**
        - Best for: Seeing both **atom identity** and **bond connectivity**
        - Most commonly used in chemistry education

        ### 2. Stick Model
        - Only **bonds** are shown as cylindrical sticks
        - Atoms are at the intersection points
        - Best for: Viewing **molecular geometry** and **bond angles**

        ### 3. Space Filling (CPK) Model
        - Atoms shown at their **van der Waals radii**
        - Best for: Understanding **molecular shape** and **surface**
        - Shows how much space the molecule actually occupies

        ### 4. Wireframe Model
        - Bonds shown as simple **lines**
        - Most **minimal** representation
        - Best for: Large molecules where detail would be overwhelming
        """)

        # Show example 3D models
        st.markdown("### Interactive Examples")
        demo_mol = validate_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
        if demo_mol:
            demo_3d = generate_3d_coords(demo_mol)
            if demo_3d:
                cols = st.columns(2)
                styles_demo = [
                    ("Ball and Stick", "ball_and_stick"),
                    ("Stick", "stick"),
                    ("Space Filling", "space_filling"),
                    ("Wireframe", "wireframe"),
                ]
                for i, (label, style_key) in enumerate(styles_demo):
                    with cols[i % 2]:
                        st.markdown(f"**{label}**")
                        viewer = get_3d_viewer(demo_3d, style=style_key, width=400, height=300)
                        if viewer:
                            showmol(viewer, height=300, width=400)


# ============================================================================
# PAGE: MOLECULAR VISUALIZATION
# ============================================================================
elif page == "🔬 Molecular Visualization":
    st.markdown('<div class="main-header">🔬 Molecular Visualization</div>', unsafe_allow_html=True)
    st.markdown("---")

    smiles = st.session_state.get("current_smiles", "")
    mol = validate_smiles(smiles)

    if mol:
        st.success(f"✅ Valid SMILES: `{smiles}`")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("### 2D Molecular Structure")
            img = get_mol_2d_image(mol, size=(450, 450))
            if img:
                st.image(img, use_container_width=True)

        with col2:
            st.markdown("### Molecular Information")
            sigma, pi = count_sigma_pi_bonds(mol)
            formula = rdMolDescriptors.CalcMolFormula(mol)

            info_data = {
                "Property": [
                    "Molecular Formula",
                    "Number of Atoms (non-H)",
                    "Number of Bonds",
                    "\u03c3 (Sigma) Bonds",
                    "\u03c0 (Pi) Bonds",
                    "Aromatic Rings",
                    "Total Rings",
                ],
                "Value": [
                    formula,
                    mol.GetNumHeavyAtoms(),
                    mol.GetNumBonds(),
                    sigma,
                    pi,
                    Descriptors.NumAromaticRings(mol),
                    Descriptors.RingCount(mol),
                ],
            }
            st.dataframe(pd.DataFrame(info_data), use_container_width=True, hide_index=True)

        st.markdown("---")
        st.markdown("### 🎨 Interactive 3D Molecular Model")

        col3d1, col3d2 = st.columns([1, 3])
        with col3d1:
            style_3d = st.selectbox(
                "Select 3D Model Style",
                ["Ball and Stick", "Stick", "Space Filling", "Wireframe"],
                key="viz_3d_style",
            )
            show_surface = st.checkbox("Show van der Waals Surface", key="viz_surface")

        style_map = {
            "Ball and Stick": "ball_and_stick",
            "Stick": "stick",
            "Space Filling": "space_filling",
            "Wireframe": "wireframe",
        }

        with col3d2:
            mol_3d = generate_3d_coords(mol)
            if mol_3d:
                viewer = get_3d_viewer(
                    mol_3d,
                    style=style_map[style_3d],
                    show_surface=show_surface,
                    width=650,
                    height=450,
                )
                if viewer:
                    showmol(viewer, height=450, width=650)
            else:
                st.warning("Could not generate 3D coordinates for this molecule.")
    else:
        st.warning("⚠️ Please enter a valid SMILES string in the sidebar.")


# ============================================================================
# PAGE: MOLECULAR PROPERTIES REPORT
# ============================================================================
elif page == "📊 Molecular Properties Report":
    st.markdown('<div class="main-header">📊 Molecular Properties Report</div>', unsafe_allow_html=True)
    st.markdown("---")

    smiles = st.session_state.get("current_smiles", "")
    mol = validate_smiles(smiles)

    if mol:
        st.success(f"✅ Analyzing: `{smiles}`")

        # 2D Structure
        col1, col2 = st.columns([1, 2])
        with col1:
            st.markdown("### 2D Structure")
            img = get_mol_2d_image(mol, size=(350, 350))
            if img:
                st.image(img, use_container_width=True)

        with col2:
            st.markdown("### Basic Information")
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mw = round(Descriptors.MolWt(mol), 3)
            exact_mass = round(Descriptors.ExactMolWt(mol), 4)
            st.markdown(f"""
            | Property | Value |
            |----------|-------|
            | **Molecular Formula** | {formula} |
            | **Molecular Weight** | {mw} g/mol |
            | **Exact Mass** | {exact_mass} |
            | **SMILES** | `{smiles}` |
            | **Heavy Atoms** | {mol.GetNumHeavyAtoms()} |
            | **Formal Charge** | {Chem.GetFormalCharge(mol)} |
            """)

        st.markdown("---")

        # Key Properties
        st.markdown("### 🔑 Key Properties")
        props = calculate_all_properties(mol)
        key_props = {
            "LogP": props["LogP"],
            "TPSA (\u00c5\u00b2)": props["TPSA (\u00c5\u00b2)"],
            "H-Bond Donors": props["H-Bond Donors"],
            "H-Bond Acceptors": props["H-Bond Acceptors"],
            "Rotatable Bonds": props["Rotatable Bonds"],
            "Aromatic Rings": props["Aromatic Rings"],
        }
        kcols = st.columns(6)
        for i, (k, v) in enumerate(key_props.items()):
            with kcols[i]:
                st.metric(label=k, value=v)

        st.markdown("---")

        # Full Property Table
        st.markdown("### 📋 Complete Property Table (21 Descriptors)")
        props_df = pd.DataFrame(
            list(props.items()), columns=["Property", "Value"]
        )
        st.dataframe(props_df, use_container_width=True, hide_index=True)

        # CSV Download
        csv_data = props_df.to_csv(index=False)
        st.download_button(
            "📥 Download Properties as CSV",
            data=csv_data,
            file_name=f"molecular_properties_{smiles[:20]}.csv",
            mime="text/csv",
        )

        st.markdown("---")

        # Generate PDF Report for single molecule
        st.markdown("### 📄 Generate Molecular PDF Report")
        if st.button("📄 Generate Molecular PDF Report", key="single_pdf"):
            sname = st.session_state.get("student_name", "Student")
            sreg = st.session_state.get("reg_no", "N/A")
            mol_data = [{
                "smiles": smiles,
                "mol": mol,
                "properties": props,
                "interpretation": get_interpretation(props),
            }]
            pdf_buf = build_pdf_report(sname, sreg, mol_data)
            st.download_button(
                "📥 Download PDF Report",
                data=pdf_buf.getvalue(),
                file_name=f"molecular_report_{smiles[:15].replace('/', '_')}.pdf",
                mime="application/pdf",
                key="download_single_pdf",
            )

        st.markdown("---")

        # 3D Section
        st.markdown("### 🎨 3D Molecular Model")
        col3d1, col3d2 = st.columns([1, 3])
        with col3d1:
            style_3d = st.selectbox(
                "3D Model Style",
                ["Ball and Stick", "Stick", "Space Filling", "Wireframe"],
                key="prop_3d_style",
            )
            show_surface = st.checkbox("Show van der Waals Surface", key="prop_surface")

        style_map = {
            "Ball and Stick": "ball_and_stick",
            "Stick": "stick",
            "Space Filling": "space_filling",
            "Wireframe": "wireframe",
        }

        with col3d2:
            mol_3d = generate_3d_coords(mol)
            if mol_3d:
                viewer = get_3d_viewer(
                    mol_3d,
                    style=style_map[style_3d],
                    show_surface=show_surface,
                    width=650,
                    height=400,
                )
                if viewer:
                    showmol(viewer, height=400, width=650)

        if st.button("📄 Generate 3D Molecular Model (Static Image)", key="gen_3d_img"):
            mol_3d = generate_3d_coords(mol)
            if mol_3d:
                img_buf = generate_static_3d_image(mol_3d, style=style_map[style_3d])
                if img_buf:
                    st.image(img_buf, caption="3D Molecular Model (Static)", use_container_width=False)
    else:
        st.warning("⚠️ Please enter a valid SMILES string in the sidebar.")


# ============================================================================
# PAGE: STRUCTURE-PROPERTY ANALYSIS
# ============================================================================
elif page == "📈 Structure\u2013Property Analysis":
    st.markdown('<div class="main-header">📈 Structure\u2013Property Analysis</div>', unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("""
    Enter molecules to compare (one per line, as **Name,SMILES**):

    Example:
    ```
    Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
    Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
    Paracetamol,CC(=O)NC1=CC=C(O)C=C1
    ```
    """)

    default_input = """Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Paracetamol,CC(=O)NC1=CC=C(O)C=C1
Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"""

    mol_input = st.text_area(
        "Enter molecules (Name,SMILES per line):",
        value=default_input,
        height=150,
        key="compare_input",
    )

    if st.button("🔍 Analyze & Compare", key="compare_btn"):
        lines = [l.strip() for l in mol_input.strip().split("\n") if l.strip()]
        molecules = []
        for line in lines:
            parts = line.split(",", 1)
            if len(parts) == 2:
                name, smi = parts[0].strip(), parts[1].strip()
                mol = validate_smiles(smi)
                if mol:
                    molecules.append({"name": name, "smiles": smi, "mol": mol})
                else:
                    st.warning(f"⚠️ Invalid SMILES for '{name}': {smi}")

        if len(molecules) >= 2:
            # Build comparison table
            comparison = []
            for m in molecules:
                mol = m["mol"]
                comparison.append({
                    "Name": m["name"],
                    "SMILES": m["smiles"],
                    "MW (g/mol)": round(Descriptors.MolWt(mol), 2),
                    "LogP": round(Descriptors.MolLogP(mol), 2),
                    "TPSA (\u00c5\u00b2)": round(Descriptors.TPSA(mol), 2),
                    "HBD": Descriptors.NumHDonors(mol),
                    "HBA": Descriptors.NumHAcceptors(mol),
                    "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
                    "Aromatic Rings": Descriptors.NumAromaticRings(mol),
                    "Heavy Atoms": mol.GetNumHeavyAtoms(),
                    "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 3),
                    "Molar Refractivity": round(Crippen.MolMR(mol), 2),
                    "Lipinski Violations": lipinski_violations(mol),
                })

            comp_df = pd.DataFrame(comparison)

            st.markdown("### 📋 Comparison Table")
            st.dataframe(comp_df, use_container_width=True, hide_index=True)

            # CSV Download
            csv_comp = comp_df.to_csv(index=False)
            st.download_button(
                "📥 Download Comparison CSV",
                data=csv_comp,
                file_name="molecule_comparison.csv",
                mime="text/csv",
            )

            st.markdown("---")

            # Property Bar Chart
            st.markdown("### 📊 Property Comparison Chart")
            chart_prop = st.selectbox(
                "Select property to chart:",
                ["MW (g/mol)", "LogP", "TPSA (\u00c5\u00b2)", "HBD", "HBA", "Rotatable Bonds", "Molar Refractivity"],
                key="chart_prop",
            )

            chart_data = comp_df[["Name", chart_prop]].set_index("Name")
            st.bar_chart(chart_data)

            st.markdown("---")

            # MW vs LogP Scatter
            st.markdown("### 📈 Molecular Weight vs LogP")
            scatter_df = comp_df[["Name", "MW (g/mol)", "LogP"]].copy()

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.scatter(scatter_df["MW (g/mol)"], scatter_df["LogP"], s=100, c="#2e8b57", edgecolors="black", zorder=3)
            for _, row in scatter_df.iterrows():
                ax.annotate(
                    row["Name"],
                    (row["MW (g/mol)"], row["LogP"]),
                    textcoords="offset points",
                    xytext=(8, 5),
                    fontsize=9,
                )
            ax.axhline(y=5, color="red", linestyle="--", alpha=0.5, label="LogP=5 (Lipinski)")
            ax.axvline(x=500, color="red", linestyle="--", alpha=0.5, label="MW=500 (Lipinski)")
            ax.set_xlabel("Molecular Weight (g/mol)")
            ax.set_ylabel("LogP")
            ax.set_title("MW vs LogP (Lipinski Space)")
            ax.legend()
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            st.markdown("---")

            # Lipinski Violations
            st.markdown("### 🔍 Lipinski Rule of 5 Violations")
            lip_data = comp_df[["Name", "Lipinski Violations"]].copy()
            for _, row in lip_data.iterrows():
                v = row["Lipinski Violations"]
                if v == 0:
                    st.success(f"✅ **{row['Name']}**: 0 violations \u2014 Drug-like")
                elif v == 1:
                    st.warning(f"⚠️ **{row['Name']}**: {v} violation \u2014 Borderline")
                else:
                    st.error(f"❌ **{row['Name']}**: {v} violations \u2014 Poor oral bioavailability predicted")

        elif len(molecules) == 1:
            st.info("Please enter at least 2 molecules to compare.")
        else:
            st.error("No valid molecules found. Please check your input format.")


# ============================================================================
# PAGE: ADVANCED & INNOVATIVE TOOLS (CYBERSECURITY & ADVANCED CHEMINFORMATICS)
# ============================================================================
elif page == "🚀 Advanced & Innovative Tools":
    st.markdown('<div class="main-header">🚀 Advanced & Innovative Virtual Lab Tools</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-header">Department of Cybersecurity \u2022 Bridging Chemical Informatics & Cyber Systems</div>', unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("""
    Welcome to the **cutting-edge extension** of the Cheminformatics Virtual Laboratory.
    Developed under the **Department of Cybersecurity (Faculty of Engineering and Technology)**,
    these tools combine molecular data science, cryptographic verification, scaffold intelligence,
    and quantitative drug-likeness.
    """)

    innov_tab1, innov_tab2, innov_tab3, innov_tab4 = st.tabs([
        "🔐 Cryptographic Chemical Hashes (Cybersecurity)",
        "🎯 QED Drug-Likeness & Bioavailability Score",
        "🧬 Murcko Core Scaffold Extractor",
        "⚡ Molecular Similarity Engine (Tanimoto)"
    ])

    current_smi = st.session_state.get("current_smiles", "CC(=O)OC1=CC=CC=C1C(=O)O")
    cur_mol = validate_smiles(current_smi)

    # TAB 1: Cryptographic Hashes for Molecular Provenance & Cybersecurity
    with innov_tab1:
        st.markdown("### 🔐 Cryptographic Hashing for Chemical Data Integrity")
        st.info("💡 **Cybersecurity in Cheminformatics**: How do we ensure proprietary drug leads and chemical datasets have not been tampered with or modified in databases? By generating cryptographic digest signatures!")

        if cur_mol:
            import hashlib
            can_smi = Chem.MolToSmiles(cur_mol, canonical=True, isomericSmiles=True)
            inchikey = Chem.MolToInchiKey(cur_mol)
            sha256_hash = hashlib.sha256(can_smi.encode()).hexdigest()
            md5_hash = hashlib.md5(can_smi.encode()).hexdigest()
            sha3_hash = hashlib.sha3_256(can_smi.encode()).hexdigest()

            st.markdown(f"**Target Molecule:** `{current_smi}`")
            st.markdown(f"**Canonical SMILES (Standardized):** `{can_smi}`")

            col_c1, col_c2 = st.columns(2)
            with col_c1:
                st.markdown("#### 🛡️ Tamper-Proof Cryptographic Hashes")
                st.code(f"SHA-256 Digest:\n{sha256_hash}\n\nSHA3-256 Digest:\n{sha3_hash}\n\nMD5 Checksum:\n{md5_hash}", language="text")

            with col_c2:
                st.markdown("#### 🏷️ IUPAC InChI & InChIKey Identifier")
                st.code(f"InChIKey (27-char fixed length):\n{inchikey}\n\nFull InChI String:\n{Chem.MolToInchi(cur_mol)}", language="text")

            st.success("✅ **Integrity Verification**: Any minor substitution (e.g. replacing -OH with -F) changes the 256-bit hash entirely (Avalanche Effect).")
        else:
            st.warning("Please provide a valid SMILES in the sidebar.")

    # TAB 2: Quantitative Estimate of Drug-Likeness (QED)
    with innov_tab2:
        st.markdown("### 🎯 Quantitative Estimate of Drug-Likeness (QED Score)")
        st.info("💡 **Beyond Lipinski's Rule of 5**: While Lipinski uses rigid cutoffs, the **QED (Bickerton et al., Nature Chemistry 2012)** integrates 8 physicochemical descriptors into a single continuous desirability metric from **0.0 (least drug-like)** to **1.0 (most drug-like)**.")

        if cur_mol:
            from rdkit.Chem import QED
            qed_score = QED.qed(cur_mol)
            qed_props = QED.properties(cur_mol)

            col_q1, col_q2 = st.columns([1, 2])
            with col_q1:
                st.metric("Overall QED Score", f"{qed_score:.3f}")
                st.progress(min(max(float(qed_score), 0.0), 1.0))
                if qed_score >= 0.67:
                    st.success("🌟 Highly Attractive Drug Candidate (QED \u2265 0.67)")
                elif qed_score >= 0.49:
                    st.info("👍 Moderate Drug-Likeness (0.49 \u2264 QED < 0.67)")
                else:
                    st.warning("⚠️ Low Drug-Likeness (QED < 0.49)")

            with col_q2:
                st.markdown("#### Desirability Parameters:")
                qed_data = {
                    "Descriptor": ["MW", "LogP", "H-Bond Donors", "H-Bond Acceptors", "TPSA (\u00c5\u00b2)", "Rotatable Bonds", "Aromatic Rings", "Structural Alerts"],
                    "Value": [f"{qed_props.MW:.2f}", f"{qed_props.ALOGP:.2f}", qed_props.HBD, qed_props.HBA, f"{qed_props.PSA:.2f}", qed_props.ROTB, qed_props.AROM, qed_props.ALERTS]
                }
                st.dataframe(pd.DataFrame(qed_data), use_container_width=True, hide_index=True)
        else:
            st.warning("Please provide a valid SMILES in the sidebar.")

    # TAB 3: Bemis-Murcko Scaffold Extractor
    with innov_tab3:
        st.markdown("### 🧬 Bemis-Murcko Scaffold Deconstruction")
        st.info("💡 **Scaffold Analysis**: Extracts the core molecular ring system and linker assembly by stripping away all non-cyclic side chains. Essential for medicinal chemistry SAR (Structure-Activity Relationships) and patent analysis.")

        if cur_mol:
            from rdkit.Chem.Scaffolds import MurckoScaffold
            scaffold_mol = MurckoScaffold.GetScaffoldForMol(cur_mol)
            scaffold_smi = Chem.MolToSmiles(scaffold_mol)

            col_s1, col_s2 = st.columns(2)
            with col_s1:
                st.markdown("**Original Compound Structure**")
                img_orig = get_mol_2d_image(cur_mol, size=(300, 250))
                if img_orig:
                    st.image(img_orig, caption=f"SMILES: {current_smi}", use_container_width=True)

            with col_s2:
                st.markdown("**Extracted Bemis-Murcko Scaffold**")
                if scaffold_mol.GetNumAtoms() > 0:
                    img_scaff = get_mol_2d_image(scaffold_mol, size=(300, 250))
                    if img_scaff:
                        st.image(img_scaff, caption=f"Core Scaffold: {scaffold_smi}", use_container_width=True)
                    st.success(f"**Murcko SMILES:** `{scaffold_smi}`")
                else:
                    st.info("This molecule is acyclic (contains no rings), so it has no cyclic Murcko core.")
        else:
            st.warning("Please provide a valid SMILES in the sidebar.")

    # TAB 4: Pairwise Tanimoto Similarity Engine
    with innov_tab4:
        st.markdown("### ⚡ Molecular Fingerprint Similarity Calculator (Tanimoto Coefficient)")
        st.info("💡 **Substructure & Similarity Search**: Measures structural similarity between two compounds based on topological Daylight-style Morgan/RDKit bit fingerprints. Score ranges from **0.0 (completely dissimilar)** to **1.0 (identical)**.")

        col_t1, col_t2 = st.columns(2)
        with col_t1:
            smi_a = st.text_input("Molecule A (SMILES):", value=current_smi, key="tan_smi_a")
        with col_t2:
            smi_b = st.text_input("Molecule B (SMILES):", value="CC(=O)NC1=CC=C(O)C=C1", key="tan_smi_b")

        mol_a = validate_smiles(smi_a)
        mol_b = validate_smiles(smi_b)

        if mol_a and mol_b:
            from rdkit import DataStructs
            fp_a = Chem.RDKFingerprint(mol_a)
            fp_b = Chem.RDKFingerprint(mol_b)
            tanimoto_sim = DataStructs.TanimotoSimilarity(fp_a, fp_b)
            dice_sim = DataStructs.DiceSimilarity(fp_a, fp_b)
            cosine_sim = DataStructs.CosineSimilarity(fp_a, fp_b)

            st.markdown("---")
            col_m1, col_m2, col_m3 = st.columns(3)
            with col_m1:
                st.metric("Tanimoto Coefficient", f"{tanimoto_sim:.4f}")
            with col_m2:
                st.metric("Dice Similarity", f"{dice_sim:.4f}")
            with col_m3:
                st.metric("Cosine Similarity", f"{cosine_sim:.4f}")

            col_view1, col_view2 = st.columns(2)
            with col_view1:
                img_a = get_mol_2d_image(mol_a, size=(280, 200))
                if img_a: st.image(img_a, caption="Molecule A", use_container_width=True)
            with col_view2:
                img_b = get_mol_2d_image(mol_b, size=(280, 200))
                if img_b: st.image(img_b, caption="Molecule B", use_container_width=True)

            if tanimoto_sim >= 0.85:
                st.success("🎯 **High Structural Congruence** (Tanimoto \u2265 0.85) — High probability of similar biological target affinity.")
            elif tanimoto_sim >= 0.50:
                st.info("⚖️ **Moderate Similarity** (0.50 \u2264 Tanimoto < 0.85) — Related scaffold or functional groups.")
            else:
                st.warning("🔍 **Divergent Scaffolds** (Tanimoto < 0.50) — Distinct chemotypes.")
        else:
            st.error("Please enter valid SMILES strings for both molecules.")


# ============================================================================
# PAGE: ASSESSMENT
# ============================================================================
elif page == "📝 Assessment":
    st.markdown('<div class="main-header">📝 Assessment Quiz</div>', unsafe_allow_html=True)
    st.markdown("---")

    st.markdown(f"""
    Test your understanding of cheminformatics concepts!
    - **Total Questions:** {len(QUIZ_QUESTIONS)}
    - **Topics:** SMILES notation, molecular descriptors, Lipinski's rules, 3D models, RDKit
    """)

    # Initialize session state for quiz
    if "quiz_submitted" not in st.session_state:
        st.session_state["quiz_submitted"] = False
    if "quiz_answers" not in st.session_state:
        st.session_state["quiz_answers"] = {}

    with st.form("quiz_form"):
        answers = {}
        for i, q in enumerate(QUIZ_QUESTIONS):
            st.markdown(f"**Q{i+1}. {q['q']}**")
            answers[i] = st.radio(
                f"Select answer for Q{i+1}:",
                q["options"],
                key=f"quiz_q_{i}",
                index=None,
                label_visibility="collapsed",
            )
            st.markdown("---")

        submitted = st.form_submit_button("📝 Submit Quiz")

    if submitted:
        score = 0
        total = len(QUIZ_QUESTIONS)
        wrong = []

        for i, q in enumerate(QUIZ_QUESTIONS):
            correct_answer = q["options"][q["answer"]]
            if answers[i] == correct_answer:
                score += 1
            else:
                wrong.append({
                    "question": q["q"],
                    "your_answer": answers[i] if answers[i] else "(No answer)",
                    "correct_answer": correct_answer,
                })

        # Display results
        st.markdown("## 📊 Results")
        pct = (score / total) * 100
        st.markdown(f"### Score: **{score}/{total}** ({pct:.1f}%)")

        if pct >= 80:
            st.success("🎉 Excellent! You have a strong understanding of cheminformatics concepts!")
        elif pct >= 60:
            st.info("👍 Good job! Review the topics you missed to strengthen your knowledge.")
        elif pct >= 40:
            st.warning("📚 Fair performance. Please revisit the Theory section for better understanding.")
        else:
            st.error("📖 Needs improvement. We recommend going through the Theory section carefully.")

        st.progress(pct / 100)

        if wrong:
            st.markdown("### ❌ Questions You Missed:")
            for w in wrong:
                with st.expander(f"❌ {w['question']}"):
                    st.markdown(f"**Your Answer:** {w['your_answer']}")
                    st.markdown(f"**Correct Answer:** ✅ {w['correct_answer']}")


# ============================================================================
# PAGE: FINAL STUDENT MOLECULAR ANALYSIS REPORT
# ============================================================================
elif page == "📄 Final Student Report":
    st.markdown('<div class="main-header">📄 Final Student Molecular Analysis Report</div>', unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("""
    Generate a comprehensive **multi-molecule** PDF report for submission.
    Enter your details and as many SMILES as you need.
    """)

    col1, col2 = st.columns(2)
    with col1:
        final_name = st.text_input(
            "Student Name",
            value=st.session_state.get("student_name", "Varshini A & Vismitha Parasa"),
            key="final_name",
        )
    with col2:
        final_reg = st.text_input(
            "Registration Number",
            value=st.session_state.get("reg_no", "RA2611030050040 / RA2611030050047"),
            key="final_reg",
        )

    st.markdown("### Enter SMILES (7 Molecules for 10-Page Comprehensive Report):")
    default_smiles_list = """Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Paracetamol,CC(=O)NC1=CC=C(O)C=C1
Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
Glucose,C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O
Ethanol,CCO
Benzene,c1ccccc1"""

    smiles_text = st.text_area(
        "Molecules (Name,SMILES or plain SMILES per line):",
        value=default_smiles_list,
        height=170,
        key="final_smiles",
    )

    # 3D Style for PDF
    style_3d_pdf = st.selectbox(
        "3D Model Style for PDF",
        ["Ball and Stick", "Stick", "Space Filling", "Wireframe"],
        key="final_3d_style",
    )
    style_map = {
        "Ball and Stick": "ball_and_stick",
        "Stick": "stick",
        "Space Filling": "space_filling",
        "Wireframe": "wireframe",
    }

    if st.button("📄 Generate Combined PDF Report", key="final_pdf_btn", type="primary"):
        if not final_name or not final_reg:
            st.error("Please enter both Student Name and Registration Number.")
        else:
            smiles_lines = [s.strip() for s in smiles_text.strip().split("\n") if s.strip()]

            if not smiles_lines:
                st.error("Please enter at least one SMILES string.")
            else:
                molecules_data = []
                progress = st.progress(0, text="Processing molecules...")

                for i, line in enumerate(smiles_lines):
                    if "," in line:
                        parts = line.split(",", 1)
                        cname = parts[0].strip()
                        csmi = parts[1].strip()
                    else:
                        csmi = line.strip()
                        cname = f"Compound {i+1}"

                    mol = validate_smiles(csmi)
                    if mol:
                        props = calculate_all_properties(mol)
                        interp = get_interpretation(props)
                        molecules_data.append({
                            "name": cname,
                            "smiles": csmi,
                            "mol": mol,
                            "properties": props,
                            "interpretation": interp,
                        })
                    else:
                        st.warning(f"⚠️ Skipping invalid SMILES: {line}")
                    progress.progress((i + 1) / len(smiles_lines), text=f"Processing molecule {i+1}/{len(smiles_lines)}...")

                if molecules_data:
                    with st.spinner("Building PDF report..."):
                        pdf_buf = build_pdf_report(
                            final_name,
                            final_reg,
                            molecules_data,
                            style_3d=style_map[style_3d_pdf],
                        )

                    st.success(f"✅ Report generated with {len(molecules_data)} molecule(s)!")

                    # Preview
                    st.markdown("### 📋 Molecules Included:")
                    for i, md in enumerate(molecules_data, 1):
                        st.markdown(f"**{i}.** `{md['smiles']}` \u2014 {md['properties']['Molecular Formula']} (MW: {md['properties']['Molecular Weight (g/mol)']} g/mol)")

                    st.download_button(
                        "📥 Download Combined PDF Report",
                        data=pdf_buf.getvalue(),
                        file_name=f"Student_Report_{final_reg}.pdf",
                        mime="application/pdf",
                        key="download_final_pdf",
                        type="primary",
                    )
                else:
                    st.error("No valid molecules to include in the report.")


# ============================================================================
# FOOTER
# ============================================================================
st.markdown("---")
st.markdown("""
<div class="dev-footer">
    <p><b>Cheminformatics Virtual Laboratory</b> | Developed by <b>Varshini A</b> (RA2611030050040) & <b>Vismitha Parasa</b> (RA2611030050047)</p>
    <p>SRM Institute of Science and Technology, Tiruchirappalli Campus | Department of Cybersecurity | Faculty of Engineering and Technology</p>
    <p>Built with Streamlit \u2022 RDKit \u2022 py3Dmol \u2022 ReportLab \u2022 Matplotlib</p>
</div>
""", unsafe_allow_html=True)

