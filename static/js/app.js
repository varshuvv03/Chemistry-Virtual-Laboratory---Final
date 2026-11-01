// ============================================================================
// Cheminformatics Virtual Laboratory Client Application (Non-Streamlit)
// Developed by Varshini A & Vismitha Parasa — SRMIST Trichy
// ============================================================================

const DEFAULT_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"; // Aspirin

// Application State
const state = {
  currentPage: "home",
  studentName: localStorage.getItem("studentName") || "Varshini A & Vismitha Parasa",
  regNo: localStorage.getItem("regNo") || "RA2611030050040 / RA2611030050047",
  studentEmail: localStorage.getItem("studentEmail") || "va3375@srmist.edu.in",
  googleVerified: localStorage.getItem("googleVerified") === "true",
  currentSmiles: localStorage.getItem("currentSmiles") || DEFAULT_SMILES,
  molSource: "example",
  mol3dCache: {},
  theory3dLoaded: false,
  lastQuizResult: null,
};

// DOM ready initialization
document.addEventListener("DOMContentLoaded", () => {
  initSidebar();
  initNavigation();
  initPages();
  syncInputs();
  loadCurrentMolecule();
});

// Initialize Sidebar
function initSidebar() {
  const nameInput = document.getElementById("sidebar-name");
  const regInput = document.getElementById("sidebar-reg");
  const molSelect = document.getElementById("sidebar-example-select");
  const customSmilesInput = document.getElementById("sidebar-custom-smiles");
  const smilesCode = document.getElementById("sidebar-smiles-code");

  nameInput.value = state.studentName;
  regInput.value = state.regNo;
  customSmilesInput.value = state.currentSmiles;
  smilesCode.innerText = state.currentSmiles;

  nameInput.addEventListener("input", (e) => {
    state.studentName = e.target.value;
    localStorage.setItem("studentName", state.studentName);
    const finalName = document.getElementById("final-name");
    if (finalName) finalName.value = state.studentName;
  });

  regInput.addEventListener("input", (e) => {
    state.regNo = e.target.value;
    localStorage.setItem("regNo", state.regNo);
    const finalReg = document.getElementById("final-reg");
    if (finalReg) finalReg.value = state.regNo;
  });

  // Toggle input method (Select Example vs Enter SMILES)
  const btnExample = document.getElementById("btn-source-example");
  const btnCustom = document.getElementById("btn-source-custom");
  const groupExample = document.getElementById("group-source-example");
  const groupCustom = document.getElementById("group-source-custom");

  btnExample.addEventListener("click", () => {
    btnExample.classList.add("active");
    btnCustom.classList.remove("active");
    groupExample.style.display = "block";
    groupCustom.style.display = "none";
    state.molSource = "example";
    updateSmiles(molSelect.value);
  });

  btnCustom.addEventListener("click", () => {
    btnCustom.classList.add("active");
    btnExample.classList.remove("active");
    groupExample.style.display = "none";
    groupCustom.style.display = "block";
    state.molSource = "custom";
    updateSmiles(customSmilesInput.value);
  });

  molSelect.addEventListener("change", (e) => {
    updateSmiles(e.target.value);
  });

  customSmilesInput.addEventListener("input", (e) => {
    updateSmiles(e.target.value);
  });
}

function updateSmiles(smiles) {
  state.currentSmiles = (smiles || "").trim();
  localStorage.setItem("currentSmiles", state.currentSmiles);
  document.getElementById("sidebar-smiles-code").innerText = state.currentSmiles;
  loadCurrentMolecule();
}

// Navigation between 8 pages
function initNavigation() {
  const navItems = document.querySelectorAll(".nav-radio-item");
  navItems.forEach((item) => {
    item.addEventListener("click", () => {
      const pageKey = item.dataset.page;
      showPage(pageKey);
    });
  });
}

function showPage(pageKey) {
  state.currentPage = pageKey;

  // Update navigation highlighting
  document.querySelectorAll(".nav-radio-item").forEach((item) => {
    if (item.dataset.page === pageKey) {
      item.classList.add("active");
    } else {
      item.classList.remove("active");
    }
  });

  // Update page visibility
  document.querySelectorAll(".page-view").forEach((view) => {
    if (view.id === `page-${pageKey}`) {
      view.style.display = "block";
    } else {
      view.style.display = "none";
    }
  });

  // Trigger page-specific loaders
  if (pageKey === "theory" && !state.theory3dLoaded) {
    loadTheory3DModels();
    state.theory3dLoaded = true;
  } else if (pageKey === "visualization") {
    loadVisualizationPage();
  } else if (pageKey === "properties") {
    loadPropertiesPage();
  } else if (pageKey === "advanced") {
    loadAdvancedTools();
  } else if (pageKey === "assessment") {
    renderQuizAuthState();
  } else if (pageKey === "final-report") {
    syncFinalReportInputs();
  }
}

function syncInputs() {
  const finalName = document.getElementById("final-name");
  const finalReg = document.getElementById("final-reg");
  const sideName = document.getElementById("sidebar-name");
  const sideReg = document.getElementById("sidebar-reg");
  const quizName = document.getElementById("quiz-input-name");
  const quizReg = document.getElementById("quiz-input-reg");
  const quizEmail = document.getElementById("quiz-input-email");

  if (finalName) finalName.value = state.studentName;
  if (finalReg) finalReg.value = state.regNo;
  if (sideName) sideName.value = state.studentName;
  if (sideReg) sideReg.value = state.regNo;
  if (quizName) quizName.value = state.studentName;
  if (quizReg) quizReg.value = state.regNo;
  if (quizEmail && state.studentEmail) quizEmail.value = state.studentEmail;
}

function syncFinalReportInputs() {
  syncInputs();
}

// Triggered when SMILES changes
function loadCurrentMolecule() {
  if (state.currentPage === "visualization") {
    loadVisualizationPage();
  } else if (state.currentPage === "properties") {
    loadPropertiesPage();
  } else if (state.currentPage === "advanced") {
    loadAdvancedTools();
  }
}

// ============================================================================
// 3D MOLECULAR VIEWER (3Dmol.js)
// ============================================================================
function render3DModel(containerId, molBlock, style, showSurface) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = "";

  if (!molBlock) {
    container.innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#888;">Generating 3D Coordinates...</div>';
    return;
  }

  try {
    const config = { backgroundColor: "white" };
    const viewer = $3Dmol.createViewer(container, config);
    viewer.addModel(molBlock, "mol");

    if (style === "stick") {
      viewer.setStyle({ stick: { radius: 0.15 } });
    } else if (style === "ball_and_stick") {
      viewer.setStyle({ stick: { radius: 0.1 }, sphere: { scale: 0.25 } });
    } else if (style === "space_filling") {
      viewer.setStyle({ sphere: { scale: 1.0 } });
    } else if (style === "wireframe") {
      viewer.setStyle({ line: { linewidth: 2 } });
    } else {
      viewer.setStyle({ stick: {} });
    }

    if (showSurface) {
      viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.7, color: "lightblue" });
    }

    viewer.zoomTo();
    viewer.render();
  } catch (err) {
    console.error("3Dmol render error:", err);
    container.innerHTML = `<div style="padding:20px;color:red;">Error rendering 3D model: ${err.message}</div>`;
  }
}

async function fetchMol3DBlock(smiles) {
  if (state.mol3dCache[smiles]) {
    return state.mol3dCache[smiles];
  }
  try {
    const res = await fetch("/api/mol3d", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles }),
    });
    const data = await res.json();
    if (data.success && data.mol_block) {
      state.mol3dCache[smiles] = data.mol_block;
      return data.mol_block;
    }
    return null;
  } catch (e) {
    console.error("fetchMol3DBlock error:", e);
    return null;
  }
}

// ============================================================================
// PAGE: HOME & INITIALIZATION
// ============================================================================
function initPages() {
  initHomePage();
  initTheoryTabs();
  initVisualizationControls();
  initPropertiesControls();
  initComparisonPage();
  initAdvancedTabs();
  initQuizPage();
  initFinalReportPage();
}

function initHomePage() {
  const pills = document.querySelectorAll("#mol-category-pills .cat-pill");
  const cards = document.querySelectorAll("#mol-cards-container .enhanced-mol-card");
  const searchInput = document.getElementById("mol-library-search");

  function filterCards() {
    const activePill = document.querySelector("#mol-category-pills .cat-pill.active");
    const cat = activePill ? activePill.dataset.cat : "all";
    const q = (searchInput ? searchInput.value : "").trim().toLowerCase();

    cards.forEach((card) => {
      const cardCat = card.dataset.category;
      const cardName = (card.dataset.name || "").toLowerCase();
      const cardSmiles = (card.dataset.smiles || "").toLowerCase();

      const catMatches = cat === "all" || cardCat === cat;
      const searchMatches = !q || cardName.includes(q) || cardSmiles.includes(q) || (cardCat && cardCat.toLowerCase().includes(q));

      if (catMatches && searchMatches) {
        card.style.display = "flex";
      } else {
        card.style.display = "none";
      }
    });
  }

  pills.forEach((pill) => {
    pill.addEventListener("click", () => {
      pills.forEach((p) => p.classList.remove("active"));
      pill.classList.add("active");
      filterCards();
    });
  });

  if (searchInput) {
    searchInput.addEventListener("input", () => filterCards());
  }

  // Select Molecule Buttons
  document.querySelectorAll(".btn-select-mol").forEach((btn) => {
    btn.addEventListener("click", () => {
      const smi = btn.dataset.smiles;
      const molSelect = document.getElementById("sidebar-example-select");
      if (molSelect) {
        molSelect.value = smi;
      }
      updateSmiles(smi);
      showPage("visualization");
    });
  });
}

function initTheoryTabs() {
  const tabs = document.querySelectorAll("#theory-tabs .tab-btn");
  tabs.forEach((tab) => {
    tab.addEventListener("click", () => {
      tabs.forEach((t) => t.classList.remove("active"));
      tab.classList.add("active");
      const targetPane = tab.dataset.pane;
      document.querySelectorAll("#theory-tabs-content .tab-pane").forEach((p) => {
        if (p.id === targetPane) {
          p.classList.add("active");
        } else {
          p.classList.remove("active");
        }
      });
    });
  });
}

async function loadTheory3DModels() {
  const aspirinSmiles = "CC(=O)OC1=CC=CC=C1C(=O)O";
  const molBlock = await fetchMol3DBlock(aspirinSmiles);
  if (molBlock) {
    render3DModel("theory-3d-ball-stick", molBlock, "ball_and_stick", false);
    render3DModel("theory-3d-stick", molBlock, "stick", false);
    render3DModel("theory-3d-space-filling", molBlock, "space_filling", false);
    render3DModel("theory-3d-wireframe", molBlock, "wireframe", false);
  }
}

// ============================================================================
// PAGE: MOLECULAR VISUALIZATION
// ============================================================================
function initVisualizationControls() {
  const styleSelect = document.getElementById("viz-3d-style");
  const surfaceCheck = document.getElementById("viz-surface-check");

  styleSelect.addEventListener("change", () => updateViz3D());
  surfaceCheck.addEventListener("change", () => updateViz3D());
}

async function loadVisualizationPage() {
  const smiles = state.currentSmiles;
  const validBanner = document.getElementById("viz-valid-banner");
  const img2d = document.getElementById("viz-2d-img");
  const infoTbody = document.getElementById("viz-info-tbody");

  if (!smiles) {
    validBanner.className = "warning-box";
    validBanner.innerHTML = "⚠️ Please enter a valid SMILES string in the sidebar.";
    img2d.src = "";
    infoTbody.innerHTML = "";
    return;
  }

  // Validate and get info
  try {
    const res = await fetch("/api/validate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles }),
    });
    const data = await res.json();
    if (!data.valid) {
      validBanner.className = "warning-box";
      validBanner.innerHTML = `⚠️ Invalid SMILES string: <code>${escapeHtml(smiles)}</code>`;
      img2d.src = "";
      infoTbody.innerHTML = "";
      document.getElementById("viz-3d-container").innerHTML = '<div style="display:flex;align-items:center;justify-content:center;height:100%;color:#888;">Invalid molecule</div>';
      return;
    }

    validBanner.className = "success-box";
    validBanner.innerHTML = `✅ Valid SMILES: <code>${escapeHtml(smiles)}</code>`;

    // 2D Image
    img2d.src = `/api/mol2d?smiles=${encodeURIComponent(smiles)}&w=450&h=450`;

    // Info Table
    infoTbody.innerHTML = `
      <tr><td><b>Molecular Formula</b></td><td>${data.formula}</td></tr>
      <tr><td><b>Number of Atoms (non-H)</b></td><td>${data.heavy_atoms}</td></tr>
      <tr><td><b>Number of Bonds</b></td><td>${data.bonds}</td></tr>
      <tr><td><b>σ (Sigma) Bonds</b></td><td>${data.sigma_bonds}</td></tr>
      <tr><td><b>π (Pi) Bonds</b></td><td>${data.pi_bonds}</td></tr>
      <tr><td><b>Aromatic Rings</b></td><td>${data.aromatic_rings}</td></tr>
      <tr><td><b>Total Rings</b></td><td>${data.total_rings}</td></tr>
    `;

    updateViz3D();
  } catch (e) {
    console.error("loadVisualizationPage error:", e);
  }
}

async function updateViz3D() {
  const smiles = state.currentSmiles;
  const styleSelect = document.getElementById("viz-3d-style");
  const surfaceCheck = document.getElementById("viz-surface-check");
  const style = styleSelect.value;
  const showSurface = surfaceCheck.checked;

  const molBlock = await fetchMol3DBlock(smiles);
  render3DModel("viz-3d-container", molBlock, style, showSurface);
}

// ============================================================================
// PAGE: MOLECULAR PROPERTIES REPORT
// ============================================================================
function initPropertiesControls() {
  const styleSelect = document.getElementById("prop-3d-style");
  const surfaceCheck = document.getElementById("prop-surface-check");
  const staticBtn = document.getElementById("prop-gen-static-btn");

  styleSelect.addEventListener("change", () => updateProp3D());
  surfaceCheck.addEventListener("change", () => updateProp3D());

  staticBtn.addEventListener("click", () => {
    const smiles = state.currentSmiles;
    const style = styleSelect.value;
    const staticImgContainer = document.getElementById("prop-static-img-container");
    staticImgContainer.style.display = "block";
    staticImgContainer.innerHTML = `
      <div style="margin-top:14px; text-align:center;">
        <p style="font-weight:600; margin-bottom:8px; color:#1b4d32;">3D Molecular Model (Static Matplotlib Render)</p>
        <img src="/api/mol3d_static?smiles=${encodeURIComponent(smiles)}&style=${style}&t=${Date.now()}" style="max-width:450px; border:1px solid #d1e7dd; border-radius:6px; box-shadow:0 2px 6px rgba(0,0,0,0.06);" />
      </div>
    `;
  });

  // CSV download button
  document.getElementById("prop-csv-btn").addEventListener("click", () => {
    downloadPropertiesCSV();
  });

  // Single PDF download button
  document.getElementById("prop-pdf-btn").addEventListener("click", () => {
    downloadSinglePDF();
  });
}

async function loadPropertiesPage() {
  const smiles = state.currentSmiles;
  const banner = document.getElementById("prop-banner");
  const img2d = document.getElementById("prop-2d-img");

  if (!smiles) {
    banner.className = "warning-box";
    banner.innerHTML = "⚠️ Please enter a valid SMILES string in the sidebar.";
    return;
  }

  try {
    const res = await fetch("/api/properties", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles }),
    });
    const data = await res.json();
    if (!data.success) {
      banner.className = "warning-box";
      banner.innerHTML = `⚠️ Invalid SMILES string: <code>${escapeHtml(smiles)}</code>`;
      return;
    }

    banner.className = "success-box";
    banner.innerHTML = `✅ Analyzing: <code>${escapeHtml(smiles)}</code>`;

    img2d.src = `/api/mol2d?smiles=${encodeURIComponent(smiles)}&w=350&h=350`;

    const props = data.properties;

    // Basic Info table
    document.getElementById("prop-basic-tbody").innerHTML = `
      <tr><td><b>Molecular Formula</b></td><td>${props["Molecular Formula"]}</td></tr>
      <tr><td><b>Molecular Weight</b></td><td>${props["Molecular Weight (g/mol)"]} g/mol</td></tr>
      <tr><td><b>Exact Mass</b></td><td>${props["Exact Mass"]}</td></tr>
      <tr><td><b>SMILES</b></td><td><code>${escapeHtml(smiles)}</code></td></tr>
      <tr><td><b>Heavy Atoms</b></td><td>${props["Number of Heavy Atoms"]}</td></tr>
      <tr><td><b>Formal Charge</b></td><td>${props["Formal Charge"]}</td></tr>
    `;

    // 6 Key metrics
    document.getElementById("prop-metric-logp").innerText = props["LogP"];
    document.getElementById("prop-metric-tpsa").innerText = props["TPSA (Å²)"];
    document.getElementById("prop-metric-hbd").innerText = props["H-Bond Donors"];
    document.getElementById("prop-metric-hba").innerText = props["H-Bond Acceptors"];
    document.getElementById("prop-metric-rb").innerText = props["Rotatable Bonds"];
    document.getElementById("prop-metric-ar").innerText = props["Aromatic Rings"];

    // Complete 21 Descriptors Table
    let tableHtml = "";
    for (const [propName, propVal] of Object.entries(props)) {
      tableHtml += `<tr><td><b>${escapeHtml(propName)}</b></td><td>${escapeHtml(String(propVal))}</td></tr>`;
    }
    document.getElementById("prop-full-tbody").innerHTML = tableHtml;

    updateProp3D();
  } catch (e) {
    console.error("loadPropertiesPage error:", e);
  }
}

async function updateProp3D() {
  const smiles = state.currentSmiles;
  const style = document.getElementById("prop-3d-style").value;
  const showSurface = document.getElementById("prop-surface-check").checked;
  const molBlock = await fetchMol3DBlock(smiles);
  render3DModel("prop-3d-container", molBlock, style, showSurface);
}

async function downloadPropertiesCSV() {
  const smiles = state.currentSmiles;
  const res = await fetch("/api/properties_csv", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles }),
  });
  if (res.ok) {
    const blob = await res.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `molecular_properties_${smiles.slice(0, 15)}.csv`;
    a.click();
    window.URL.revokeObjectURL(url);
  }
}

async function downloadSinglePDF() {
  const smiles = state.currentSmiles;
  const style = document.getElementById("prop-3d-style").value;
  const btn = document.getElementById("prop-pdf-btn");
  btn.innerText = "⏳ Generating PDF...";
  btn.disabled = true;

  try {
    const res = await fetch("/api/single_pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        student_name: state.studentName,
        reg_no: state.regNo,
        smiles: smiles,
        style_3d: style,
      }),
    });
    if (res.ok) {
      const blob = await res.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `molecular_report_${smiles.slice(0, 15).replace(/\//g, "_")}.pdf`;
      a.click();
      window.URL.revokeObjectURL(url);
    }
  } finally {
    btn.innerText = "📄 Generate Molecular PDF Report";
    btn.disabled = false;
  }
}

// ============================================================================
// PAGE: STRUCTURE-PROPERTY ANALYSIS
// ============================================================================
let comparisonBarChart = null;
let comparisonScatterChart = null;
let lastComparisonData = [];

function initComparisonPage() {
  const analyzeBtn = document.getElementById("compare-analyze-btn");
  const chartPropSelect = document.getElementById("compare-chart-prop");
  const csvBtn = document.getElementById("compare-csv-btn");

  const presetMap = {
    nsaids: `Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Paracetamol,CC(=O)NC1=CC=C(O)C=C1
Ibuprofen,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
Naproxen,COc1ccc2cc(C(C)C(=O)O)ccc2c1
Diclofenac,O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl`,
    cns: `Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Dopamine,NCCc1ccc(O)c(O)c1
Serotonin,NCCc1c[nH]c2ccc(O)cc12
Nicotine,CN1CCCC1c2cccnc2
Fluoxetine,CNCCC(Oc1ccc(C(F)(F)F)cc1)c2ccccc2`,
    antibiotics: `Penicillin G,CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O
Amoxicillin,CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O
Ciprofloxacin,O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O
Sulfamethoxazole,Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1`,
    cardio: `Atenolol,CC(C)NCC(O)COc1ccc(CC(N)=O)cc1
Losartan,CCCC1=NC(=C(N1Cc2ccc(cc2)c3ccccc3c4nnn[nH]4)Cl)CO
Amlodipine,CCOC(=O)C1=C(COCCN)NC(=C(C1c2ccccc2Cl)C(=O)OC)C
Metformin,CN(C)C(=N)NC(=N)N`,
    lipinski: `Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Atorvastatin,CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O
Cholesterol,CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C`,
    diverse: `Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Penicillin G,CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O
Atenolol,CC(C)NCC(O)COc1ccc(CC(N)=O)cc1
Glucose,OCC1OC(O)C(O)C(O)C1O
Cholesterol,CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C
Atorvastatin,CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O
Benzene,c1ccccc1`,
  };

  document.querySelectorAll(".btn-preset-load").forEach((btn) => {
    btn.addEventListener("click", () => {
      const pKey = btn.dataset.preset;
      if (presetMap[pKey]) {
        document.getElementById("compare-input").value = presetMap[pKey];
        runComparison();
      }
    });
  });

  analyzeBtn.addEventListener("click", () => runComparison());
  chartPropSelect.addEventListener("change", () => renderComparisonBarChart());
  csvBtn.addEventListener("click", () => downloadComparisonCSV());
}

async function runComparison() {
  const input = document.getElementById("compare-input").value;
  const resultsArea = document.getElementById("compare-results-area");
  const invalidBox = document.getElementById("compare-invalid-box");

  const res = await fetch("/api/compare", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ input_text: input }),
  });
  const data = await res.json();
  const molecules = data.molecules || [];

  if (data.invalid_lines && data.invalid_lines.length > 0) {
    invalidBox.style.display = "block";
    invalidBox.innerHTML = `⚠️ Invalid SMILES skipped:<br/>` + data.invalid_lines.map(escapeHtml).join("<br/>");
  } else {
    invalidBox.style.display = "none";
  }

  if (molecules.length < 2) {
    resultsArea.style.display = "none";
    alert("Please enter at least 2 valid molecules to compare.");
    return;
  }

  resultsArea.style.display = "block";
  lastComparisonData = molecules;

  // Render Comparison Table
  let tbodyHtml = "";
  molecules.forEach((m) => {
    tbodyHtml += `
      <tr>
        <td><b>${escapeHtml(m.name)}</b></td>
        <td><code>${escapeHtml(m.smiles)}</code></td>
        <td>${m.mw}</td>
        <td>${m.logp}</td>
        <td>${m.tpsa}</td>
        <td>${m.hbd}</td>
        <td>${m.hba}</td>
        <td>${m.rotatable_bonds}</td>
        <td>${m.aromatic_rings}</td>
        <td>${m.heavy_atoms}</td>
        <td>${m.fraction_csp3}</td>
        <td>${m.molar_refractivity}</td>
        <td>${m.lipinski_violations === 0 ? "0 (Pass)" : m.lipinski_violations + " Viol."}</td>
      </tr>
    `;
  });
  document.getElementById("compare-table-tbody").innerHTML = tbodyHtml;

  // Render Charts
  renderComparisonBarChart();
  renderComparisonScatterChart(molecules);

  // Render Lipinski Violations Cards
  let lipCardsHtml = "";
  molecules.forEach((m) => {
    const v = m.lipinski_violations;
    if (v === 0) {
      lipCardsHtml += `<div class="success-box" style="margin:6px 0;">✅ <b>${escapeHtml(m.name)}</b>: 0 violations — Drug-like</div>`;
    } else if (v === 1) {
      lipCardsHtml += `<div class="warning-box" style="margin:6px 0;">⚠️ <b>${escapeHtml(m.name)}</b>: 1 violation — Borderline</div>`;
    } else {
      lipCardsHtml += `<div class="info-box" style="margin:6px 0; border-left-color:#e11d48; background:#fff1f2;">❌ <b>${escapeHtml(m.name)}</b>: ${v} violations — Poor oral bioavailability predicted</div>`;
    }
  });
  document.getElementById("compare-lipinski-cards").innerHTML = lipCardsHtml;
}

function renderComparisonBarChart() {
  const prop = document.getElementById("compare-chart-prop").value;
  const propMap = {
    "MW (g/mol)": "mw",
    "LogP": "logp",
    "TPSA (Å²)": "tpsa",
    "HBD": "hbd",
    "HBA": "hba",
    "Rotatable Bonds": "rotatable_bonds",
    "Molar Refractivity": "molar_refractivity",
  };
  const key = propMap[prop] || "mw";

  const labels = lastComparisonData.map((m) => m.name);
  const values = lastComparisonData.map((m) => m[key]);

  const ctx = document.getElementById("compare-bar-chart").getContext("2d");
  if (comparisonBarChart) comparisonBarChart.destroy();

  comparisonBarChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels: labels,
      datasets: [
        {
          label: prop,
          data: values,
          backgroundColor: "#2e8b57",
          borderColor: "#1b4d32",
          borderWidth: 1,
          borderRadius: 4,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        title: { display: true, text: `Comparison: ${prop}`, color: "#1b4d32", font: { size: 14, weight: "bold" } },
      },
      scales: {
        y: { beginAtZero: key !== "logp" },
      },
    },
  });
}

function renderComparisonScatterChart(molecules) {
  const ctx = document.getElementById("compare-scatter-chart").getContext("2d");
  if (comparisonScatterChart) comparisonScatterChart.destroy();

  const points = molecules.map((m) => ({ x: m.mw, y: m.logp, name: m.name }));

  comparisonScatterChart = new Chart(ctx, {
    type: "scatter",
    data: {
      datasets: [
        {
          label: "Molecules",
          data: points,
          backgroundColor: "#2e8b57",
          borderColor: "#1b4d32",
          borderWidth: 1.5,
          pointRadius: 8,
          pointHoverRadius: 10,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: { display: false },
        title: { display: true, text: "MW vs LogP (Lipinski Chemical Space)", color: "#1b4d32", font: { size: 14, weight: "bold" } },
        tooltip: {
          callbacks: {
            label: function (ctx) {
              const p = ctx.raw;
              return `${p.name}: MW=${p.x}, LogP=${p.y}`;
            },
          },
        },
      },
      scales: {
        x: {
          title: { display: true, text: "Molecular Weight (g/mol)" },
          suggestedMax: 600,
        },
        y: {
          title: { display: true, text: "LogP (Lipophilicity)" },
          suggestedMax: 7,
        },
      },
    },
  });
}

async function downloadComparisonCSV() {
  if (!lastComparisonData.length) return;
  const res = await fetch("/api/compare_csv", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ molecules: lastComparisonData }),
  });
  if (res.ok) {
    const blob = await res.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "molecule_comparison.csv";
    a.click();
    window.URL.revokeObjectURL(url);
  }
}

// ============================================================================
// PAGE: ADVANCED & INNOVATIVE TOOLS
// ============================================================================
function initAdvancedTabs() {
  const tabs = document.querySelectorAll("#advanced-tabs .tab-btn");
  tabs.forEach((tab) => {
    tab.addEventListener("click", () => {
      tabs.forEach((t) => t.classList.remove("active"));
      tab.classList.add("active");
      const targetPane = tab.dataset.pane;
      document.querySelectorAll("#advanced-tabs-content .tab-pane").forEach((p) => {
        if (p.id === targetPane) {
          p.classList.add("active");
        } else {
          p.classList.remove("active");
        }
      });
    });
  });

  document.getElementById("tan-calc-btn").addEventListener("click", () => {
    runSimilarityCalculator();
  });
}

async function loadAdvancedTools() {
  const smilesA = state.currentSmiles;
  const smilesB = document.getElementById("tan-smi-b").value || "CC(=O)NC1=CC=C(O)C=C1";
  document.getElementById("tan-smi-a").value = smilesA;

  try {
    const res = await fetch("/api/advanced", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles_a: smilesA, smiles_b: smilesB }),
    });
    const data = await res.json();

    if (!data.mol_a_valid) {
      document.getElementById("adv-hash-target").innerText = "Invalid molecule";
      return;
    }

    // Hashes
    const h = data.hashes;
    document.getElementById("adv-hash-target").innerText = smilesA;
    document.getElementById("adv-hash-canonical").innerText = h.canonical_smiles;
    document.getElementById("adv-hash-sha256").innerText = h.sha256;
    document.getElementById("adv-hash-sha3").innerText = h.sha3_256;
    document.getElementById("adv-hash-md5").innerText = h.md5;
    document.getElementById("adv-hash-inchikey").innerText = h.inchikey;
    document.getElementById("adv-hash-inchi").innerText = h.inchi;

    // QED
    const q = data.qed;
    document.getElementById("adv-qed-score").innerText = q.score.toFixed(3);
    document.getElementById("adv-qed-progress").style.width = `${Math.min(Math.max(q.score * 100, 0), 100)}%`;

    const qedEval = document.getElementById("adv-qed-eval");
    if (q.score >= 0.67) {
      qedEval.className = "success-box";
      qedEval.innerHTML = "🌟 Highly Attractive Drug Candidate (QED ≥ 0.67)";
    } else if (q.score >= 0.49) {
      qedEval.className = "info-box";
      qedEval.innerHTML = "👍 Moderate Drug-Likeness (0.49 ≤ QED < 0.67)";
    } else {
      qedEval.className = "warning-box";
      qedEval.innerHTML = "⚠️ Low Drug-Likeness (QED < 0.49)";
    }

    document.getElementById("adv-qed-tbody").innerHTML = `
      <tr><td>MW</td><td>${q.mw}</td></tr>
      <tr><td>LogP</td><td>${q.logp}</td></tr>
      <tr><td>H-Bond Donors</td><td>${q.hbd}</td></tr>
      <tr><td>H-Bond Acceptors</td><td>${q.hba}</td></tr>
      <tr><td>TPSA (Å²)</td><td>${q.tpsa}</td></tr>
      <tr><td>Rotatable Bonds</td><td>${q.rotatable_bonds}</td></tr>
      <tr><td>Aromatic Rings</td><td>${q.aromatic_rings}</td></tr>
      <tr><td>Structural Alerts</td><td>${q.alerts}</td></tr>
    `;

    // Scaffold
    document.getElementById("adv-scaff-orig-img").src = `/api/mol2d?smiles=${encodeURIComponent(smilesA)}&w=300&h=250`;
    document.getElementById("adv-scaff-orig-smi").innerText = `SMILES: ${smilesA}`;

    const scaffImg = document.getElementById("adv-scaff-img");
    const scaffBanner = document.getElementById("adv-scaff-banner");
    if (data.scaffold.has_scaffold) {
      scaffImg.src = `/api/mol2d?smiles=${encodeURIComponent(data.scaffold.scaffold_smiles)}&w=300&h=250`;
      scaffImg.style.display = "block";
      scaffBanner.className = "success-box";
      scaffBanner.innerHTML = `<b>Murcko SMILES:</b> <code>${escapeHtml(data.scaffold.scaffold_smiles)}</code>`;
    } else {
      scaffImg.style.display = "none";
      scaffBanner.className = "info-box";
      scaffBanner.innerHTML = "This molecule is acyclic (contains no rings), so it has no cyclic Murcko core.";
    }

    // Similarity
    updateSimilarityView(data);
  } catch (e) {
    console.error("loadAdvancedTools error:", e);
  }
}

async function runSimilarityCalculator() {
  const smilesA = document.getElementById("tan-smi-a").value;
  const smilesB = document.getElementById("tan-smi-b").value;

  const res = await fetch("/api/advanced", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ smiles_a: smilesA, smiles_b: smilesB }),
  });
  const data = await res.json();
  updateSimilarityView(data);
}

function updateSimilarityView(data) {
  const smiA = document.getElementById("tan-smi-a").value;
  const smiB = document.getElementById("tan-smi-b").value;

  if (!data.mol_a_valid || !data.mol_b_valid || !data.similarity) {
    alert("Please provide valid SMILES for both molecules.");
    return;
  }

  const s = data.similarity;
  document.getElementById("tan-metric-tanimoto").innerText = s.tanimoto.toFixed(4);
  document.getElementById("tan-metric-dice").innerText = s.dice.toFixed(4);
  document.getElementById("tan-metric-cosine").innerText = s.cosine.toFixed(4);

  document.getElementById("tan-img-a").src = `/api/mol2d?smiles=${encodeURIComponent(smiA)}&w=280&h=200`;
  document.getElementById("tan-img-b").src = `/api/mol2d?smiles=${encodeURIComponent(smiB)}&w=280&h=200`;

  const interpBox = document.getElementById("tan-interpretation");
  if (s.tanimoto >= 0.85) {
    interpBox.className = "success-box";
    interpBox.innerHTML = "🎯 <b>High Structural Congruence</b> (Tanimoto ≥ 0.85) — High probability of similar biological target affinity.";
  } else if (s.tanimoto >= 0.50) {
    interpBox.className = "info-box";
    interpBox.innerHTML = "⚖️ <b>Moderate Similarity</b> (0.50 ≤ Tanimoto < 0.85) — Related scaffold or functional groups.";
  } else {
    interpBox.className = "warning-box";
    interpBox.innerHTML = "🔍 <b>Divergent Scaffolds</b> (Tanimoto < 0.50) — Distinct chemotypes.";
  }
}

// ============================================================================
// PAGE: ASSESSMENT (24 QUESTIONS & GOOGLE VERIFICATION)
// ============================================================================
function initQuizPage() {
  const submitBtn = document.getElementById("quiz-submit-btn");
  const btnVerify = document.getElementById("quiz-btn-verify-signin");
  const btnVarshini = document.getElementById("quiz-btn-fill-varshini");
  const btnSwitch = document.getElementById("quiz-btn-switch-account");
  const btnCert = document.getElementById("quiz-btn-download-cert");
  const btnRetake = document.getElementById("quiz-btn-retake");

  renderQuizAuthState();
  initGoogleIdentityServices();

  if (btnVerify) {
    btnVerify.addEventListener("click", () => {
      const name = (document.getElementById("quiz-input-name").value || "").trim();
      const reg = (document.getElementById("quiz-input-reg").value || "").trim();
      const email = (document.getElementById("quiz-input-email").value || "").trim();

      if (!name || !reg || !email) {
        alert("Please enter your Student Name, Registration Number, and Google Account Email.");
        return;
      }
      if (!email.includes("@")) {
        alert("Please enter a valid Google email address.");
        return;
      }

      setVerifiedStudent(name, reg, email);
    });
  }

  if (btnVarshini) {
    btnVarshini.addEventListener("click", () => {
      setVerifiedStudent("Varshini A", "RA2611030050040", "va3375@srmist.edu.in");
    });
  }

  if (btnSwitch) {
    btnSwitch.addEventListener("click", () => {
      state.googleVerified = false;
      localStorage.removeItem("googleVerified");
      renderQuizAuthState();
    });
  }

  if (submitBtn) {
    submitBtn.addEventListener("click", () => submitQuiz());
  }

  if (btnCert) {
    btnCert.addEventListener("click", () => downloadAssessmentCertificate());
  }

  if (btnRetake) {
    btnRetake.addEventListener("click", () => {
      document.querySelectorAll("#quiz-form input[type='radio']").forEach((r) => (r.checked = false));
      document.getElementById("quiz-results-area").style.display = "none";
      const qCard = document.getElementById("quiz-auth-card");
      if (qCard) qCard.scrollIntoView({ behavior: "smooth" });
    });
  }
}

function initGoogleIdentityServices() {
  if (typeof google !== "undefined" && google.accounts && google.accounts.id) {
    try {
      google.accounts.id.initialize({
        client_id: "928374928374-placeholder.apps.googleusercontent.com",
        callback: handleGoogleCredentialResponse,
        auto_select: false,
      });
      const mount = document.getElementById("g_id_signin_mount");
      if (mount) {
        google.accounts.id.renderButton(mount, {
          theme: "outline",
          size: "medium",
          text: "signin_with",
          shape: "pill",
        });
      }
    } catch (e) {
      console.log("GIS initialization note:", e);
    }
  }
}

function handleGoogleCredentialResponse(response) {
  try {
    const payload = parseJwt(response.credential);
    const name = payload.name || "Google Student";
    const email = payload.email || "";
    const reg = (document.getElementById("quiz-input-reg") ? document.getElementById("quiz-input-reg").value.trim() : "") || state.regNo || "RA2611030050040";
    setVerifiedStudent(name, reg, email);
  } catch (err) {
    console.error("JWT parse error:", err);
  }
}

function parseJwt(token) {
  try {
    const base64Url = token.split(".")[1];
    const base64 = base64Url.replace(/-/g, "+").replace(/_/g, "/");
    const jsonPayload = decodeURIComponent(
      atob(base64)
        .split("")
        .map((c) => "%" + ("00" + c.charCodeAt(0).toString(16)).slice(-2))
        .join("")
    );
    return JSON.parse(jsonPayload);
  } catch (e) {
    return {};
  }
}

function setVerifiedStudent(name, reg, email) {
  state.studentName = name;
  state.regNo = reg;
  state.studentEmail = email;
  state.googleVerified = true;

  localStorage.setItem("studentName", name);
  localStorage.setItem("regNo", reg);
  localStorage.setItem("studentEmail", email);
  localStorage.setItem("googleVerified", "true");

  syncInputs();
  renderQuizAuthState();
}

function renderQuizAuthState() {
  const unverifiedBox = document.getElementById("quiz-auth-unverified");
  const verifiedBanner = document.getElementById("quiz-auth-verified");

  if (!unverifiedBox || !verifiedBanner) return;

  if (state.googleVerified && state.studentEmail) {
    unverifiedBox.style.display = "none";
    verifiedBanner.style.display = "block";

    const dName = document.getElementById("quiz-display-name");
    const dReg = document.getElementById("quiz-display-reg");
    const dEmail = document.getElementById("quiz-display-email");
    const dRecipient = document.getElementById("quiz-display-recipient");
    const avatar = document.getElementById("quiz-avatar-circle");

    if (dName) dName.innerText = state.studentName;
    if (dReg) dReg.innerText = state.regNo;
    if (dEmail) dEmail.innerText = state.studentEmail;
    if (dRecipient) dRecipient.innerText = state.studentEmail;
    if (avatar) avatar.innerText = (state.studentName || "S").charAt(0).toUpperCase();
  } else {
    unverifiedBox.style.display = "block";
    verifiedBanner.style.display = "none";

    const qName = document.getElementById("quiz-input-name");
    const qReg = document.getElementById("quiz-input-reg");
    const qEmail = document.getElementById("quiz-input-email");

    if (qName) qName.value = state.studentName;
    if (qReg) qReg.value = state.regNo;
    if (qEmail) qEmail.value = state.studentEmail;
  }
}

async function submitQuiz() {
  if (!state.googleVerified || !state.studentEmail) {
    alert("⚠️ Please authenticate with your Google Account and enter your Registration Number before submitting the quiz.");
    const authCard = document.getElementById("quiz-auth-card");
    if (authCard) {
      authCard.scrollIntoView({ behavior: "smooth" });
      authCard.style.outline = "2px solid #2e8b57";
      setTimeout(() => (authCard.style.outline = "none"), 2500);
    }
    return;
  }

  const answers = {};
  for (let i = 0; i < 24; i++) {
    const selected = document.querySelector(`input[name="quiz_q_${i}"]:checked`);
    if (selected) {
      answers[String(i)] = parseInt(selected.value);
    }
  }

  const submitBtn = document.getElementById("quiz-submit-btn");
  if (submitBtn) {
    submitBtn.disabled = true;
    submitBtn.innerHTML = "⏳ Submitting &amp; Dispatching Records...";
  }

  try {
    const res = await fetch("/api/quiz/submit", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        answers: answers,
        student_name: state.studentName,
        reg_no: state.regNo,
        student_email: state.studentEmail,
        google_verified: true,
      }),
    });
    const data = await res.json();
    state.lastQuizResult = data;

    const resultsArea = document.getElementById("quiz-results-area");
    resultsArea.style.display = "block";

    document.getElementById("quiz-score-text").innerText = `${data.score}/${data.total} (${data.percentage}%)`;
    document.getElementById("quiz-progress-bar").style.width = `${data.percentage}%`;

    const feedbackBox = document.getElementById("quiz-feedback-box");
    const pct = data.percentage;
    if (pct >= 80) {
      feedbackBox.className = "success-box";
      feedbackBox.innerHTML = `🎉 <b>Outstanding Performance!</b> Score: <b>${data.score}/24 (${pct}%)</b>.<br/>Official scorecard has been sent to <b>${escapeHtml(data.student_email)}</b> and examination records forwarded to <b>Varshini A (va3375@srmist.edu.in)</b>.`;
    } else if (pct >= 60) {
      feedbackBox.className = "info-box";
      feedbackBox.innerHTML = `👍 <b>Good Job!</b> Score: <b>${data.score}/24 (${pct}%)</b>.<br/>Scorecard delivered to <b>${escapeHtml(data.student_email)}</b>; examination record archived with <b>Varshini A (va3375@srmist.edu.in)</b>.`;
    } else {
      feedbackBox.className = "warning-box";
      feedbackBox.innerHTML = `📚 <b>Assessment Completed.</b> Score: <b>${data.score}/24 (${pct}%)</b>.<br/>Results sent to <b>${escapeHtml(data.student_email)}</b> and logged for <b>Varshini A (va3375@srmist.edu.in)</b>.`;
    }

    // Update Email Status Badges
    const studentBadge = document.getElementById("quiz-student-email-result");
    if (studentBadge) {
      studentBadge.innerText = `${data.student_name} <${data.student_email}>`;
    }

    const missedContainer = document.getElementById("quiz-missed-list");
    if (data.wrong && data.wrong.length > 0) {
      document.getElementById("quiz-missed-header").style.display = "block";
      missedContainer.innerHTML = data.wrong
        .map(
          (w) => `
          <div class="missed-question-box">
            <b>❌ Q${w.index}. ${escapeHtml(w.question)}</b><br/>
            <span style="color:#64748b;">Your Answer:</span> ${escapeHtml(w.your_answer)}<br/>
            <span style="color:#1b4d32; font-weight:600;">Correct Answer:</span> ✅ ${escapeHtml(w.correct_answer)}
          </div>
        `
        )
        .join("");
    } else {
      document.getElementById("quiz-missed-header").style.display = "none";
      missedContainer.innerHTML = "";
    }

    resultsArea.scrollIntoView({ behavior: "smooth" });
  } catch (err) {
    console.error("Quiz submission error:", err);
    alert("Submission error: " + err.message);
  } finally {
    if (submitBtn) {
      submitBtn.disabled = false;
      submitBtn.innerHTML = "📝 Submit Quiz &amp; Dispatch Results";
    }
  }
}

async function downloadAssessmentCertificate() {
  const result = state.lastQuizResult || {
    score: 24,
    total: 24,
    percentage: 100,
    student_name: state.studentName,
    reg_no: state.regNo,
    student_email: state.studentEmail,
  };

  const btn = document.getElementById("quiz-btn-download-cert");
  const origText = btn ? btn.innerHTML : "";
  if (btn) btn.innerHTML = "⏳ Generating Official Certificate PDF...";

  try {
    const res = await fetch("/api/quiz/certificate_pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        student_name: state.studentName,
        reg_no: state.regNo,
        student_email: state.studentEmail,
        score: result.score,
        total: result.total,
        percentage: result.percentage,
      }),
    });

    if (!res.ok) throw new Error("Failed to generate certificate PDF");

    const blob = await res.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `Assessment_Certificate_${state.studentName.replace(/\s+/g, "_")}.pdf`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(url);
  } catch (e) {
    alert("Could not download certificate: " + e.message);
  } finally {
    if (btn) btn.innerHTML = origText || "🏆 Download Official PDF Certificate & Scorecard";
  }
}

// ============================================================================
// PAGE: FINAL STUDENT REPORT (10-PAGE PDF)
// ============================================================================
function initFinalReportPage() {
  const genBtn = document.getElementById("final-pdf-btn");
  genBtn.addEventListener("click", () => generateFinalPDFReport());
}

async function generateFinalPDFReport() {
  const sName = document.getElementById("final-name").value.trim();
  const sReg = document.getElementById("final-reg").value.trim();
  const smilesText = document.getElementById("final-smiles").value.trim();
  const style3d = document.getElementById("final-3d-style").value;
  const statusMsg = document.getElementById("final-pdf-status");
  const genBtn = document.getElementById("final-pdf-btn");

  if (!sName || !sReg) {
    alert("Please enter both Student Name and Registration Number.");
    return;
  }
  if (!smilesText) {
    alert("Please enter at least one SMILES string.");
    return;
  }

  genBtn.disabled = true;
  genBtn.innerText = "⏳ Building 10-Page Laboratory Report...";
  statusMsg.style.display = "block";
  statusMsg.className = "info-box";
  statusMsg.innerHTML = "⏳ Processing 7 molecular conformers and compiling full 10-page institutional report...";

  try {
    const res = await fetch("/api/final_pdf", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        student_name: sName,
        reg_no: sReg,
        smiles_text: smilesText,
        style_3d: style3d,
      }),
    });

    if (!res.ok) {
      const err = await res.text();
      statusMsg.className = "warning-box";
      statusMsg.innerHTML = `⚠️ Failed to build PDF: ${escapeHtml(err)}`;
      return;
    }

    const blob = await res.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `Student_Report_${sReg.replace(/[\/\s\\]/g, "_")}.pdf`;
    a.click();
    window.URL.revokeObjectURL(url);

    statusMsg.className = "success-box";
    statusMsg.innerHTML = "✅ Combined 10-Page Student Molecular Analysis Report generated and downloaded successfully!";
  } catch (err) {
    console.error("PDF generation error:", err);
    statusMsg.className = "warning-box";
    statusMsg.innerHTML = `⚠️ Error generating report: ${err.message}`;
  } finally {
    genBtn.disabled = false;
    genBtn.innerText = "📄 Generate Combined PDF Report";
  }
}

// Helper: Escape HTML
function escapeHtml(str) {
  if (str === null || str === undefined) return "";
  return String(str)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#039;");
}

