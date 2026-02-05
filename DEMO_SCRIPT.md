# BioScribe Hackathon Demo Script

## 🎯 2-3 Minute Live Demo Script for Judges

---

## 1️⃣ Opening (15 seconds)

**Say:** "BioScribe accelerates drug discovery from **years to days** using AI. Let me show you how."

**Do:** Point to the hackathon badge on the landing page.

---

## 2️⃣ Problem Statement (20 seconds)

**Say:** "Drug discovery today takes 10-15 years and costs $2.6 billion with a 90% failure rate. We're changing that."

**Do:** Click "Try Live Demo" to enter the platform.

---

## 3️⃣ Complete Pipeline Demo (60 seconds)

### Step 1: Enter Protein Sequence
**Say:** "Let's discover drugs for SARS-CoV-2. I'll paste the protease sequence here."

**Do:** Navigate to `/pipeline/new`, paste example sequence:
```
SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ
```

### Step 2: AI Analysis Results
**Say:** "In seconds, our AI identifies binding sites, druggability scores, and molecular properties."

**Do:** Show the analysis results with animated visualizations.

### Step 3: Drug Candidate Generation
**Say:** "Now watch as we generate novel drug candidates optimized for this target."

**Do:** Run drug generation, highlight the 3D molecular viewer.

### Step 4: Docking Visualization
**Say:** "This 3D visualization shows exactly how our top candidate binds to the target protein."

**Do:** Rotate the 3D molecular docking view.

---

## 4️⃣ Clinical Safety Suite (30 seconds)

**Say:** "What sets us apart is our comprehensive safety analysis."

**Do:** Navigate to `/clinical-safety`, show:
- **Drug Interactions** - CYP450 matrix
- **Adverse Events** - Organ toxicity map
- **Pharmacogenomics** - Genetic variant dosing

**Say:** "This is real clinical decision support that reduces patient harm."

---

## 5️⃣ Responsible AI Showcase (20 seconds)

**Say:** "Our platform is built on responsible AI principles."

**Do:** Point to:
- ✅ Disclaimer banner at top
- ✅ XAI explanations for every prediction
- ✅ FAIR data compliance
- ✅ Blockchain audit trails

**Say:** "Every prediction is explainable, traceable, and validated."

---

## 6️⃣ Closing (15 seconds)

**Say:** "BioScribe: From protein to drug candidate in minutes, not years. Questions?"

---

## 🎬 Demo Shortcuts

| Feature | URL |
|---------|-----|
| Landing Page | http://localhost:3000/landing |
| Dashboard | http://localhost:3000/dashboard |
| New Pipeline | http://localhost:3000/pipeline/new |
| Clinical Safety | http://localhost:3000/clinical-safety |
| AtomNet XAI | http://localhost:3000/atomnet |

---

## ⚠️ Key Talking Points

1. **Healthcare Workflow Impact** - End-to-end drug discovery pipeline
2. **Responsible AI** - XAI, FAIR data, blockchain audit trails
3. **Clinical Relevance** - DDI checker, pharmacogenomics, safety signals
4. **Scientific Accuracy** - Real algorithms (BioPython, RDKit, ESOL)
5. **Hackathon Compliance** - Synthetic data only, clear disclaimers
