# AutoGrow4 Modification: Dual-Energy Scoring for Molecular Glues

**What this is:** A fork of AutoGrow4 that generates small molecules bridging protein-protein interfaces (e.g., CRBN-CK1α), with a custom scoring function.

**Core modification:**
- **Input:** Takes *both* a docked protein-protein structure (e.g., CRBN-CK1α complex) *and* the apo protein structure (CRBN alone).
- **Scoring:** New function combines binding energies from **both** structures to prioritize candidates that (1) fit the interface and (2) don’t destabilize the apo state.
- **Inspiration:** Methodology conceptually inspired by [Chem. Sci., 2019, 10, 2869], but independently implemented and adapted.

**Status:** Functional prototype. Generated 20+ candidate molecules for CRBN-CK1α system. No experimental validation. Code is research-grade, not production-ready.

<img width="265" height="189" alt="image" src="https://github.com/user-attachments/assets/81560502-3bd5-498c-bfbf-e9af12245d99" />

### Computational Pipeline Overview

The complete design workflow consists of four phases:

1. **Rosetta Relaxation**: Energy minimization of apo-CRBN and apo-CK1α structures
2. **Protein-Protein Docking**: Global docking search to generate CRBN-CK1α complex poses  
3. **Interface Selection**: Clustering and selection of top energetically favorable PPI pose
4. **De Novo Generation**: Fragment-based growth of candidate stabilizers at the selected interface using a custom fitness function

**Current Repository**: This repo contains **Phase 4 only**—the AutoGrow4 molecular generation engine with custom thermodynamic scoring.

---

### De Novo Generation with Cooperative Scoring

The AutoGrow4 implementation uses a **custom fitness function** that evaluates candidates based on their predicted **cooperative stabilization potential**, inspired by the thermodynamic framework from Brunsveld *et al.*, *Chem. Sci.*, 2019.

#### Thermodynamic Model Implementation

# Custom fitness function (from this repo: vina.py)

r = 0.001987      # Gas constant kcal/mol·K
t = 298.0         # Temperature (K)

kd1 = np.exp(affinity / (r*t))           # KD_II: affinity for apo-CRBN

kd2 = np.exp(affinity_complex / (r*t))   # KD_II/α: affinity for CRBN-CK1α complex

alpha = kd2 / kd1                         # Cooperativity factor (α < 1 = positive cooperativity)

# Scoring prioritizes high affinity to complex AND strong cooperativity
score = -((-np.log(kd2)) * (1 - 1/(alpha + 1)))
