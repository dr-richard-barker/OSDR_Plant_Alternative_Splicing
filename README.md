# Differential Alternative Splicing in Plant Spaceflight Transcriptomes

A coordinated multi-study alternative splicing (AS) analysis of plant spaceflight transcriptomes. By leveraging 199 samples from 7 Arabidopsis and 1 Brassica studies in the NASA Open Science Data Repository (OSDR), this analysis reveals post-transcriptional splicing networks responding dynamically to microgravity, cosmic radiation, and lunar substrates.

## Interactive GitHub Pages Landing Page

This repository is deployed to GitHub Pages using the **COSE theme/layout**.
You can view the full interactive dashboard, including high-resolution figures, UMAP embeddings, cross-study overlap matrices, and reconstructed splice graphs here:
👉 **[Interactive Landing Page](https://[user].github.io/plant-spaceflight-as/)** (hosted on GitHub Pages)

### Deploying the Landing Page to GitHub Pages

To deploy or update this website under your own GitHub account:

#### Method A: Fully Automated GitHub Actions Deployment (Recommended)
We have pre-configured a GitHub Pages Action in `.github/workflows/pages.yml`. When you push to the `main` or `master` branch:
1. GitHub Actions automatically starts the workflow.
2. It packages only the `docs/` folder.
3. It deploys the compiled landing page styled with the **COSE theme** directly to GitHub Pages.
*Note: Make sure to navigate to **Settings > Pages** in your GitHub repository, and under **Build and deployment**, set the **Source** to **GitHub Actions**.*

#### Method B: Native Deploy via docs/ Folder
If you prefer standard branch-based deployment:
1. Navigate to your repository's **Settings > Pages**.
2. Under **Build and deployment**, set the **Source** to **Deploy from a branch**.
3. Select your active branch (e.g., `main` or `master`), and choose the `/docs` folder from the dropdown menu.
4. Click **Save**.

The `.nojekyll` files included in the root and `/docs` directories bypass default Jekyll compilation, ensuring all custom styling sheets, layout templates, and high-resolution assets render natively with optimal performance.

## Contents

- `docs/` - Complete GitHub Pages static site (styled with the custom BARKER-COSE layout)
  - `docs/index.html` - Interactive landing page with results and summaries
  - `docs/assets/` - Full high-resolution figure suites (SVG + PNG)
  - `docs/cose/` - Barker Lab's custom lightweight styling and site-nav registry
- `multi_study_analysis/` - Multi-study raw feature data and metadata
- `multi_study_analysis_v4/` - Autoritative final multi-study results (7-study pooled matrices, machine learning, cross-study recurrence)
- `osd*_analysis/` - Individual OSD study results (rMATS outputs, NMD predictions, GO enrichment, SpliceGrapher diagrams)
  - `osd314_analysis/` - Root altered gravity vs 1g control
  - `osd120_analysis/` - ISS spaceflight root transcriptomes
  - `osd37_analysis/` - Root spaceflight response across 4 accessions (56 samples)
  - `osd251_analysis/` - Fractional gravity under blue light
  - `osd476_analysis/` - Root transcriptomes of plants grown in authentic Apollo lunar regolith
  - `osd658_analysis/` - Simulated galactic cosmic ray (GCR) radiation response
  - `osd678_analysis/` - ISS spaceflight root transcriptomes with phyD photoreceptor mutants
- `fair_deposit/` - Packaged Zenodo and RO-Crate metadata for deposit
  - `fair_deposit/README.md` - Repository overview and Zenodo guide
  - `fair_deposit/ro-crate-metadata.json` - FAIR-compliant RO-Crate manifest (v1.1)
  - `fair_deposit/zenodo.json` - Zenodo repository publication metadata
  - `fair_deposit/data_dictionary.json` - Splicing-specific column annotations
  - `fair_deposit/CITATION.cff` - Citation metadata (Citation File Format)
  - `fair_deposit/LICENSE` - Code and documentation licenses
- `npj_manuscript.md` - Complete submission-ready manuscript text (Markdown)
- `npj_manuscript.tex` - LaTeX source manuscript formatted in the Nature portfolio template

## Key Scientific Findings

1. **Retained Intron Dominance:** Retained Intron (RI) events represent 57–67% of significant events across all 7 Arabidopsis studies.
2. **NMD-coupled Regulation:** 79–86% of significant space splicing changes are sensitive to nonsense-mediated decay, functioning as "unproductive splicing" switches to post-transcriptionally downregulate stress-responsive transcripts.
3. **Circadian Clock Regulation:** Circadian rhythm regulation is significantly enriched in **all four flight studies** (OSD-314, OSD-120, OSD-37, and OSD-678), highlighting clock splicing as a conserved post-transcriptional target of spaceflight stress.
4. **Stressor-Specific Enrichment:** Growth in Apollo lunar regolith (OSD-476, 306 events) enriches for photosynthesis light reactions and cell cycle; simulated GCR radiation (OSD-658) enriches for phospholipid metabolism; fractional gravity (OSD-251) enriches for superoxide response.
5. **Universal Splicing Biomarker (AT4G17310):** Splicing of *AT4G17310* is significantly altered in **four separate studies** spanning fractional gravity, partial gravity, lunar regolith, and GCR radiation—representing the most robust, stressor-independent splicing target under space stress.
6. **Dominant Batch Effects:** Machine learning classification on the 199-sample × 9,875-gene PSI matrix shows that study-specific batch effects overwhelmingly dominate the splicing profiles (UMAP ARI = 0.995 for study, -0.002 for treatment).

## FAIR Zenodo Deposit

This repository is fully packaged to meet FAIR data principles for its peer-reviewed Zenodo deposit. Refer to `fair_deposit/README.md` for comprehensive metadata profiles, column-by-column schemas, and data dictionaries.

## Authors & Citation

- **Author:** Richard Barker (`@dr-richard-barker`, NASA GeneLab / OSDR)
- **License:** Code and documentation are under the **MIT License**.
- **Citation File:** See `fair_deposit/CITATION.cff` or reference below:

```
Barker, R. (2026). Differential alternative splicing in plant spaceflight transcriptomes: a multi-study rMATS analysis with machine learning and ontology integration. npj Microgravity.
```
