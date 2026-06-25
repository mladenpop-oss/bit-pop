# Bit-Pop & RNA-Pop — Context za novu sesiju
*Mladen Popović — June 2026*

---

## Tko sam ja
- CNC programer, 22 godine iskustva, Karl Widenmann GmbH, Gerstetten, Njemačka
- Email: mladenpop@gmail.com
- Projekti razvijeni u slobodno vrijeme, na privatnom računalu

---

## Bit-Pop

**Što je:** Open-source multi-genome DNA read classifier u Rustu
**GitHub:** https://github.com/mladenpop-oss/bit-pop
**DOI:** 10.5281/zenodo.20043593

### Finalni benchmark rezultati

**CAMI Low (61 genom, ~1M reads, Illumina 150bp):**
| Config | Mapping | Accuracy |
|--------|---------|----------|
| k=13 + EM t=1.0 | 70.0% | 92.29% strain / ~100% species |
| k12-k15 consensus + EM | 91.1% | 90.07% |
| k13+k22 + EM | 99.48% | 89.86% |

- Species-level accuracy: ~100% u svim konfiguracijama
- Sve greške su intra-clade (nikad cross-species)
- Index: 1.4GB, build time: 17s, runtime: ~136s za 1M reads

**PacBio HiFi (69 genoma, 86k reads, simulirani):**
- Mapping rate: 99.0%
- Accuracy: 95.2%
- Runtime: ~8 minuta, laptop

**Ebola ONT benchmark (klinički metagenomics, 70% human, 15% bakt, 10% Ebola):**
| Config | Mapped | Overall | Ebola | FPR |
|--------|--------|---------|-------|-----|
| Default | 5,576 | 81.5% | 76.5% | 2.49% |
| k21, chunk 125/130 + EM | 10,012 | 99.3% | 99.5% | 0.96% |

- Zero cross-species misclassification u svim testovima
- Ebola→Human: 0%, Bact→Ebola: 0%
- 100 konzervativnih human reads uvijek mapira na Bundibugyo (EVE hipoteza)
- R9.4 i R10.4 kemija testirana

**Što postoji:**
- Desktop GUI (Tauri + Svelte, napravljen za 1 dan)
- Android app (JNI + Kotlin, ConCon multi-index, radi na telefonu)
- Docker + Bioconda
- 312+ unit testova

---

## RNA-Pop

**Što je:** RNA-seq quantification na FM-index engineu
**GitHub:** https://github.com/mladenpop-oss/RNA-Pop
**DOI:** 10.5281/zenodo.20578611

### Finalni benchmark rezultati (2M reads, 11,567 human chr19 transkripata):

| Metric | RNA-Pop | Salmon |
|--------|---------|--------|
| Spearman ρ | 0.993 | 0.989 |
| Pearson r | 0.979 | 0.973 |
| Mapping rate | 97.1% | - |
| Top-10 overlap | 7/10 | 7/10 |
| Speed | 20,500 reads/sec | - |

**Features:**
- Cancer biomarker panels (breast, lung, prostate, colorectal, pancancer)
- Differential expression analysis
- Fusion gene detection
- Clinical HTML reports
- QC metrics
- Paired-end support

---

## Akademski suradnici

### Paul Kagame MSc. (Ruanda)
- Svježe obranjen PhD (phylogenomics)
- Co-author na Bit-Pop paperu
- Status: čita draft u Google Docs
- Target: Oxford Bioinformatics, submission ~kolovoz 2026
- Kontaktiran oko Ebola outbreaka — nema direktnih veza s response timom
- Traži tjedni update, ne dnevni

### Animesh Sharma (NTNU Norway)
- Bioinformatics engineer
- Forkao Bit-Pop → **bit-pep** za proteomiku
  - 574K proteina, 123M peptida, 3 sekunde za index
  - Substack članak objavljen
- Radi na 6-frame translation Ebola genoma (viral proteomics)
- Poly-A slippage problem: 141 poly-A traktova u Ebola, ±3 fuzzy = 987 varijanti — upravljivo
- Na putu 2 tjedna (Swiss-France meetings), vraća se sljedeći vikend
- Rekao: "your work is impressive and will clearly be noticed"
- Zainteresiran za RNA-Pop i viral proteomics suradnju

---

## Status papera

- Bit-Pop draft postoji (paper_new_v2.md / .docx)
- Paul čita, čeka feedback, poslat će komentare
- RNA-Pop Zenodo preprint: DOI 10.5281/zenodo.20578611
- Redoslijed: Bit-Pop paper prvo → RNA-Pop paper → eventualno klinički metagenomics paper

---

## Pending akcije

- [ ] Bit-Pop paper revision s Paulom (Oxford Bioinformatics, kolovoz 2026)
- [ ] RNA-Pop paper (nakon Bit-Pop)
- [ ] Animesh: viral proteomics (6-frame translation + bit-pep + RNA-Pop)
- [ ] Nick Loman (ARTIC Network): poslan email s Ebola benchmark, čeka odgovor
- [ ] ARTIC forum: čeka approval

## Completed
- [x] Bit-Pop Ebola benchmark (99.5% sensitivity, 0.96% FPR)
- [x] RNA-Pop Zenodo: DOI 10.5281/zenodo.20578611
- [x] Android app (radi na telefonu)
- [x] Desktop GUI
- [x] Docker + Bioconda
- [x] Paul kontaktiran, Google Docs setup
- [x] Animesh kontaktiran, zainteresiran za suradnju
