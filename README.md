# SIRT1–GSDME cell-death switch in hepatocellular carcinoma (HCC)

Boolean-network modeling (GINsim), Monte Carlo seed consensus, in silico
epistasis/rescue, and independent TCGA-LIHC expression validation for a
MALAT1 / miR-204-5p / SIRT1 / p53 / CASP3 / GSDME logical model of
pyroptotic competence in liver cancer.

> 📄 Manuscript: *Boolean modeling and transcriptomic integration of the miR-204-5p/MALAT1/SIRT1 axis in the regulation of GSDME-dependent pyroptosis in hepatocellular carcinoma)*.

---

## Repository structure

```
GSDME-pyroptosis-HCC/
├── model/
│   └── GINsim-miR_204_GSDME_Pyroptosis.zginml   # Boolean network (GINsim)
├── scripts/
│   ├── bulkrna-seq-gsdme.R          
│   ├── tcga-lihc-gsdme.R              
├── results/                
├── LICENSE
└── README.md
```


---

## License

Code is released under the [MIT License](LICENSE). The GINsim model file
(`model/*.zginml`) is provided for reproducibility of the companion
manuscript;
