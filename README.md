# MELP
**Multimodal Evidence-based Locus Prioritization (MELP) Algorithm**

MELP is a method designed for efficient mining of key genes underlying target traits by leveraging population-scale multi‑omics data and enhancing gene signals through multi‑level detection.  

**The MELP pipeline consists of two main steps:**  

**Step 1: Signal detection across individual omics layers using conventional approaches**  
Examples include GWAS, TWAS, eQTL analysis, etc. The analysis codes used in this project are provided below:  
• GWAS (MLM model in TASSEL):  
  ``` sh 1_gwas.sh``` 
  
• TWAS (rrBLUP method):  
 ```  Rscript 2_twas.r``` 
  
• eQTL analysis (EMMAX software):  
  ``` sh 3_emmax.sh``` 
  
• Co‑expression analysis:  
 ```  Rscript 4_coexpress.r``` 
  

**Step 2: Integration of multi‑omics signals to obtain association scores for all genes**  
Run the integration script:  
``` python3 5_integration_features.py``` 
