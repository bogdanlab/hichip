H3k27ac-HiChIP in LNCaP
===================================================

The repository contains the intermediate files reported in
*H3k27ac-HiChIP in prostate cell lines identifies risk genes for prostate cancer susceptibility*
Claudia Giambartolomei* and Ji-Heui Seo*, Tommer Schwarz, Malika Kumar Freund, Ruth Dolly Johnson, Sandor Spisak, Sylvan C. Baca, Alexander Gusev, Nicholas Mancuso, Bogdan Pasaniuc and Matthew L. Freedman.
bioRxiv 2020.

Results
-------

`LNCaP_back0_loops_longrange.txt.gz` HiChIP loops in longrange format viewable in https://epigenomegateway.readthedocs.io/en/latest/tracks.html#hic. 

`paintor_1causals.txt.gz` contains the results from fine-mapping using Paintor (1 causal) across 130 regions. 
<!--- contains the results from fine-mapping using Paintor (1 causal) across 137 regions. In the analyses we removed regions that did not reach P < 5*10^-8, and we replaced the fine-mapped results in two regions in the 8q24 (127.6–129.0Mb, "chr8.127849659.128199916.rs16901979,rs183373024,rs10086908,rs12543663", "chr8.128260673.128560038.rs1447295,rs620861,rs6983267") with JAM fine-mapping (Supplementary Data 3 of Matejcic et al.).  --->

`Thib_eGenes_bonf_bestSNP` contains the cis-eQTLs association statistics in Thibodeau.

`TCGA_eGenes_bonf_bestSNP` contains the cis-eQTLs association statistics in TCGA.

`PrCa_GeneList_Used.csv` contains the list of genes somatically mutated in prostate cancer were extracted from three publications. 

`paintor_genes_small_pad5000.csv` contains the overlaps of HiChIP looing data with TSS and fine-mapped SNPs. 


Software and support
--------------------
If you have any questions or comments please contact claudia.giambartolomei@gmail.com
