rule uORFreport:
    input:
        ufcplot="uORFs/xtail_cds_fc.pdf",
        urplot="uORFs/xtail_cds_r.pdf",
        cfcplot="uORFs/xtail_cds_fc.pdf",
        crplot="uORFs/xtail_cds_r.pdf"
    output:
        cfcplot=report("report/xtail_cds_fc.jpg", caption="../report/xtail_cds_fc.rst", category="CDS"),
./report.smk:        
        crplot=report("report/xtail_cds_r.jpg", caption="../report/xtail_cds_fc.rst", category="CDS"),
        ufcplot=report("report/xtail_uORFs_fc.jpg", caption="../report/xtail_uORFs_fc.rst", category="uORFs")
        urplot=report("report/xtail_uORFs_r.pdf", caption="../report/xtail_uORFs_fc.rst", category="uORFs") 
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p report; convert -density 150 -trim {input.ufcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.ufcplot}; convert -density 150 -trim {input.urplot}  -quality 100  -flatten -sharpen 0x1.0 {output.urplot}; convert -density 150 -trim {input.cfcplot}  -quality 100  -flatten -sharpen 0x1.0 {output.cfcplot}; convert -density 150 -trim {input.crplot}  -quality 100  -flatten -sharpen 0x1.0 {output.crplot}; ")

rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf" 
    output:
        report("report/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Ribotish")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    shell: ("mkdir -p report; convert -density 150 -trim {input}  -quality 100  -flatten -sharpen 0x1.0 {output};")
