rule ribotishreport:
    input:
        "ribotish/{condition}-{replicate}-qual.pdf" 
    output:
        report("report/{condition}-{replicate}-qual.jpg", caption="../report/ribotishquality.rst", category="Ribotish")
    conda:
        "../envs/imagemagick.yaml"
    threads: 1
    log: "logs/{condition}-{replicate}_ribotishreport.log"
    shell: ("mkdir -p report; convert -density 150 -trim {input}  -quality 100  -flatten -sharpen 0x1.0 {output};")
