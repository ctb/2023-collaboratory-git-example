rule all:
    input:
        "SRR2584857_1.x.ecoli-rel606.specific.vcf",

rule download_genome:
    output: "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O {output}"

rule map_reads:
    input: ref="ecoli-rel606.fa.gz", reads="SRR2584857_1.fastq.gz"
    output: "SRR2584857_1.x.ecoli-rel606.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    input: "SRR2584857_1.x.ecoli-rel606.sam"
    output: "SRR2584857_1.x.ecoli-rel606.bam"
    shell: """
        samtools view -b -F 4 {input} > {output}
     """

rule sort_bam:
    input: "SRR2584857_1.x.ecoli-rel606.bam"
    output: "SRR2584857_1.x.ecoli-rel606.bam.sorted"
    shell: """
        samtools sort {input} > {output}
    """

rule gunzip_fa:
    input:
        "ecoli-rel606.fa.gz"
    output:
        "ecoli-rel606.fa"
    shell: """  
        gunzip -c {input} > {output}
    """

rule call_variants:
    input:
        ref="ecoli-rel606.fa",
        bamsort="SRR2584857_1.x.ecoli-rel606.bam.sorted"
    output:
        pileup="SRR2584857_1.x.ecoli-rel606.pileup",
        bcf="SRR2584857_1.x.ecoli-rel606.bcf",
        vcf="SRR2584857_1.x.ecoli-rel606.sensitive.vcf"
    shell: """
        bcftools mpileup -Ou -f {input.ref} {input.bamsort} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """

rule filter:
    input: "SRR2584857_1.x.ecoli-rel606.sensitive.vcf"
    output:"SRR2584857_1.x.ecoli-rel606.specific.vcf"
    shell:"""
        bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' {input} > {output}
    """  
    
