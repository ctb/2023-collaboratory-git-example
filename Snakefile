rule download_genome:
    shell:
        "wget https://osf.io/8sm92/download -O ecoli-rel606.fa.gz"

rule map_reads:
    shell: """
        minimap2 -ax sr ecoli-rel606.fa.gz SRR2584857_1.fastq.gz > SRR2584857_1.x.ecoli-rel606.sam
    """

rule sam_to_bam:
    shell: """
        samtools view -b -F 4 SRR2584857_1.x.ecoli-rel606.sam > SRR2584857_1.x.ecoli-rel606.bam
     """

rule sort_bam:
    shell: """
        samtools sort SRR2584857_1.x.ecoli-rel606.bam > SRR2584857_1.x.ecoli-rel606.bam.sorted
    """

rule gunzip_fa:
    shell: """  
        gunzip -c ecoli-rel606.fa.gz > ecoli-rel606.fa
    """

rule call_variants:
    shell: """
        bcftools mpileup -Ou -f ecoli-rel606.fa SRR2584857_1.x.ecoli-rel606.bam.sorted > SRR2584857_1.x.ecoli-rel606.pileup
        bcftools call -mv -Ob SRR2584857_1.x.ecoli-rel606.pileup -o SRR2584857_1.x.ecoli-rel606.bcf
        bcftools view SRR2584857_1.x.ecoli-rel606.bcf > SRR2584857_1.x.ecoli-rel606.vcf
    """

rule filter:
    shell:"""
        bcftools filter -Ov -e 'QUAL<40 || DP<10 || GT!="1/1"' SRR2584857_1.x.ecoli-rel606.sensitive.vcf > SRR2584857_1.x.ecoli-rel606.specific.vcf
    """  
    
