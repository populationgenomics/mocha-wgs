version 1.0

# ================================== #
# MoChA WGS - Preprocessing Workflow #
# ================================== #

# This workflow is designed to preprocess WGS data for MoChA analysis.
# The inputs to the workflow are a CRAM or BAM file and a gVCF.
#
# The gVCF will be processed through GATK's GenotypeGVCFs to generate a VCF.
# The VCF will first be processed through the mochatools bcftools plugin to
# add GC content annotations to the VCF.

workflow MochaWgsPreprocess {
    input {
        File alignments
        File alignments_index
        File gvcf
        File gvcf_index
        File ref_fasta
        File ref_fai
        File ref_dict
        File intervals
        Int scatter_count = 10
        File hapmap
        File hapmap_index
        File omni
        File omni_index
        File g1000
        File g1000_index
        File dbsnp
        File dbsnp_index
        File mills
        File mills_index
        File axiom_poly
        File axiom_poly_index
        Int snps_max_gaussians = 6
        Int indels_max_gaussians = 4

        # Runtime options
        String gatk_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk:4.2.1.0"
        String mochatools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/mochatools:latest"
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    scatter (subintervals in SplitIntervals.interval_files) {
        call MochaGenotypeGVCFs {
            input:
                gvcf = gvcf,
                gvcf_index = gvcf_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                intervals = subintervals,
                dbsnp = dbsnp,
                dbsnp_index = dbsnp_index,
                gatk_docker = gatk_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                gatk_cpu = gatk_cpu,
                gatk_mem = gatk_mem,
                gatk_mem_padding = gatk_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }
    }

    String merged_vcf_basename = basename(MochaGenotypeGVCFs.gt_vcf[0], ".vcf.gz")

    call MergeVcfs as MergeGenotypedVcfs {
        input:
            vcfs = MochaGenotypeGVCFs.gt_vcf,
            vcf_indexes = MochaGenotypeGVCFs.gt_vcf_index,
            output_name = merged_vcf_basename,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaVariantFiltration {
        input:
            vcf = MergeGenotypedVcfs.merged_vcf,
            vcf_index = MergeGenotypedVcfs.merged_vcf_index,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaSitesOnlyVcf {
        input:
            vcf = MochaVariantFiltration.filtered_vcf,
            vcf_index = MochaVariantFiltration.filtered_vcf_index,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaSnpRecalibrator {
        input:
            vcf = MochaSitesOnlyVcf.so_vcf,
            vcf_index = MochaSitesOnlyVcf.so_vcf_index,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            hapmap = hapmap,
            hapmap_index = hapmap_index,
            omni = omni,
            omni_index = omni_index,
            g1000 = g1000,
            g1000_index = g1000_index,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            max_gaussians = snps_max_gaussians,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaIndelRecalibrator {
        input:
            vcf = MochaSitesOnlyVcf.so_vcf,
            vcf_index = MochaSitesOnlyVcf.so_vcf_index,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            mills = mills,
            mills_index = mills_index,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            axiom_poly = axiom_poly,
            axiom_poly_index = axiom_poly_index,
            max_gaussians = indels_max_gaussians,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaApplyVqsr {
        input:
            vcf = MochaVariantFiltration.filtered_vcf,
            vcf_index = MochaVariantFiltration.filtered_vcf_index,
            snp_recal = MochaSnpRecalibrator.snp_recal,
            snp_re1cal_index = MochaSnpRecalibrator.snp_recal_index,
            snp_tranches = MochaSnpRecalibrator.snp_tranches,
            indel_recal = MochaIndelRecalibrator.indel_recal,
            indel_recal_index = MochaIndelRecalibrator.indel_recal_index,
            indel_tranches = MochaIndelRecalibrator.indel_tranches,
            intervals = intervals,
            gatk_docker = gatk_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }

    call MochaAddGcContent {
        input:
            vcf = MochaApplyVqsr.vqsr_vcf,
            vcf_index = MochaApplyVqsr.vqsr_vcf_index,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            mochatools_docker = mochatools_docker,
            preemptible = preemptible,
            max_retries = max_retries,
            gatk_cpu = gatk_cpu,
            gatk_mem = gatk_mem,
            gatk_mem_padding = gatk_mem_padding,
            disk = disk,
            boot_disk_size = boot_disk_size
    }


    output {
        File genotyped_vcf = MergeGenotypedVcfs.merged_vcf
        File genotyped_vcf_index = MergeGenotypedVcfs.merged_vcf_index
        File filtered_vcf = MochaApplyVqsr.vqsr_vcf
        File filtered_vcf_index = MochaApplyVqsr.vqsr_vcf_index
        File mocha_ready_vcf = MochaAddGcContent.gc_vcf
        File mocha_ready_vcf_index = MochaAddGcContent.gc_vcf_index
    }
}

task SplitIntervals {
    input {
        File intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000

    command <<<
        mkdir interval-files
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" SplitIntervals \
            -R ~{ref_fasta} \
            -L ~{intervals} \
            -scatter ~{scatter_count} \
            -O interval-files
        cp interval-files/*.interval_list .
    >>>

    output {
        Array[File] interval_files = glob("*.interval_list")
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaGenotypeGVCFs {
    input {
        File gvcf
        File gvcf_index
        File ref_fasta
        File ref_fai
        File ref_dict
        File intervals
        File dbsnp
        File dbsnp_index

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(gvcf, ".gz"), ".g.vcf")

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" GenotypeGVCFs \
            -R ~{ref_fasta} \
            -V ~{gvcf} \
            -O ~{vcf_basename}.vcf.gz \
            -L ~{intervals} \
            -D ~{dbsnp} \
            -G StandardAnnotation \
            --only-output-calls-starting-in-intervals \
            --use-new-qual-calculator
    >>>

    output {
        File gt_vcf = "~{vcf_basename}.vcf.gz"
        File gt_vcf_index = "~{vcf_basename}.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MergeVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        String output_name

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String output_vcf = output_name + ".vcf.gz"

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" MergeVcfs \
            -I ~{sep=" -I " vcfs} \
            -O ~{output_vcf}
    >>>

    output {
        File merged_vcf = output_vcf
        File merged_vcf_index = output_vcf + ".tbi"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }

}

task MochaVariantFiltration {
    input {
        File vcf
        File vcf_index

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")
    Float excess_het_threshold = 54.69

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" VariantFiltration \
            -V ~{vcf} \
            -O ~{vcf_basename}.filtered.vcf.gz \
            --filter-expression "ExcessHet > ~{excess_het_threshold}" \
            --filter-name ExcessHet
    >>>

    output {
        File filtered_vcf = "~{vcf_basename}.filtered.vcf.gz"
        File filtered_vcf_index = "~{vcf_basename}.filtered.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaSitesOnlyVcf {
    input {
        File vcf
        File vcf_index

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" MakeSitesOnlyVcf \
            --INPUT ~{vcf} \
            --OUTPUT ~{vcf_basename}.sites_only.vcf.gz
    >>>

    output {
        File so_vcf = "~{vcf_basename}.sites_only.vcf.gz"
        File so_vcf_index = "~{vcf_basename}.sites_only.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaSnpRecalibrator {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fai
        File hapmap
        File hapmap_index
        File omni
        File omni_index
        File g1000
        File g1000_index
        File dbsnp
        File dbsnp_index
        Int max_gaussians = 6

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")
    String snp_recal_tranche_values = "100.0,99.95,99.9,99.8,99.6,99.5,99.4,99.3,99.0,98.0,97.0,90.0"
    String snp_recal_an_values = "QD,MQRankSum,ReadPosRankSum,FS,MQ,SOR,DP"

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" VariantRecalibrator \
            -V ~{vcf} \
            -O ~{vcf_basename}.snps.recal \
            --tranches-file ~{vcf_basename}.snps.tranches \
            --trust-all-polymorphic \
            -tranche $(echo ~{snp_recal_tranche_values} | sed -E -e 's/,/ -tranche /g') \
            -an $(echo ~{snp_recal_an_values} | sed -E -e 's/,/ -an /g') \
            -mode SNP \
            --max-gaussians ~{max_gaussians} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ~{hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 ~{omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 ~{g1000} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ~{dbsnp}
    >>>

    output {
        File snp_recal = "~{vcf_basename}.snps.recal"
        File snp_recal_index = "~{vcf_basename}.snps.recal.idx"
        File snp_tranches = "~{vcf_basename}.snps.tranches"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaIndelRecalibrator {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fai
        File mills
        File mills_index
        File dbsnp
        File dbsnp_index
        File axiom_poly
        File axiom_poly_index
        Int max_gaussians = 4

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")
    String indel_recal_tranche_values = "100.0,99.95,99.9,99.5,99.0,97.0,96.0,95.0,94.0,93.5,93.0,92.0,91.0,90.0"
    String indel_recal_an_values = "FS,ReadPosRankSum,MQRankSum,QD,SOR,DP"

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" VariantRecalibrator \
            -V ~{vcf} \
            -O ~{vcf_basename}.indels.recal \
            --tranches-file ~{vcf_basename}.indels.tranches \
            --trust-all-polymorphic \
            -tranche $(echo ~{indel_recal_tranche_values} | sed -E -e 's/,/ -tranche /g') \
            -an $(echo ~{indel_recal_an_values} | sed -E -e 's/,/ -an /g') \
            -mode INDEL \
            --max-gaussians ~{max_gaussians} \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 ~{mills} \
            -resource:axiomPoly,known=false,training=true,truth=false,prior=2.0 ~{axiom_poly} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~{dbsnp}
    >>>

    output {
        File indel_recal = "~{vcf_basename}.indels.recal"
        File indel_recal_index = "~{vcf_basename}.indels.recal.idx"
        File indel_tranches = "~{vcf_basename}.indels.tranches"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaApplyVqsr {
    input {
        File vcf
        File vcf_index
        File snp_recal
        File snp_re1cal_index
        File snp_tranches
        File indel_recal
        File indel_recal_index
        File indel_tranches
        File intervals

        # Runtime options
        String gatk_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")
    Float snp_vqsr_threshold = 99.7
    Float indel_vqsr_threshold = 99.7

    command <<<
        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" ApplyVQSR \
            -V ~{vcf} \
            -O tmp.snp.recalibrated.vcf \
            --truth-sensitivity-filter-level ~{snp_vqsr_threshold} \
            --tranches-file ~{snp_tranches} \
            --recal-file ~{snp_recal} \
            --create-output-variant-index true \
            --mode SNP

        gatk --java-options "-Xmx~{command_mem}m -Xms~{command_mem - 1000}m" ApplyVQSR \
            -V tmp.snp.recalibrated.vcf \
            -O ~{vcf_basename}.vqsr.vcf.gz \
            --truth-sensitivity-filter-level ~{indel_vqsr_threshold} \
            --tranches-file ~{indel_tranches} \
            --recal-file ~{indel_recal} \
            --create-output-variant-index true \
            --mode INDEL
    >>>

    output {
        File vqsr_vcf = "~{vcf_basename}.vqsr.vcf.gz"
        File vqsr_vcf_index = "~{vcf_basename}.vqsr.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaAddGcContent {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fai

        # Runtime options
        String mochatools_docker
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (gatk_mem - gatk_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        bcftools +mochatools \
            --no-version \
            -Oz \
            -o ~{vcf_basename}.gc.vcf.gz \
            ~{vcf} \
            -- \
            -t GC \
            -f ~{ref_fasta}

        tabix -s 1 -b 2 -e 2 ~{vcf_basename}.gc.vcf.gz
    >>>

    output {
        File gc_vcf = "~{vcf_basename}.gc.vcf.gz"
        File gc_vcf_index = "~{vcf_basename}.gc.vcf.gz.tbi"
    }

    runtime {
        docker: mochatools_docker
        cpu: gatk_cpu
        memory: gatk_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}
