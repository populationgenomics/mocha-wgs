version 1.0

# ================================== #
# MoChA WGS - Preprocessing Workflow #
# ================================== #

# This workflow is designed to preprocess WGS data for MoChA analysis.
# The inputs to the workflow are a CRAM or BAM file and/or a pre-genotyped VCF.
#
# If supplying a VCF and use_existing_vcf is true,
# the VCF will be processed to be ready for MoChA analysis.
#
# If supplying a CRAM/BAM and run_bcftools is true,
# bcftools will be run to call and genotype variants, and the resulting VCF
# will be processed to be ready for MoChA analysis.
# Additionally, if a pre-genotyped VCF is also provided and
# restrict_bcftools_to_gatk_sites is true,
# the sites in the existing VCF will be used to restrict bcftools mpileup
#
# The VCF will first be processed through the mochatools bcftools plugin to
# add GC content annotations to the VCF.

workflow MochaWgsPreprocess {
    input {
        File? alignments
        File? alignments_index
        File? vcf
        File? vcf_index
        File ref_fasta
        File ref_fai
        File ref_dict
        String ref_name = "GRCh38"  # Currently only supports GRCh38 or GRCh37
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        File? intervals
        File? samples_file
        Int scatter_count = 10
        File? hapmap
        File? hapmap_index
        File? omni
        File? omni_index
        File? g1000
        File? g1000_index
        File? dbsnp
        File? dbsnp_index
        File? mills
        File? mills_index
        File? axiom_poly
        File? axiom_poly_index
        Int snps_max_gaussians = 6
        Int indels_max_gaussians = 4
        Boolean use_existing_vcf = true
        Boolean run_bcftools = true
        Boolean restrict_bcftools_to_gatk_sites = true
        Boolean filter_duplicate_reads = true
        Boolean filter_secondary_alignments = true
        Boolean filter_unmapped_reads = true

        # Runtime options
        String gatk_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/gatk:4.2.1.0"
        String mochatools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/mochatools:latest"
        String samtools_docker = "australia-southeast1-docker.pkg.dev/pb-dev-312200/somvar-images/samtools:latest"
        Int preemptible = 2
        Int max_retries = 2
        Int gatk_cpu = 4
        Int gatk_mem = 10
        Int gatk_mem_padding = 1
        Int bcftools_cpu = 4
        Int bcftools_mem = 10
        Int bcftools_mem_padding = 1
        Int samtools_cpu = 4
        Int samtools_mem = 10
        Int samtools_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
        Float cram_to_bam_multiplier = 6.0
    }

    if (use_existing_vcf) {
        File r_vcf = select_first([vcf])
        File r_vcf_index = select_first([vcf_index])

        call MochaAddGcContent as MochaAddGcContentOriginal {
            input:
                vcf = r_vcf,
                vcf_index = r_vcf_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                mochatools_docker = mochatools_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                bcftools_cpu = bcftools_cpu,
                bcftools_mem = bcftools_mem,
                bcftools_mem_padding = bcftools_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }

        call MochaFilterVcf as MochaFilterVcfOriginal {
            input:
                vcf = MochaAddGcContentOriginal.gc_vcf,
                vcf_index = MochaAddGcContentOriginal.gc_vcf_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                bcftools_docker = mochatools_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                bcftools_cpu = bcftools_cpu,
                bcftools_mem = bcftools_mem,
                bcftools_mem_padding = bcftools_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }
    }

    if (run_bcftools) {
        File r_alignments = select_first([alignments])
        File r_alignments_index = select_first([alignments_index])

        if (filter_duplicate_reads || filter_secondary_alignments || filter_unmapped_reads) {
            Boolean is_cram = basename(basename(r_alignments, ".bam"), ".cram") == basename(r_alignments, ".cram")
            Float multipler = if (is_cram) then cram_to_bam_multiplier else 2.5
            Int alignment_size = ceil(size(r_alignments, "GB") + size(r_alignments_index, "GB"))
            Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB"))
            Int samtools_disk = ceil((alignment_size * multipler) + ref_size + disk)
            call FilterBam {
                input:
                    alignments = r_alignments,
                    alignments_index = r_alignments_index,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    filter_duplicate_reads = filter_duplicate_reads,
                    filter_secondary_alignments = filter_secondary_alignments,
                    filter_unmapped_reads = filter_unmapped_reads,
                    samtools_docker = samtools_docker,
                    preemptible = preemptible,
                    max_retries = max_retries,
                    samtools_cpu = samtools_cpu,
                    samtools_mem = samtools_mem,
                    samtools_mem_padding = samtools_mem_padding,
                    samtools_disk = samtools_disk,
                    boot_disk_size = boot_disk_size
            }
        }

        scatter(chrom in chromosomes) {
            call BcftoolsMpileup {
                input:
                    vcf = vcf,
                    vcf_index = vcf_index,
                    alignments = select_first([FilterBam.filtered_bam, r_alignments]),
                    alignments_index = select_first([FilterBam.filtered_bam_index, r_alignments_index]),
                    samples = samples_file,
                    regions = chrom,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_name = ref_name,
                    restrict_bcftools_to_gatk_sites = restrict_bcftools_to_gatk_sites,
                    bcftools_docker = mochatools_docker,
                    preemptible = preemptible,
                    max_retries = max_retries,
                    bcftools_cpu = bcftools_cpu,
                    bcftools_mem = bcftools_mem,
                    bcftools_mem_padding = bcftools_mem_padding,
                    disk = disk,
                    boot_disk_size = boot_disk_size
            }
        }

        Array[File] mpileup_vcfs = select_all(BcftoolsMpileup.bcftools_vcf)
        Array[File] mpileup_vcf_indexes = select_all(BcftoolsMpileup.bcftools_vcf_index)
        
        call ConcatMpileupChrVCFs {
            input:
                vcfs = mpileup_vcfs,
                vcf_indexes = mpileup_vcf_indexes,
                bcftools_docker = mochatools_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                bcftools_cpu = bcftools_cpu,
                bcftools_mem = bcftools_mem,
                bcftools_mem_padding = bcftools_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }

        call MochaAddGcContent as MochaAddGcContentMpileup {
            input:
                vcf = ConcatMpileupChrVCFs.concat_vcf,
                vcf_index = ConcatMpileupChrVCFs.concat_vcf_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                mochatools_docker = mochatools_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                bcftools_cpu = bcftools_cpu,
                bcftools_mem = bcftools_mem,
                bcftools_mem_padding = bcftools_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }

        call MochaFilterVcf as MochaFilterVcfMpileup {
            input:
                vcf = MochaAddGcContentMpileup.gc_vcf,
                vcf_index = MochaAddGcContentMpileup.gc_vcf_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                bcftools_docker = mochatools_docker,
                preemptible = preemptible,
                max_retries = max_retries,
                bcftools_cpu = bcftools_cpu,
                bcftools_mem = bcftools_mem,
                bcftools_mem_padding = bcftools_mem_padding,
                disk = disk,
                boot_disk_size = boot_disk_size
        }
    }

    output {
        File? mocha_ready_gatk_vcf = MochaFilterVcfOriginal.mocha_filtered_vcf
        File? mocha_ready_gatk_vcf_index = MochaFilterVcfOriginal.mocha_filtered_vcf_index
        File? mocha_ready_mpileup_vcf = MochaFilterVcfMpileup.mocha_filtered_vcf
        File? mocha_ready_mpileup_vcf_index = MochaFilterVcfMpileup.mocha_filtered_vcf_index
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
        Int bcftools_cpu = 4
        Int bcftools_mem = 10
        Int bcftools_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (bcftools_mem - bcftools_mem_padding) * 1000
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
        cpu: bcftools_cpu
        memory: bcftools_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task FilterBam {
    input {
        File alignments
        File alignments_index
        File ref_fasta
        File ref_fai
        File ref_dict
        Boolean filter_duplicate_reads = true
        Boolean filter_secondary_alignments = true
        Boolean filter_unmapped_reads = true

        # Runtime options
        String samtools_docker
        Int preemptible = 2
        Int max_retries = 2
        Int samtools_cpu = 4
        Int samtools_mem = 10
        Int samtools_mem_padding = 1
        Int samtools_disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (samtools_mem - samtools_mem_padding) * 1000
    String sample_name = basename(basename(alignments, ".cram"), ".bam")
    Int exclude_flags = (if (filter_duplicate_reads) then 1024 else 0) + (if (filter_secondary_alignments) then 256 else 0) + (if (filter_unmapped_reads) then 4 else 0)
    String exclude_flags_param = (if (exclude_flags > 0) then "-F ~{exclude_flags}" else "")
    Int nthreads = if (samtools_cpu < 2) then 0 else samtools_cpu - 1

    command <<<
        samtools view \
            -h \
            -b \
            ~{exclude_flags_param} \
            -@ ~{nthreads} \
            -T ~{ref_fasta} \
            -o "~{sample_name}.filtered.bam" \
            ~{alignments}
        samtools index \
            -b \
            -@ ~{nthreads} \
            "~{sample_name}.filtered.bam" \
            "~{sample_name}.filtered.bam.bai"
    >>>

    output {
        File filtered_bam = "~{sample_name}.filtered.bam"
        File filtered_bam_index = "~{sample_name}.filtered.bam.bai"
    }

    runtime {
        docker: samtools_docker
        cpu: samtools_cpu
        memory: samtools_mem + " GB"
        disks: "local-disk " + samtools_disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task BcftoolsMpileup {
    input {
        File? vcf
        File? vcf_index
        File alignments
        File alignments_index
        File? samples
        String? regions
        File ref_fasta
        File ref_fai
        String ref_name = "GRCh38"  # Currently only supports GRCh38 or GRCh37
        Boolean restrict_bcftools_to_gatk_sites
        Int min_mapq = 10

        # Runtime options
        String bcftools_docker
        Int preemptible = 2
        Int max_retries = 2
        Int bcftools_cpu = 4
        Int bcftools_mem = 10
        Int bcftools_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (bcftools_mem - bcftools_mem_padding) * 1000
    String sample_name = basename(basename(alignments, ".cram"), ".bam")
    String ploidy = if (ref_name == "GRCh38" || ref_name == "GRCh37") then "--ploidy ~{ref_name}" else ""
    String regions_param = if (defined(regions)) then "-r ~{regions}" else ""
    Boolean use_vcf = (restrict_bcftools_to_gatk_sites && defined(vcf))
    String use_vcf_str = if (use_vcf) then "TRUE" else "FALSE"
    String regions_vcf = "regions.sites_only.vcf.gz"
    String regions_vcf_param = if (use_vcf) then ("-R " + regions_vcf) else ""

    command <<<
        if [ "~{use_vcf_str}" = "TRUE" ]
        then
            # Optional: generate a sites-only VCF from the original VCF
            bcftools view \
                -G \
                ~{regions_param} \
                ~{vcf} | \
            bcftools annotate \
                -x INFO \
                -Oz \
                -o ~{regions_vcf}
            tabix -s 1 -b 2 -e 2 ~{regions_vcf}
            # Only proceed if there are variant sites in the VCF
            NOTEMPTY="$(bcftools view -H ~{regions_vcf} | head -n 1 | wc -l | sed -E -e 's/^\s+//g')"
            if [ ! "$NOTEMPTY" -eq "1" ]
            then
                exit 0
            fi
        fi

        # Run bcftools mpileup and call to generate GT and AD fields
        bcftools mpileup \
            -d 8000 \
            -a "FORMAT/DP,FORMAT/AD" \
            -q ~{min_mapq} \
            -f ~{ref_fasta} \
            ~{regions_vcf_param} \
            ~{alignments} | \
        bcftools call \
            -mv \
            -f GQ \
            ~{ploidy} \
            ~{"--samples-file " + samples} \
            -Oz \
            -o ~{sample_name}.mpileup.unnorm.vcf.gz
        tabix -s 1 -b 2 -e 2 ~{sample_name}.mpileup.unnorm.vcf.gz
        bcftools norm \
            --fasta-ref ~{ref_fasta} \
            ~{sample_name}.mpileup.unnorm.vcf.gz \
            -Oz \
            -o ~{sample_name}.mpileup.vcf.gz
        tabix -s 1 -b 2 -e 2 ~{sample_name}.mpileup.vcf.gz
    >>>

    output {
        File? bcftools_vcf = "~{sample_name}.mpileup.vcf.gz"
        File? bcftools_vcf_index = "~{sample_name}.mpileup.vcf.gz.tbi"
    }

    runtime {
        docker: bcftools_docker
        cpu: bcftools_cpu
        memory: bcftools_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task ConcatMpileupChrVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes

        # Runtime options
        String bcftools_docker
        Int preemptible = 2
        Int max_retries = 2
        Int bcftools_cpu = 4
        Int bcftools_mem = 10
        Int bcftools_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (bcftools_mem - bcftools_mem_padding) * 1000
    String vcf_basename = basename(basename(basename(vcfs[0], ".gz"), ".vcf"), ".mpileup")

    command <<<
        bcftools concat \
            -a \
            -Oz \
            -o ~{vcf_basename}.mpileup.vcf.gz \
            ~{sep=" " vcfs}
        tabix -s 1 -b 2 -e 2 ~{vcf_basename}.mpileup.vcf.gz
    >>>

    output {
        File concat_vcf = "~{vcf_basename}.mpileup.vcf.gz"
        File concat_vcf_index = "~{vcf_basename}.mpileup.vcf.gz.tbi"
    }

    runtime {
        docker: bcftools_docker
        cpu: bcftools_cpu
        memory: bcftools_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}

task MochaFilterVcf {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fai

        # Runtime options
        String bcftools_docker
        Int preemptible = 2
        Int max_retries = 2
        Int bcftools_cpu = 4
        Int bcftools_mem = 10
        Int bcftools_mem_padding = 1
        Int disk = 100
        Int boot_disk_size = 12
    }

    Int command_mem = (bcftools_mem - bcftools_mem_padding) * 1000
    String vcf_basename = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        # Filter MoChA VCF
        bcftools view \
            --no-version \
            -h \
            ~{vcf} | \
        sed 's/^\(##FORMAT=<ID=AD,Number=\)\./\1R/' | \
        bcftools reheader \
            -h \
            /dev/stdin \
            ~{vcf} | \
        bcftools filter \
            --no-version \
            -Ou \
            -e "FMT/DP<10 | FMT/GQ<20" \
            --set-GT . | \
        bcftools annotate \
            --no-version \
            -Ou \
            -x ID,QUAL,^INFO/GC,^FMT/GT,^FMT/AD | \
        bcftools norm \
            --no-version \
            -Ou \
            -m -any \
            --keep-sum AD | \
        bcftools norm \
            --no-version \
            -Ob \
            -o ~{vcf_basename}.filtered.minimal.bcf \
            -f ~{ref_fasta} \
            --write-index
    >>>

    output {
        File mocha_filtered_vcf = "~{vcf_basename}.filtered.minimal.bcf"
        File mocha_filtered_vcf_index = "~{vcf_basename}.filtered.minimal.bcf.csi"
    }

    runtime {
        docker: bcftools_docker
        cpu: bcftools_cpu
        memory: bcftools_mem + " GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        bootDiskSizeGb: boot_disk_size
    }
}
