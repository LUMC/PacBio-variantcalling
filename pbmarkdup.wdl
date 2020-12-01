version 1.0

task Pbmarkdup {
  input {
    File in_file
    String out_file
    Boolean? cross_library
    Boolean? rmdup
    Boolean? dup_file
    Boolean? clobber
    String? log_file
    String log_level = "INFO"

    Int num_threads = 8
    String memory = "20G"
    Int timeMinutes = 1440
    String dockerImage = "quay.io/biocontainers/pbmarkdup:1.0.2--0"
  }

  command {
    pbmarkdup \
      ~{in_file} \
      ~{out_file} \
      ~{true="--cross-library" false="" cross_library} \
      ~{true="--rmdup" false="" rmdup} \
      ~{true="--dup-file" false="" dup_file} \
      ~{true="--clobber" false="" clobber} \
      ~{true="--num-threads" false="" num_threads} \
      ~{true="--log-level" false="" log_level} \
      ~{true="--log-file" false="" log_file}
  }

  runtime {
        cpu: num_threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }


  parameter_meta {
    cross_library: "Identify duplicates across sequencing libraries (LB tag in read group)."
    rmdup: "Exclude duplicates from OUTFILE. Redundant when --dup-file is provided."
    dup_file: "Write duplicates to this file instead of OUTFILE."
    clobber: "Overwrite out_file if it exists."
    num_threads: "Number of threads to use, 0 means autodetection."
    log_level: "Set log level. Valid choices: TRACE, DEBUG, INFO, WARN, FATAL."
    log_file: "Log to a file, instead of stderr."
    in_file: "Input file as BAM, DATASET.XML, FASTA[.GZ], FASTQ[.GZ], or FOFN."
    out_file: "Output file as BAM, DATASET.XML, FASTA.GZ or FASTQ.GZ."
  }
}
