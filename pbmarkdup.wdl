version 1.0

task Pbmarkdup {
  input {
    Array[File] in_file
    String out_file
    String log_file

    Int num_threads = 8
    String memory = "20G"
    Int timeMinutes = 1440
    String dockerImage = "quay.io/biocontainers/pbmarkdup:1.0.2--0"
  }

  command {
    pbmarkdup \
      ~{"--num-threads " + num_threads} \
      ~{sep=' ' in_file} \
      ~{out_file} > ~{log_file}
  }

  runtime {
        cpu: num_threads
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

  output {
      File outputBam = out_file
      File outputBamIndex = out_file + ".pbi"
      File logFile = log_file
    }

  parameter_meta {
    in_file: "Input file as BAM, DATASET.XML, FASTA[.GZ], FASTQ[.GZ], or FOFN."
    out_file: "Output file as BAM, DATASET.XML, FASTA.GZ or FASTQ.GZ."
  }
}
