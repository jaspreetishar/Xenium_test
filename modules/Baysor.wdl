version 1.0

import "../modules/BaysorConfig" as BaysorConfig_Module

workflow Baysor {
  input {
    File ch_xenium_output
    File ch_nuclear_segmentation
  }

  call BaysorConfig_Module.BaysorConfig as ch_baysor_config {
    input: 
      ch_xenium_output = ch_xenium_output
  }

  call tile_Xenium as ch_tile_xenium {
    input:
      transcripts = ch_nuclear_segmentation
  }

  call getNumberOfTranscripts as ch_get_number_of_transcripts {
    input:
      transcripts = ch_tile_xenium.out_files[0]
  }

  call runBaysor as ch_run_baysor {
    input:
      TRANSCRIPTS = ch_get_number_of_transcripts.TRANSCRIPTS
      transcripts = ch_tile_xenium.out_files[0]
      config = ch_baysor_config.config_file
  }
  call mergeTiles as ch_merge_tiles {
    input:
      transcripts = ch_run_baysor.segmentation
  }
  output {
    File transcripts = ch_merge_tiles.merged_transcripts
    File config = ch_baysor_config.config_file
  }
}


task tile_Xenium {
  input {
    File transcripts
  }
  
  command <<<
    tile-xenium ${transcripts} \
      --width ${tile.width} \
      --height ${tile.height} \
      --overlap ${tile.overlap} \
      --min-qv ${tile.qv} \
      --out-dir out/ \
      --minimal-transcripts ${tile.minimal_transcripts}
  >>>

  output {
    Array[File] out_files = "out_files"
  }
  
  runtime {
    docker: "maximilianheeg/tile-xenium:v0.1.2"
    cpu: 8
    memory: "10 GB * task.attempt"
    timeout: "2h * task.attempt"
    retryStrategy: {
      attempts: 3
    }
  }
}

task getNumberOfTranscripts {
  input {
    File transcripts
  }
  
  command <<<
    TRANSCRIPTS=$(cat ${transcripts} | wc -l)
    echo ${TRANSCRIPTS} > transcripts_count
  >>>
  
  output {
    Int TRANSCRIPTS
    File transcripts_count
  }

  runtime {
    docker: "ubuntu:latest"
    cpu: 1
    memory: "1 GB"
  }
}

task runBaysor {
  input {
    Int TRANSCRIPTS
    File transcripts
    File config
  }
  
  command <<<
    JULIA_NUM_THREADS=${task.cpus} baysor run \
      -c ${config} \
      -o out/ \
      -p \
      ${transcripts} \
      :cell_id
  >>>

  output {
    File segmentation
  }

  runtime {
    docker: "docker://maximilianheeg/baysor:v0.6.2"
    cpu: 8
    memory: "10 GB + (1 GB * round(${TRANSCRIPTS} / 1000000 * 20) * task.attempt)"
    timeout: "12h * task.attempt"
    maxRetries: 3
  }
}

task mergeTiles {
  input {
    Array[File] transcripts
  }
  
  command {
    merge-baysor ${transcripts} \
      --threshold ${merge.iou_threshold} \
      --additional-columns x \
      --additional-columns y \
      --additional-columns z \
      --additional-columns qv \
      --additional-columns overlaps_nucleus \
      --additional-columns gene \
      --outfile transcripts.csv
  }

  output {
    File transcripts.csv = "transcripts.csv"
  }

  runtime {
    docker: "docker://maximilianheeg/merge-baysor:v0.1.1"
    cpu: 8
    memory: "30 GB"
    
    # no equivalent of "time" in WDL, removing the following line; time: "2h * ${task.attempt}"
    # setting different resource requirements due to the removal of the "time" parameter, to indirectly influence how long a task is allowed to run.
    # in the nextlfow script, "memory" was set as "10 GB * ${task.attempt}". Here, setting it as "30 GB * ${task.attempt}" as a test. 
    # as no "maxRetries" is set, the default is applied, 1. because of this, the variable "${task.attempt}" has been removed.
  }
}

