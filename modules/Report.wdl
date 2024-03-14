workflow Report {

  input {
    Array[File] ch_xenium_output
    Array[File] ch_baysor_segmentation
    Array[File] ch_nuclear_segmentation
    Array[File] ch_nuclear_segmentation_notebook
    Array[File] ch_baysor_config
  }
  
  call dumpParameters
  
  call diagnostics {
    input: 'notebook.ipynb' = ch_xenium_output[0]
  }
  
  call evaluation {
    input: 'notebook.ipynb' = ch_xenium_output[0]
  }
  
  call scanpy {
    input: 'notebook.ipynb' = ch_xenium_output[0]
  }
  
  call boundaries {
    input: 'notebook.ipynb' = ch_xenium_output[0]
  }
  
  call build {
    input:
      'parameter.md',
      'segmentation.py',
      'baysor.toml',
      'diagnostics.py',
      'evaluation.py',
      'scanpy.py',
      'boundaries.py'
  }

}

task dumpParameters {
  
  command <<<
    echo "# Parameters" > parameter.md
    echo '```json' >> parameter.md
    json_str=$(echo '${JsonOutput.toJson(params)}')
    echo "$json_str" | python -m json.tool >> parameter.md
    echo '```' >> parameter.md
  >>>
  
  output {
    File "parameter.md"
  }

}

task diagnostics {
  
  input {
    File 'notebook.ipynb'
    File 'data/xenium'
    File "data/transcripts.csv"
    File "data/transcripts_cellpose.csv"
  }
  
  command <<<
    jupyter nbconvert --to notebook --execute notebook.ipynb
  >>>
  
  output {
    File 'notebook.nbconvert.ipynb'
  }
  
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
    cpu: 8
    memory: "20 GB * task.attempt"
    maxRetries: 3
  }

}

task evaluation {
  
  input {
    File 'notebook.ipynb'
    File "data/transcripts.csv"
    File "data/transcripts_cellpose.csv"
  }
  
  command <<<
    jupyter nbconvert --to notebook --execute notebook.ipynb
  >>>
  
  output {
    File 'notebook.nbconvert.ipynb'
  }
  
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
    cpu: 8
    memory: "60 GB"
    maxRetries: 3
  }

}

task scanpy {
  
  input {
    File 'notebook.ipynb'
    File 'data/xenium'
    File "data/transcripts.csv"
  }
  
  command <<<
    export WIDTH=${params.report.width}
    export HEIGHT=${params.report.height}
    export X_OFFSET=${params.report.x_offset}
    export Y_OFFSET=${params.report.y_offset}
    jupyter nbconvert --to notebook --execute notebook.ipynb
  >>>
  
  output {
    File 'notebook.nbconvert.ipynb' as notebook
    File 'anndata.h5ad'
  }
  
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
    cpu: 8
    memory: "20 GB * task.attempt"
    maxRetries: 3
  }

}

task boundaries {
  
  input {
    File 'notebook.ipynb'
    File 'data/xenium'
    File "data/transcripts.csv"
    File "data/transcripts_cellpose.csv"
  }
  
  command <<<
    export WIDTH=${params.report.width}
    export HEIGHT=${params.report.height}
    export X_OFFSET=${params.report.x_offset}
    export Y_OFFSET=${params.report.y_offset}
    jupyter nbconvert --to notebook --execute notebook.ipynb
  >>>
  
  output {
    File 'notebook.nbconvert.ipynb'
  }
  
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
    cpu: 8
    memory: "20 GB * task.attempt"
    maxRetries: 3
  }

}

task build {

  input {
    File 'parameter.md'
    File 'segmentation.ipynb'
    File "baysor.toml"
    File 'diagnostics.ipynb'
    File 'evaluation.ipynb'
    File 'scanpy.ipynb'
    File 'boundaries.ipynb'
  }

  command <<<
    cp -r ${parameter.md} ${segmentation.ipynb} ${baysor.toml} ${diagnostics.ipynb} ${evaluation.ipynb} ${scanpy.ipynb} ${boundaries.ipynb} .
    echo "# Baysor config \n\n\`\`\`toml" > baysor_config.md
    cat baysor.toml >> baysor_config.md
    echo "\`\`\`" >> baysor_config.md
    jupyter-book build .
    mkdir report
    cp -r _build/html/* report/
  >>>
 
  output {
    Array[File] 'report'
  }
 
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
  }

}