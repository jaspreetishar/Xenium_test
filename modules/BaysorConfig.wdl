task getMedianTranscriptsPerCell {
  input {
    File data
  }
  command <<<
    #!/usr/bin/env python
    
    import pandas as pd
    df = pd.read_csv("${data}")
    result = df[~(df['cell_id'] == "UNASSIGNED")].groupby("cell_id").size().median()
    print(round(result))
  >>>
  output {
    String stdout
  }
  runtime {
    docker: "maximilianheeg/docker-scanpy:v1.9.5"
  }
}

task Create {
  input {
    Int ch_min_mollecules_per_cell
    Int ch_min_molecules_per_segment
    String ch_scale
    String ch_scale_std
    Float ch_prior_segmentation_confidence
    Float ch_new_component_weight
    Int ch_n_clusters
  }
  command <<<
    cat > config.toml << EOF
    [data]
    x = "x_location"
    y = "y_location"
    z = "z_location"
    gene = "feature_name"
    min_molecules_per_cell = ${ch_min_mollecules_per_cell}
    min_molecules_per_segment = ${ch_min_molecules_per_segment}
    exclude_genes = "${params.baysor.exclude_genes}"

    [segmentation]
    scale = "${ch_scale}"
    scale_std = "${ch_scale_std}"
    n_clusters = ${ch_n_clusters}
    nuclei_genes = "${params.baysor.nuclei_genes}"
    cyto_genes   = "${params.baysor.cyto_genes}"
    prior_segmentation_confidence = ${ch_prior_segmentation_confidence}
    new_component_weight = ${ch_new_component_weight}
    EOF
  >>>
  output {
    File "config.toml"
  }
}

workflow estimateMinMoleculesPerCell {
  input {
    File ch_xenium_output
  }
  call getMedianTranscriptsPerCell { input: data = ch_xenium_output }
  scatter (tpc in getMedianTranscriptsPerCell.stdout) {
    rounded_tpc = Math.round(tpc.toInt() * params.baysor.min_molecules_per_cell_fraction.toFloat())
  }
  output {
    Int rounded_tpc
  }
}

workflow estimateMinMoleculesPerSegment {
  input {
    Int ch_min_molecules_per_cell
  }
  Int percent = Int(params.baysor.min_molecules_per_segment.replace("%", ""))
  scatter (mmps in ch_min_molecules_per_cell) {
    rounded_mmps = Math.round(mmps.toInt() / 100 * percent)
  }
  output {
    Int rounded_mmps
  }
}

workflow BaysorConfig {
  input {
    File ch_xenium_output
  }
  call estimateMinMoleculesPerCell as min_molecules_per_cell_estimation { input: ch_xenium_output = ch_xenium_output }
  call estimateMinMoleculesPerSegment as min_molecules_per_segment_estimation { input: ch_min_molecules_per_cell = min_molecules_per_cell_estimation.rounded_tpc }
  
  Int ch_min_molecules_per_cell = params.baysor.min_molecules_per_cell.toInt() < 0 ? min_molecules_per_cell_estimation.rounded_tpc : params.baysor.min_molecules_per_cell.toInt()
  
  Int ch_min_molecules_per_segment = params.baysor.min_molecules_per_segment.endsWith("%") ? min_molecules_per_segment_estimation.rounded_mmps : params.baysor.min_molecules_per_segment.toInt()
  
  String ch_scale = params.baysor.scale
  String ch_scale_std = params.baysor.scale_std
  Float ch_prior_segmentation_confidence = params.baysor.prior_segmentation_confidence.toFloat()
  Float ch_new_component_weight = params.baysor.new_component_weight.toFloat()
  Int ch_n_clusters = params.baysor.n_clusters.toInt()
  
  call Create {
    input:
      ch_min_mollecules_per_cell = ch_min_molecules_per_cell,
      ch_min_molecules_per_segment = ch_min_molecules_per_segment,
      ch_scale = ch_scale,
      ch_scale_std = ch_scale_std,
      ch_prior_segmentation_confidence = ch_prior_segmentation_confidence,
      ch_new_component_weight = ch_new_component_weight,
      ch_n_clusters = ch_n_clusters
  }
  output {
    File baysor_config = Create.config
  }
}
