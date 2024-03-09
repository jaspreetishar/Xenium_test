version 1.0

task getImageSize {
  input {
    File data
  }

  command <<<
    python <<EOF
    import tifffile

    tif = tifffile.TiffFile('data/morphology_mip.ome.tif')
    page = tif.pages[0]  # get shape and dtype of image in first page
    pixel = page.shape[0] * page.shape[1]
    bytes = pixel * 16
    print(bytes)
    EOF
  >>>

  output {
    Float bytes
  }

  runtime {
    docker: "docker://maximilianheeg/docker-cellpose:v2.2.3"
    cpu: 1
    memory: "1 GB"
  }
}

task cellpose {
  input {
    File nuclear_segmentation
    File data
    File 'models/DAPI'
    Float bytes
  }
  
  command {
    python nuclear_segmentation.py
  }

  output {
    File transcripts
    File notebook
  }

  runtime {
    docker: "docker://maximilianheeg/docker-cellpose:v2.2.3"
    cpu: 12
    memory: "10 GB + (1 GB * round(${bytes.toLong()} / 1000 / 1000 / 1000 * 7))"
    timeout: "8h"
  }
}

workflow NuclearSegmentation {
  input {
    File ch_xenium_output
    File ch_cellpose_model
  }
  call getImageSize as ch_size {
    input:
      data = ch_xenium_output
  }
  File segmentation_notebook = "nuclear_segmentation.ipynb"
  call cellpose as ch_process_nuclear_segmentation {
    input:
      nuclear_segmentation_ipynb = segmentation_notebook,
      data = ch_xenium_output,
      dapi_model = ch_cellpose_model,
      bytes = ch_size.bytes
  }
  output {
    File transcripts = ch_process_nuclear_segmentation.transcripts
    File notebook = ch_process_nuclear_segmentation.notebook
  }
}
