// WDL pipeline for 10x segmentation

import "./modules/Logo.wdl" as Logo
import "./modules/Baysor.wdl" as Baysor
import "./modules/NuclearSegmentation.wdl" as NuclearSegmentation
import "./modules/Report.wdl" as Report

workflow {

  input {
  
      File xenium_path
      File cellpose_model = "$baseDir/models/DAPI"

  }

  call Logo.Logo {}

  call NuclearSegmentation.NuclearSegmentation {
    input:
    xenium_output = xenium_path,
    cellpose_model = cellpose_model
  }

  call Baysor.Baysor {
    input:
    xenium_output = xenium_path,
    nuclear_segmentation = NuclearSegmentation.NuclearSegmentation.transcripts
  }

  call Report.Report {
    input:
    xenium_output = xenium_path,
    baysor_segmentation = Baysor.Baysor.transcripts,
    nuclear_segmentation = NuclearSegmentation.NuclearSegmentation.transcripts,
    nuclear_segmentation_notebook = NuclearSegmentation.NuclearSegmentation.notebook,
    baysor_config = Baysor.Baysor.config
  }
}