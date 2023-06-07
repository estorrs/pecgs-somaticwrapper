$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 0
  prefix: --run-dir
  valueFrom: run
baseCommand:
- python
- /pecgs-somaticwrapper/sommaticwrapper/somaticwrapper.py
class: CommandLineTool
cwlVersion: v1.0
id: sommaticwrapper
inputs:
- id: sample
  inputBinding:
    position: '1'
  type: string
- id: tumor_bam
  inputBinding:
    position: '2'
  secondaryFiles:
  - .bai
  type: File
- id: normal_bam
  inputBinding:
    position: '3'
  secondaryFiles:
  - .bai
  type: File
- id: reference
  inputBinding:
    position: '0'
    prefix: --reference
  secondaryFiles:
  - .fai
  type: File
- id: rescue_genes
  inputBinding:
    position: '0'
    prefix: --rescue-genes
  type: File
- default: /miniconda/envs/somaticwrapper/bin:$PATH
  id: environ_PATH
  type: string?
label: sommaticwrapper
outputs:
- id: dnp_annotated_maf
  outputBinding:
    glob: run.dnp.annotated.maf
  type: File
- id: withmutect.maf
  outputBinding:
    glob: run.withmutect.maf
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/pecgs-somaticwrapper:0.0.1
- class: ResourceRequirement
  ramMin: 100000
- class: EnvVarRequirement
  envDef:
    PATH: $(inputs.environ_PATH)
