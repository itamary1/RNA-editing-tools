nextflow.enable.dsl=2


process  disjoint_bed_file{
  input:
    path bed6file_regions
  output:
    // path 'splited_files_dir'
    path "disjointed_files/*bed"
  script:
    out_disjointed_files_dir = params.disjounted_bed_dir ? params.disjounted_bed_dir : "disjointed_files"
    """
    bed6file_regions_no_suff=$bed6file_regions
    bed6file_regions_no_suff=\${bed6file_regions_no_suff%%.bed*}
    mkdir disjointed_files
    $params.which_disjoint_script_command -i $bed6file_regions -o ${out_disjointed_files_dir}/\${bed6file_regions_no_suff}.%s.bed3.bed --strand
    """

}