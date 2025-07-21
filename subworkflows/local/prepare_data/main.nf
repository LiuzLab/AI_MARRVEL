include {
    RESTRUCTURE_LOCAL_DATA_EXCEPT_VEP;
    RESTRUCTURE_LOCAL_DATA_ONLY_VEP;
    STORE_S3_BUCKET_DATA_EXCEPT_VEP;
    STORE_S3_BUCKET_DATA_ONLY_VEP;
    SPLIT_DATA;
    BUILD_REFERENCE_INDEX;
} from "../../../modules/local/prepare_data"

include {
    GENERATE_MANIFEST_JSON
} from "../../../modules/local/singleton"

workflow PREPARE_DATA {
    if (params.ref_dir) {
        GENERATE_MANIFEST_JSON(file(params.ref_dir))
        data_except_vep = RESTRUCTURE_LOCAL_DATA_EXCEPT_VEP(file(params.ref_dir))
        data_only_vep = RESTRUCTURE_LOCAL_DATA_ONLY_VEP(file(params.ref_dir))
    } else {
        error "Not Implemented"
        // data_except_vep = STORE_S3_BUCKET_DATA_EXCEPT_VEP()
        // data_only_vep = STORE_S3_BUCKET_DATA_ONLY_VEP()
    }

    (
        chrmap_file,
        ensembl_to_location_file,
        ensembl_to_symbol_file,
        ref_filter_bed_file,

        ref_annot_dir,
        ref_var_tier_dir,
        ref_merge_expand_dir,
        ref_mod5_diffusion_dir,
        ref_predict_new_dir,
        ref_model_inputs_dir,

        phrank_tuple,
        omim_tuple,
        gnomad_tuple,
        vep_tuple
    ) = SPLIT_DATA(
        data_except_vep,
        data_only_vep,
        params.bed_filter ? file(params.bed_filter) : moduleDir.resolve("assets/NO_FILE"),
        params.exome_filter,
    )

    fasta_tuple = BUILD_REFERENCE_INDEX()

    emit:
    data = chrmap_file.map{["chrmap_file", it]}.concat(
        ensembl_to_location_file.map{["ensembl_to_location_file", it]},
        ensembl_to_symbol_file.map{["ensembl_to_symbol_file", it]},
        ref_filter_bed_file.map{["ref_filter_bed_file", it]},

        ref_annot_dir.map{["ref_annot_dir", it]},
        ref_var_tier_dir.map{["ref_var_tier_dir", it]},
        ref_merge_expand_dir.map{["ref_merge_expand_dir", it]},
        ref_mod5_diffusion_dir.map{["ref_mod5_diffusion_dir", it]},
        ref_predict_new_dir.map{["ref_predict_new_dir", it]},
        ref_model_inputs_dir.map{["ref_model_inputs_dir", it]},

        fasta_tuple.map{["fasta_tuple", it]},
        gnomad_tuple.map{["gnomad_tuple", it]},
        omim_tuple.map{["omim_tuple", it]},
        phrank_tuple.map{["phrank_tuple", it]},
        vep_tuple.map{["vep_tuple", it]}
    ).toList().map{ it.collectEntries { k, v -> [(k): v]} }
}