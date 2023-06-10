base_script_dir = "/home/cyverse/"
templates_dir = base_script_dir + "templates"
schemas_dir = base_script_dir + "schemas/"
private_key_path = base_script_dir + ".ssh/id_rsa_ncbi"

# ASCP command and default options
ascp_cmd = [
    base_script_dir + ".aspera/connect/bin/ascp",
    "-k2",
    "-L", "logs",
    "--file-manifest-path=logs",
    "--file-manifest=text"
]

submission_schema_path = schemas_dir + "submission.xsd"
bioproject_schema_path = schemas_dir + "bioproject.xsd"
biosample_schema_path = schemas_dir + "biosample.xsd"

ncbi_user = "asp-bioci"
ncbi_host = "upload.ncbi.nlm.nih.gov"
ncbi_sumbit_path = "submit/Production"
# ncbi_sumbit_path = "submit/Test"

bio_sample_reserved_attributes = {'sra_sample_id',
                                  'sra_bio_sample_package'}
library_reserved_attributes =    {'sra_sample_id',
                                  'library_id'}
bio_sample_dup_attributes =      {'organism',
                                  'sample_title'}

compressed_content_types = {'application/gzip',
                            'application/x-bzip2',
                            'application/zip'}
