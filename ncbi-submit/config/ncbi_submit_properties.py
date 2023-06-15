base_script_dir = "/root/"
templates_dir = base_script_dir + "templates"
schemas_dir = base_script_dir + "schemas/"
private_key_path = base_script_dir + ".ssh/id_rsa_ncbi"

# The URL to use when validating the submission XML.
validation_url = "https://www.ncbi.nlm.nih.gov/projects/biosample/validate/"

# ASCP command and default options
ascp_cmd = [
    base_script_dir + ".aspera/connect/bin/ascp",
    "-k2",
    "-L", "logs",
    "--file-manifest-path=logs",
    "--file-manifest=text"
]

schema_paths = {
    'submission': schemas_dir + 'submission.xsd',
    'bioproject': schemas_dir + 'bioproject.xsd',
    'biosample':  schemas_dir + 'biosample.xsd',
    'genome':     schemas_dir + 'genome.xsd'
}

ncbi_user = "asp-bioci"
ncbi_host = "upload.ncbi.nlm.nih.gov"
ncbi_sumbit_path = "submit/Production"
# ncbi_sumbit_path = "submit/Test"

bio_sample_reserved_attributes = {'sample_id',
                                  'bio_sample_package'}
library_reserved_attributes =    {'sample_id',
                                  'library_id',
                                  'data_type',
                                  'genome_coverage',
                                  'sequencing_technology',
                                  'assembly_method',
                                  'assembly_method_version',
                                  'assembly_name',
                                  'genome_representation',
                                  'expected_final_version',
                                  'publication_status',
                                  'publication_id',
                                  'publication_db_type',
                                  'publication_title',
                                  'publication_journal',
                                  'publication_year',
                                  'publication_volume',
                                  'publication_issue',
                                  'publication_pages_from',
                                  'publication_pages_to',
                                  'annotate_with_pgap',
                                  'notes',
                                  'gap_min_size',
                                  'gap_unknown_size',
                                  'gap_linkage_evidence'}
library_repeated_attributes =    {'authors'}
bio_sample_dup_attributes =      {'organism',
                                  'sample_title'}

compressed_content_types = {'application/gzip',
                            'application/x-bzip2',
                            'application/zip'}

# Categorized attributes are useful for optional XML element attributes. Search for
# "genome_representation_attrs" in genome-metadata.xml for an example.
library_categorized_attributes = {
    'genome_representation_description': 'genome_representation_attrs',
    'author_first_name': 'authors',
    'author_last_name': 'authors',
    'author_middle_name': 'authors',
    'author_name_suffix': 'authors'
}
