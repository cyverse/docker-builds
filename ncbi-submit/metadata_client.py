__author__ = 'Paul Sarando'

import config.ncbi_submit_properties

import os
import json

class MetadataClient:
    """
    A class for transforming DE Data Store Metadata into NCBI XML template metadata.
    The get_metadata method expects the DE metadata file to contain JSON in the following format:
    {
        "id": "642ad94c-bd2a-11e4-891f-6abdce5a08d5",
        "label": "BioProject",
        "metadata": [
            {
                "attr": "BioProject-attribute1",
                "value": "BioProject-value1"
            },
            {
                "attr": "BioProject-attribute2",
                "value": "BioProject-value2"
            },
            ...
        ],
        "folders": [
            {
                "id": "6e80bd08-bd2a-11e4-891f-6abdce5a08d5",
                "label": "bio_sample_1",
                "metadata": [
                    {
                        "attr": "sample_id",
                        "value": "12345.biosample"
                    },
                    {
                        "attr": "bio_sample-reserved-attribute",
                        "value": "bio_sample-reserved-attribute"
                    },
                    {
                        "attr": "bio_sample-attribute1",
                        "value": "bio_sample-value1"
                    },
                    {
                        "attr": "bio_sample-attribute2",
                        "value": "bio_sample-value2"
                    },
                    ...
                ],
                "folders": [
                    {
                        "id": "7707eb7c-bd2a-11e4-891f-6abdce5a08d5",
                        "label": "library1",
                        "metadata": [
                            {
                                "attr": "library-attribute1",
                                "value": "library-value1"
                            },
                            {
                                "attr": "library-attribute2",
                                "value": "library-value2"
                            },
                            ...
                        ],
                        "files": [
                            {
                                "id": "7b6bc394-bd31-11e4-891f-6abdce5a08d5",
                                "label": "fasta-file1.gz",
                                "md5": "064848455cab44dade17bc5a3414d8b1"
                            },
                            {
                                "id": "cedd728e-bd31-11e4-891f-6abdce5a08d5",
                                "label": "fasta-file2.fgz",
                                "md5": "86913eb52f5f67aaaf711f4c32c2b0c6"
                            },
                            ...
                        ]
                    },
                    ...
                ]
            },
            ...
        ]
    }

    The get_metadata method will return NCBI XML template metadata in the following format:
    {
        "BioProject-attribute1": "BioProject-value1",
        "BioProject-attribute2": "BioProject-value2",
        "bio_samples": [
            {
                "name": "bio_sample_name",
                "sample_id": "12345.biosample",
                "bio_sample-reserved-attribute": "bio_sample-reserved-attribute",
                "attributes": [
                    {
                        "name": "bio_sample-attribute1",
                        "value": "bio_sample-value1"
                    },
                    {
                        "name": "bio_sample-attribute2",
                        "value": "bio_sample-value2"
                    }
                ]
            }
        ],
        "libraries": [
            {
                "name": "library_name",
                "sample_id": "12345.biosample",
                "library-reserved-attribute": "library-reserved-attribute",
                "attributes": [
                    {
                        "name": "library-attribute1",
                        "value": "library-value1"
                    },
                    {
                        "name": "library-attribute2",
                        "value": "library-value2"
                    }
                ],
                "files": [
                    {
                        "filename": "fasta-file1.gz",
                        "md5": "064848455cab44dade17bc5a3414d8b1"
                    },
                    {
                        "filename": "fasta-file2.fgz",
                        "md5": "86913eb52f5f67aaaf711f4c32c2b0c6"
                    }
                ]
            }
        ]
    }
    """

    def __init__(self, target_database='NONE'):
        self.target_database = target_database
        self.bio_sample_reserved_attributes = config.ncbi_submit_properties.bio_sample_reserved_attributes
        self.bio_sample_dup_attributes = config.ncbi_submit_properties.bio_sample_dup_attributes
        self.library_reserved_attributes = config.ncbi_submit_properties.library_reserved_attributes
        self.library_categorized_attributes = config.ncbi_submit_properties.library_categorized_attributes
        self.compressed_content_types = config.ncbi_submit_properties.compressed_content_types

    def get_metadata(self, json_file):
        with open(json_file) as file:
            metadata = json.load(file)

        if not metadata.get('metadata'):
            raise Exception("Could not load Bio Project metadata")
        if not metadata.get('folders'):
            raise Exception("Could not find Bio Project folder metadata")

        bio_project = {"object_id": metadata['id']}
        project_metadata = metadata['metadata']
        for attribute in project_metadata:
            if not (attribute.get('attr') and attribute.get('value')):
                continue

            attr = attribute['attr']
            value = attribute['value']
            if attr in bio_project:
                raise Exception("Duplicate '{0}' attribute found in Bio Project folder metadata.\nValues:\n{1}\n{2}".format(attr, bio_project[attr], value))

            bio_project[attr] = value

        bio_project['target_database'] = self.target_database
        bio_project['bio_samples'] = [self._parse_folder_metadata(folder)
                                      for folder in metadata['folders']]
        bio_project['libraries'] = [library
                                    for bio_sample in bio_project['bio_samples']
                                    for library in bio_sample['libraries']]

        return bio_project

    def _parse_folder_metadata(self, bio_sample_folder):
        if not bio_sample_folder.get('path'):
            raise Exception("Could not find Bio Sample folder path")

        bio_sample_name = os.path.basename(bio_sample_folder['path'])

        if not bio_sample_folder.get('id'):
            raise Exception("Could not find Bio Sample ID: {0}".format(bio_sample_name))
        if not bio_sample_folder.get('metadata'):
            raise Exception("Could not load Bio Sample metadata: {0}".format(bio_sample_name))
        if not bio_sample_folder.get('folders'):
            raise Exception("Could not find Bio Sample folder metadata: {0}".format(bio_sample_name))

        bio_sample = {"sample_id": bio_sample_folder['id'],
                      "name": bio_sample_name,
                      "attributes": []}
        metadata = bio_sample_folder['metadata']

        for attribute in metadata:
            if not (attribute.get('attr') and attribute.get('value')):
                continue

            attr = attribute['attr']
            value = attribute['value']
            if attr in self.bio_sample_reserved_attributes:
                if attr in bio_sample:
                    raise Exception("Duplicate '{0}' attribute found in Bio Sample metadata.\nValues:\n{1}\n{2}".format(attr, bio_sample[attr], value))

                bio_sample[attr] = value
            else:
                bio_sample['attributes'].append({'name': attr, 'value': value})
                if attr in self.bio_sample_dup_attributes:
                    # This attribute needs to be duplicated outside the attributes node as well.
                    bio_sample[attr] = value

        bio_sample['libraries'] = [self._parse_library_metadata(bio_sample['sample_id'], bio_sample['name'], library)
                                   for library in bio_sample_folder['folders']]

        return bio_sample

    def _parse_library_metadata(self, sample_id, bio_sample_name, library_folder):
        if not library_folder.get('path'):
            raise Exception("Could not find Bio Sample Library folder path")

        library_name = os.path.join(bio_sample_name, os.path.basename(library_folder['path']))

        if not library_folder.get('id'):
            raise Exception("Could not find Bio Sample Library ID: {0}".format(library_name))
        if not library_folder.get('metadata'):
            raise Exception("Could not load Bio Sample Library metadata: {0}".format(library_name))
        if not library_folder.get('files'):
            raise Exception("Could not find Bio Sample Library file metadata: {0}".format(library_name))

        library = {"library_id": library_folder['id'],
                   "name": library_name,
                   "sample_id": sample_id,
                   "attributes": []}
        metadata = library_folder['metadata']

        # Initialize categorized attributes.
        for category_name in set(self.library_categorized_attributes.values()):
            library[category_name] = {}

        for attribute in metadata:
            if not (attribute.get('attr') and attribute.get('value')):
                continue

            attr = attribute['attr']
            value = attribute['value']
            if attr in self.library_reserved_attributes:
                if attr in library:
                    raise Exception("Duplicate '{0}' attribute found in Bio Sample Library metadata.\nValues:\n{1}\n{2}".format(attr, library[attr], value))

                library[attr] = value
            elif attr in self.library_categorized_attributes:
                category_name = self.library_reserved_attributes[attr]
                library[category_name][attr] = value
            else:
                library['attributes'].append({'name': attr, 'value': value})

        library['files'] = [self._parse_file_metadata(library_name, file)
                            for file in library_folder['files']]

        return library

    def _parse_file_metadata(self, library_name, file_metadata):
        if not file_metadata.get('path'):
            raise Exception("Could not find Bio Sample Library file path: {0}".format(library_name))

        filename = os.path.basename(file_metadata['path'])

        if not file_metadata.get('content-type'):
            raise Exception("Could not find Bio Sample Library file content-type metadata: {0}".format(os.path.join(library_name, filename)))
        if file_metadata['content-type'] not in self.compressed_content_types:
            raise Exception("Bio Sample Library file does not appear to be compressed: {0}".format(os.path.join(library_name, filename)))

        if not file_metadata.get('md5'):
            raise Exception("Could not find Bio Sample Library file md5 metadata: {0}".format(os.path.join(library_name, filename)))

        return {"filename": filename,
                "md5":      file_metadata['md5']}

    def get_bio_project_file_paths(self, bio_project_metadata, bio_project_folder):
        return [os.path.join(bio_project_folder, library['name'], file['filename'])
                for library in bio_project_metadata['libraries']
                for file in library['files']]
