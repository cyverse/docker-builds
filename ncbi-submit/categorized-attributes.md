# Categorized Attributes

Categorized attributes provide a way to include optional attributes in an XML element that is associated with a
library. This feature has not been implemented for higher-level attributes yet. To use them, first add one or more
attributes to a category in your configuration file. For example:

``` python
library_categorized_attributes = {
    'foo': 'technical_terms,
    'bar': 'technical_terms',
    'baz': 'technical_terms'
}
```

Once you have the categorized attributes defined, You can use `py:attrs` to instruct Genshi to add the optional
attributes.

``` xml
<tech-terms py:attrs="library.technical_terms" />
```

With this configuration, if your library metadata contains any of the categorized attributes then they will be included
in the generated XML. Suppose, for example, that the metadata looks like this:

``` json
{
    "id": "642ad94c-bd2a-11e4-891f-6abdce5a08d5",
    "label": "BioProject",
    "metadata": [
        {
            "attr": "BioProject-attribute1",
            "value": "BioProject-value1"
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
                            "attr": "foo",
                            "value": "quux"
                        },
                        {
                            "attr": "bar",
                            "value": "blargh"
                        },
                        {
                            "attr": "baz",
                            "value": "blrfl"
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
```

When the template is rendered, the resulting XML will look like this:

``` xml
<tech-terms foo="quux" bar="blargh" baz="blrfl" \>
```
