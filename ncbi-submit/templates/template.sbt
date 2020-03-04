{% with library = metadata.libraries[0] %}\
Submit-block ::= {
  contact {
    contact {
      name name {
        last "${metadata.contact_last_name}",
        first "${metadata.contact_first_name}",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "${metadata.organization_address_institution}",
        div "${metadata.organization_address_department}",
        city "${metadata.organization_address_city}",
        sub "${metadata.get('organization_address_state', '')}",
        country "${metadata.organization_address_country}",
        street "${metadata.organization_address_street}",
        email "${metadata.contact_email}",
        postal-code "${metadata.organization_address_postal_code}"
      }
    }
  },
  cit {
    authors {
      names std {
        {
{% for author in library.authors %}
          name name {
            last "${author.author_last_name}",
            first "${author.author_first_name}",
            middle "${author.get('author_middle_name', '')}",
            initials "",
            suffix "${author.get('author_name_suffix', '')}",
            title ""
          }
{% end %}
        }
      },
      affil std {
        affil "${metadata.organization_address_institution}",
        div "${metadata.organization_address_department}",
        city "${metadata.organization_address_city}",
        sub "${metadata.get('organization_address_state', '')}",
        country "${metadata.organization_address_country}",
        street "${metadata.organization_address_street}",
        postal-code "${metadata.organization_address_postal_code}"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
{% if library.publication_status == 'eUnpublished' %}\
    gen {
      cit "unpublished",
      authors {
        names std {
          {
{% for author in library.authors %}
            name name {
              last "${author.author_last_name}",
              first "${author.author_first_name}",
              middle "${author.get('author_middle_name', '')}",
              initials "",
              suffix "${author.get('author_name_suffix', '')}",
              title ""
            }
{% end %}
          }
        }
      },
      title "${library.get('publication_title', '')}"
    }
{% end %}\
{% if library.publication_status != 'eUnpublished' %}\
    article {
      title {
        name "${library.get('publication_title', '')}"
      },
      authors {
        names std {
          {
{for author in library.authors %}
            name name {
              last "${author.author_last_name}",
              first "${author.author_first_name}",
              middle "${author.get('author_middle_name', '')}",
              initials "",
              suffix "${author.get('author_name_suffix', '')}",
              title ""
            }
{% end %}
          }
        }
      },
      from journal {
        title {
          jta "${library.get('publication_journal', '')}"
        },
        imp {
          date std {
            year ${library.get('publication_year', '')}"
          },
          volume "${library.get('publication_volume', '')}",
          issue "${library.get('publication_issue', '')}",
          pages "${library.get('publication_pages_from', '')}-${library.get('publication_pages_to', '')}"{% if library.publication_status == 'eInPress' %},
          prepub in-press
{% end %}\
        }
      }
    }
{% end %}\
  }
}
{% end %}\
