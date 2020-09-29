# Upgrade Guide

## Upgrading from OpenGeode-GeosciencesIO v2.x.x to v3.0.0

### Motivations

This new major release is to formalize PyPI support for OpenGeode.

### Breaking Changes

- **CMake**: CMake minimum requirement has been upgraded to 3.14 to ease import of PyPI modules.


## Upgrading from OpenGeode-GeosciencesIO v1.x.x to v2.0.0

### Motivations

Homogenize CMake project and macro names with repository name.

### Breaking Changes

- Replace `OpenGeode_GeosciencesIO` by `OpenGeode-GeosciencesIO`. For example:
`OpenGeode_GeosciencesIO::geosciences` is replaced by `OpenGeode-GeosciencesIO::geosciences`.
