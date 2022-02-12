# GEMmaker: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.0 - [date]

This release is a major reconfiguration of GEMmaker to meet the nf-core standards [nf-core](https://nf-co.re/). It also includes multiple bug fixes

### `Added`

- Support for nf-core compatibility and functions. With nf-core compatibility come an easier configuration and execution step. Therefore configuration is much different from earlier versions of GEMmaker.
- Added a new failed runs report.  Runs that fail to download or dump from the SRA format to fastq format no longer cause the pipeline to fail, but are skipped and information about the failure is stored in a new HTML run failure report.
- Fully updated the documentation to be compatible with the 2.1.0 version.

### `Fixed`

- Fixed bugs in downloading SRA files using Aspera. Previously any SRA files that required the `fasp` protocol would fail.
- Added additional checks to make sure input files are present and to prevent GEMmaker from running if all dependencies are not met.

### `Dependencies`

- No changes to dependencies were made for this release.

### `Deprecated`

- The end to end (e2e) version of the workflow is no longer available. This version of GEMmaker was provided in as separate workflow file. It was not compatible with nf-core requirements and was no longer needed given recent improvements to the workflow.

## v1.1 - 2021-01-10

Provides a variety of bug fixes and updated tools.

## v1.0 - 2019-02-06

The first stable relase.

## v0.9 - 2018-06-05

This is a release just prior to v1.0 to allow for additional testing for minor bug fixes and other additions that may be needed for a stable v1.0 release.
