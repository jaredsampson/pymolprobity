# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!-- Possible subheadings:
### Added
### Changed
### Deprecated
### Fixed
### Removed
### Security
-->


## Unreleased

### Added
- Config file for use with `bumpversion` utility.

### Changed
- Format of this changelog file.

### Fixed
- Improved instructions for incentive build installation and PATH setup.



## v0.1.4

### Added
- Citation instructions and release DOIs via GitHub/Zenodo.
- This changelog file.



## v0.1.3

### Added
- Python3 compatibility by Thomas Holder ([@speleo3][]).
- Instructions for symlinking to Phenix-installed Reduce and Probe.

### Changed
- Installation instructions to use Github master zip file by default.

### Fixed
- Several bugs related to Tk interface, by Thomas Holder ([@speleo3][]).

### Removed
- `make` targets "updateenv" and "run".



## v0.1.2

### Fixed
- Parsing of Flipkin output for some atoms ("invalid literal for float" error).



## v0.1.1

### Changed
- Redirect some low-level output to logger instead of console.

### Fixed
- Provide a more informative error message when required executables are not found in PATH.
- An error that occurred when trying to load a non-existant object.



## v0.1.0

- Initial release


<!-- Contributor GitHub Profile Links -->
  [@speleo3]: https://github.com/speleo3
