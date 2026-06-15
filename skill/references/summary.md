This file is a merged representation of a subset of the codebase, containing specifically included files and files not matching ignore patterns, combined into a single document by Repomix.
The content has been processed where content has been compressed (code blocks are separated by ⋮---- delimiter).

# Summary

## Purpose

This is a reference codebase organized into multiple files for AI consumption.
It is designed to be easily searchable using grep and other text-based tools.

## File Structure

This skill contains the following reference files:

| File | Contents |
|------|----------|
| `project-structure.md` | Directory tree with line counts per file |
| `files.md` | All file contents (search with `## File: <path>`) |
| `tech-stack.md` | Languages, frameworks, and dependencies |
| `summary.md` | This file - purpose and format explanation |

## Usage Guidelines

- This file should be treated as read-only. Any changes should be made to the
  original repository files, not this packed version.
- When processing this file, use the file path to distinguish
  between different files in the repository.
- Be aware that this file may contain sensitive information. Handle it with
  the same level of security as you would the original repository.

## Notes

- Some files may have been excluded based on .gitignore rules and Repomix's configuration
- Binary files are not included in this packed representation. Please refer to the Repository Structure section for a complete list of file paths, including binary files
- Only files matching these patterns are included: core/**/*.py, camel_agents/*.py, camel_agents/tooling/*.py
- Files matching these patterns are excluded: **/__pycache__/**, **/test*/**, deps/**, scripts/**, output/**, examples/**
- Files matching patterns in .gitignore are excluded
- Files matching default ignore patterns are excluded
- Content has been compressed - code blocks are separated by ⋮---- delimiter
- Files are sorted by Git change count (files with more changes are at the bottom)

## Statistics

65 files | 19,273 lines

| Language | Files | Lines |
|----------|------:|------:|
| Python | 65 | 19,273 |

**Largest files:**
- `core/pathway/barrier_predictor.py` (1,114 lines)
- `camel_agents/tooling/neb.py` (1,084 lines)
- `camel_agents/workflow.py` (1,012 lines)
- `camel_agents/tooling/geometry.py` (808 lines)
- `core/kmc/catmap_predictor.py` (804 lines)
- `core/pathway/fairchem_predictor.py` (787 lines)
- `camel_agents/tooling/persistence.py` (783 lines)
- `camel_agents/gas_solid_digital_twin.py` (704 lines)
- `core/reconstruction/electrochemical_surface_predictor.py` (637 lines)
- `core/reconstruction/vssr_mc_predictor.py` (634 lines)