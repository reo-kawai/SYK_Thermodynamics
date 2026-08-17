# AGENTS.md - SYK Thermodynamics Research Project

## Project Overview
This directory is dedicated to research on statistical mechanical applications of the SYK (Sachdev-Ye-Kitaev) model, including fluctuation theorems, thermalization, and the Mpemba effect.

## Communication Guidelines

### Language Rules
- **Japanese**: Use for all research conversations, reports, and content discussions
- **English**: Use exclusively for project management files (AGENTS.md, configuration files, code comments, and technical documentation)
- **Plan Mode**: Explain plans in Japanese for clarity and understanding

### Role Description
Assist with:
- Creating, running, and debugging numerical calculation code
- Organizing and visualizing computational results
- Writing and assisting with TeX manuscript preparation

## Directory Structure

Maintain strict separation of work by process:

```
SYK_Thermodynamics/
├── notes_md/              # Research plans and summaries (Markdown format)
├── notes_tex/             # Research summaries and generated PDFs (TeX format)
├── presentation/
│   ├── slide/             # Oral presentation materials (PDF, TeX, PPTX)
│   └── poster/            # Poster presentation materials (PDF, TeX, PPTX)
├── code/                  # Numerical calculation code
│   └── out/               # Computational output and results
├── references/            # Prior research (PDFs, TeX files)
└── AGENTS.md              # This file
```

### Directory Management Rules

1. **New Files**: Always check existing structure before creating new files. Place files in appropriate directories.

2. **Subdirectory Creation**: As files accumulate, create thematic subdirectories:
   - Organize by research topic/theme
   - Keep cross-cutting materials at parent level
   - Example: `notes_md/fluctuation_theorems/`, `notes_md/thermalization/`

3. **Naming Convention**:
   - Use descriptive names with hyphens for word separation
   - Include date prefix for time-sensitive materials: `YYYY-MM-DD_description.md`

## Code Management

### Code Organization (`/code`)
- Store all numerical computation code here
- Use Python, unless specified otherwise
- Include docstrings and comments in English
- Organize by research topic when file count grows

### Output Management (`/code/out`)
- Store raw computational results here
- Include metadata files describing computation parameters
- Organize by experiment/calculation date

### Code Standards
- Python: Use meaningful variable names, type hints where appropriate
- Include parameter documentation at file head
- Output files should include timestamp and parameter metadata

## File Naming Conventions

- **TeX/PDF**: descriptive names, e.g., `research_summary_fluctuation_theorems.tex`
- **Code**: descriptive names with underscores, e.g., `syk_dynamics_simulation.py`
- **Results**: date-prefixed where applicable, e.g., `2026-08-10_simulation_results.csv`
- **Notes**: date-prefixed for chronological tracking, e.g., `2026-08-10_notes_mpemba_effect.md`

## Workflow

1. **Research Planning**: Capture plans and notes in `/notes_md`
2. **Code Development**: Develop in `/code`, test locally
3. **Computation**: Run simulations, save results to `/code/out`
4. **Manuscript Preparation**: Draft TeX in `/notes_tex`, generate PDFs
5. **Presentation Preparation**: Organize slides/posters in `/presentation/`

## Questions and Clarification

Always ask for clarification before proceeding if:
- Task scope is unclear
- Multiple implementation approaches exist
- Existing patterns might be reused
- Directory placement is ambiguous

Maintain focus on code quality, clear documentation, and organized file management.

## LaTeX Writing in /notes_tex

When working on LaTeX files in `/notes_tex`, refer to `rules/latex-writing.md` for detailed guidelines on coding conventions, formatting, and physics content standards.
