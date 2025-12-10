# BioRxiv Manuscript Conversion

## Quick Start

Convert all markdown files to Word + PDF:

```bash
cd biorxiv
python convert_manuscript.py
```

This will automatically generate:
- `regnetagents_preprint.docx` + `.pdf` (main manuscript)
- `supplementary_material.docx` + `.pdf` (supplementary material)

## What It Does

The `convert_manuscript.py` script:
1. Reads the markdown files (`preprint_draft.md`, `supplementary_material.md`)
2. Converts them to properly formatted Word documents (.docx)
3. Automatically generates PDFs from the Word documents
4. Handles page breaks (the `<div style="page-break-before: always"></div>` directives)
5. Formats tables, headings, and text properly

## Requirements

- Python 3.8+
- Libraries: `python-docx`, `docx2pdf`
- Microsoft Word must be installed (for PDF conversion on Windows)

## Troubleshooting

**PDF conversion fails**: Make sure Microsoft Word is installed and properly licensed.

**Manual alternative**: If PDF conversion fails, you can:
1. Run the script (it will still create .docx files)
2. Open the .docx files in Word
3. File → Save As → PDF

## Files Generated

All output files are created in the `biorxiv/` folder:
- `regnetagents_preprint.docx` - Main manuscript (Word)
- `regnetagents_preprint.pdf` - Main manuscript (PDF)
- `supplementary_material.docx` - Supplementary material (Word)
- `supplementary_material.pdf` - Supplementary material (PDF)

## Git Recovery

If you need the old conversion scripts, they're in git history:
```bash
git log --oneline -- biorxiv/convert_to_pdf.py
git checkout <commit-hash> -- biorxiv/convert_to_pdf.py
```
