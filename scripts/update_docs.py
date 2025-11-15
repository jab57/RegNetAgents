#!/usr/bin/env python3
"""
Intelligent Documentation Update Tool for RegNetAgents

This script scans documentation in docs/ and biorxiv/ folders to:
1. Auto-fix mechanical issues (old names, dates, broken links)
2. Detect content accuracy issues (outdated sections, broken code examples)
3. Generate review checklists for human judgment

Usage:
    python scripts/update_docs.py --check              # Report only
    python scripts/update_docs.py --fix                # Fix mechanical issues
    python scripts/update_docs.py --fix --dates-only   # Update only dates
    python scripts/update_docs.py --review-checklist   # Generate review TODO list
"""

import os
import re
import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Set
import subprocess

# Configuration
PROJECT_ROOT = Path(__file__).parent.parent
DOCS_FOLDERS = [PROJECT_ROOT / "docs", PROJECT_ROOT / "biorxiv", PROJECT_ROOT]
CODE_FOLDERS = [PROJECT_ROOT / "regnetagents", PROJECT_ROOT / "scripts", PROJECT_ROOT / "tests"]

# Patterns to detect and fix
OLD_NAMES = {
    "GREmLN": "RegNetAgents",
    "gremln": "regnetagents",
    "GREMLN": "REGNETAGENTS",
    "GREmLN_": "RegNetAgents_",
    "gremln_": "regnetagents_",
}

# Files that were renamed
RENAMED_FILES = {
    "gremln_langgraph_mcp_server.py": "regnetagents_langgraph_mcp_server.py",
    "gremln_langgraph_workflow.py": "regnetagents_langgraph_workflow.py",
    "GREmLN_Analysis_Pipeline.md": "REGNETAGENTS_Analysis_Pipeline.md",
    "GREmLN_MCP_SETUP.md": "REGNETAGENTS_MCP_SETUP.md",
    "GREmLN_CONFERENCE_POSTER.md": "REGNETAGENTS_CONFERENCE_POSTER.md",
}


class DocumentIssue:
    """Represents an issue found in documentation"""

    def __init__(self, file_path: Path, issue_type: str, line_num: int,
                 description: str, auto_fixable: bool = False, suggestion: str = None):
        self.file_path = file_path
        self.issue_type = issue_type
        self.line_num = line_num
        self.description = description
        self.auto_fixable = auto_fixable
        self.suggestion = suggestion

    def __repr__(self):
        fix_status = "[OK] Auto-fixable" if self.auto_fixable else "[!] Needs review"
        location = f"{self.file_path.name}:{self.line_num}" if self.line_num else self.file_path.name
        return f"[{fix_status}] {location} - {self.description}"


class DocUpdateTool:
    """Main tool for analyzing and updating documentation"""

    def __init__(self, dry_run=True):
        self.dry_run = dry_run
        self.issues: List[DocumentIssue] = []
        self.stats = {
            "files_scanned": 0,
            "issues_found": 0,
            "auto_fixed": 0,
            "needs_review": 0
        }

    def scan_all_docs(self) -> List[DocumentIssue]:
        """Scan all documentation folders"""
        print("[*] Scanning documentation...\n")

        for folder in DOCS_FOLDERS:
            if folder.exists():
                self._scan_folder(folder)

        self.stats["issues_found"] = len(self.issues)
        return self.issues

    def _scan_folder(self, folder: Path):
        """Scan a single folder recursively"""
        # Special case: if scanning project root, only scan README.md
        if folder == PROJECT_ROOT:
            readme = folder / "README.md"
            if readme.exists():
                self.stats["files_scanned"] += 1
                self._scan_file(readme)
        else:
            # For other folders, scan recursively
            for file_path in folder.rglob("*.md"):
                self.stats["files_scanned"] += 1
                self._scan_file(file_path)

    def _scan_file(self, file_path: Path):
        """Scan a single markdown file"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
                lines = content.split('\n')

            # Check for various issues
            self._check_old_names(file_path, lines)
            self._check_outdated_dates(file_path, lines)
            self._check_broken_links(file_path, lines)
            self._check_renamed_files(file_path, lines)
            self._check_code_examples(file_path, lines)
            self._check_staleness(file_path)

        except Exception as e:
            self.issues.append(DocumentIssue(
                file_path, "error", 0, f"Error reading file: {e}", False
            ))

    def _check_old_names(self, file_path: Path, lines: List[str]):
        """Check for old project names"""
        for line_num, line in enumerate(lines, 1):
            for old_name, new_name in OLD_NAMES.items():
                # Skip if it's in a code block or already correct
                if old_name in line and new_name not in line:
                    # Check if it's not part of a longer word
                    pattern = r'\b' + re.escape(old_name) + r'\b'
                    if re.search(pattern, line):
                        self.issues.append(DocumentIssue(
                            file_path, "old_name", line_num,
                            f"Found old name '{old_name}' -> should be '{new_name}'",
                            auto_fixable=True,
                            suggestion=new_name
                        ))

    def _check_outdated_dates(self, file_path: Path, lines: List[str]):
        """Check for outdated 'Last Updated' dates"""
        current_year = datetime.now().year
        current_date = datetime.now().strftime("%Y-%m-%d")

        for line_num, line in enumerate(lines, 1):
            # Look for "Last Updated:" or similar patterns
            if re.search(r'Last\s+Updated:', line, re.IGNORECASE):
                # Extract date
                date_match = re.search(r'(\d{4})-(\d{2})-(\d{2})', line)
                if date_match:
                    date_str = date_match.group(0)
                    try:
                        doc_date = datetime.strptime(date_str, "%Y-%m-%d")
                        days_old = (datetime.now() - doc_date).days

                        if days_old > 30:  # More than a month old
                            self.issues.append(DocumentIssue(
                                file_path, "outdated_date", line_num,
                                f"Last updated {days_old} days ago ({date_str})",
                                auto_fixable=True,
                                suggestion=current_date
                            ))
                    except ValueError:
                        pass

    def _check_broken_links(self, file_path: Path, lines: List[str]):
        """Check for broken internal file references"""
        for line_num, line in enumerate(lines, 1):
            # Find markdown links and file references
            # [text](file.md) or `file.py` or just file.md
            links = re.findall(r'\[([^\]]+)\]\(([^)]+)\)', line)
            backticks = re.findall(r'`([^`]+\.(py|md|txt|json))`', line)

            for text, link in links:
                if not link.startswith('http') and not link.startswith('#'):
                    # Internal file reference
                    self._verify_file_exists(file_path, line_num, link)

            for match in backticks:
                file_ref = match[0]
                if '/' in file_ref or '\\' in file_ref:
                    self._verify_file_exists(file_path, line_num, file_ref)

    def _verify_file_exists(self, doc_path: Path, line_num: int, file_ref: str):
        """Verify a referenced file exists"""
        # Skip false positives
        false_positives = [
            'mailto:',  # Email links
            '%APPDATA%',  # Windows environment variables
            '~/',  # Unix home directory (not a real path to check)
            'python ',  # Command examples
            'pip ',  # Command examples
        ]

        if any(fp in file_ref for fp in false_positives):
            return  # Skip these

        # Try to resolve relative to doc location and project root
        # Also check common folders for bare filenames
        possible_paths = [
            doc_path.parent / file_ref,
            PROJECT_ROOT / file_ref,
            PROJECT_ROOT / file_ref.lstrip('./'),
            PROJECT_ROOT / "scripts" / file_ref,  # Check scripts folder
            PROJECT_ROOT / "tests" / file_ref,    # Check tests folder
        ]

        exists = any(p.exists() for p in possible_paths)

        if not exists:
            self.issues.append(DocumentIssue(
                doc_path, "broken_link", line_num,
                f"Referenced file not found: {file_ref}",
                auto_fixable=False
            ))

    def _check_renamed_files(self, file_path: Path, lines: List[str]):
        """Check for references to renamed files"""
        for line_num, line in enumerate(lines, 1):
            for old_file, new_file in RENAMED_FILES.items():
                if old_file in line:
                    self.issues.append(DocumentIssue(
                        file_path, "renamed_file", line_num,
                        f"References old filename '{old_file}' -> should be '{new_file}'",
                        auto_fixable=True,
                        suggestion=new_file
                    ))

    def _check_code_examples(self, file_path: Path, lines: List[str]):
        """Check code examples in documentation"""
        in_code_block = False
        code_lang = None
        code_start = 0
        code_lines = []

        for line_num, line in enumerate(lines, 1):
            # Detect code block start
            if line.strip().startswith('```'):
                if not in_code_block:
                    in_code_block = True
                    code_lang = line.strip()[3:].lower()
                    code_start = line_num + 1
                    code_lines = []
                else:
                    # Code block end - analyze it
                    if code_lang in ['python', 'py', 'bash', 'sh']:
                        self._analyze_code_block(file_path, code_start, code_lang, code_lines)
                    in_code_block = False
                    code_lang = None
            elif in_code_block:
                code_lines.append(line)

    def _analyze_code_block(self, file_path: Path, start_line: int, lang: str, code: List[str]):
        """Analyze a code block for potential issues"""
        code_text = '\n'.join(code)

        # Check for imports/references to renamed modules
        if lang in ['python', 'py']:
            for old_name, new_name in OLD_NAMES.items():
                if f'import {old_name}' in code_text or f'from {old_name}' in code_text:
                    self.issues.append(DocumentIssue(
                        file_path, "code_outdated", start_line,
                        f"Code example imports old module name '{old_name}'",
                        auto_fixable=False
                    ))

        # Check for references to files that don't exist
        if lang in ['bash', 'sh']:
            for line in code:
                # Look for python file.py patterns
                py_files = re.findall(r'python\s+([^\s]+\.py)', line)
                for py_file in py_files:
                    self._verify_file_exists(file_path, start_line, py_file)

    def _check_staleness(self, file_path: Path):
        """Check if documentation is stale compared to code changes"""
        try:
            # Get last modification time of the doc
            doc_mtime = os.path.getmtime(file_path)
            doc_date = datetime.fromtimestamp(doc_mtime)

            # Get git modification times for related code files
            related_code_files = self._find_related_code_files(file_path)

            for code_file in related_code_files:
                if code_file.exists():
                    code_mtime = os.path.getmtime(code_file)
                    code_date = datetime.fromtimestamp(code_mtime)

                    # If code changed more than 7 days after doc
                    if (code_date - doc_date).days > 7:
                        self.issues.append(DocumentIssue(
                            file_path, "stale_content", 0,
                            f"Code file {code_file.name} modified {code_date.strftime('%Y-%m-%d')} but doc not updated since {doc_date.strftime('%Y-%m-%d')}",
                            auto_fixable=False
                        ))
                        break  # Only report once per doc
        except Exception:
            pass  # Skip if git commands fail

    def _find_related_code_files(self, doc_path: Path) -> List[Path]:
        """Find code files that might be related to this doc"""
        doc_name = doc_path.stem.lower()
        related = []

        # Look for similarly named files
        for code_folder in CODE_FOLDERS:
            if code_folder.exists():
                for code_file in code_folder.rglob("*.py"):
                    if doc_name in code_file.stem.lower():
                        related.append(code_file)

        return related

    def fix_mechanical_issues(self, dates_only=False):
        """Auto-fix mechanical issues"""
        if self.dry_run:
            print("[*] DRY RUN: Would fix the following issues:\n")
        else:
            print("[*] Fixing mechanical issues...\n")

        # Group issues by file
        issues_by_file: Dict[Path, List[DocumentIssue]] = {}
        for issue in self.issues:
            if issue.auto_fixable:
                if dates_only and issue.issue_type != "outdated_date":
                    continue
                issues_by_file.setdefault(issue.file_path, []).append(issue)

        # Fix each file
        for file_path, file_issues in issues_by_file.items():
            self._fix_file(file_path, file_issues)

    def _fix_file(self, file_path: Path, issues: List[DocumentIssue]):
        """Fix all issues in a single file"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            fixed_count = 0

            for issue in sorted(issues, key=lambda x: x.line_num, reverse=True):
                if issue.line_num == 0:
                    continue

                line_idx = issue.line_num - 1
                original_line = lines[line_idx]

                if issue.issue_type == "old_name":
                    # Replace old name with new name
                    for old_name, new_name in OLD_NAMES.items():
                        pattern = r'\b' + re.escape(old_name) + r'\b'
                        lines[line_idx] = re.sub(pattern, new_name, lines[line_idx])

                elif issue.issue_type == "outdated_date":
                    # Update date
                    lines[line_idx] = re.sub(
                        r'\d{4}-\d{2}-\d{2}',
                        issue.suggestion,
                        lines[line_idx]
                    )

                elif issue.issue_type == "renamed_file":
                    # Replace old filename with new filename
                    for old_file, new_file in RENAMED_FILES.items():
                        lines[line_idx] = lines[line_idx].replace(old_file, new_file)

                if lines[line_idx] != original_line:
                    fixed_count += 1
                    print(f"  [OK] Fixed {file_path.name}:{issue.line_num} - {issue.issue_type}")

            if not self.dry_run and fixed_count > 0:
                with open(file_path, 'w', encoding='utf-8') as f:
                    f.writelines(lines)
                self.stats["auto_fixed"] += fixed_count

        except Exception as e:
            print(f"  [ERROR] Error fixing {file_path.name}: {e}")

    def generate_report(self):
        """Generate a summary report"""
        print("\n" + "="*70)
        print("DOCUMENTATION UPDATE REPORT")
        print("="*70 + "\n")

        print(f"Files scanned: {self.stats['files_scanned']}")
        print(f"Issues found: {self.stats['issues_found']}\n")

        if not self.issues:
            print("[OK] All documentation is up to date!\n")
            return

        # Group by issue type
        by_type: Dict[str, List[DocumentIssue]] = {}
        for issue in self.issues:
            by_type.setdefault(issue.issue_type, []).append(issue)

        # Count auto-fixable vs needs review
        auto_fixable = sum(1 for i in self.issues if i.auto_fixable)
        needs_review = len(self.issues) - auto_fixable

        print(f"Auto-fixable: {auto_fixable}")
        print(f"Needs review: {needs_review}\n")

        # Show issues by type
        for issue_type, issues in sorted(by_type.items()):
            print(f"\n>> {issue_type.upper().replace('_', ' ')} ({len(issues)} issues)")
            print("-" * 70)

            # Group by file
            by_file: Dict[Path, List[DocumentIssue]] = {}
            for issue in issues:
                by_file.setdefault(issue.file_path, []).append(issue)

            for file_path, file_issues in sorted(by_file.items()):
                print(f"\n  {file_path.relative_to(PROJECT_ROOT)}")
                for issue in sorted(file_issues, key=lambda x: x.line_num):
                    prefix = "    [OK]" if issue.auto_fixable else "    [!]"
                    location = f":{issue.line_num}" if issue.line_num else ""
                    print(f"{prefix}{location} {issue.description}")

        print("\n" + "="*70)
        print("\nRECOMMENDATIONS:")
        print("-" * 70)

        if auto_fixable > 0:
            print(f"\n1. Fix {auto_fixable} mechanical issues automatically:")
            print("   python scripts/update_docs.py --fix")

        if needs_review > 0:
            print(f"\n2. Review {needs_review} content issues manually:")
            print("   python scripts/update_docs.py --review-checklist")

        print()

    def generate_review_checklist(self, output_file: Path = None):
        """Generate a markdown checklist for manual review"""
        if output_file is None:
            output_file = PROJECT_ROOT / "docs" / "DOCUMENTATION_REVIEW_CHECKLIST.md"

        needs_review = [i for i in self.issues if not i.auto_fixable]

        if not needs_review:
            print("[OK] No items need manual review!\n")
            return

        # Generate markdown
        content = ["# Documentation Review Checklist", ""]
        content.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        content.append(f"Total items: {len(needs_review)}\n")
        content.append("---\n")

        # Group by file
        by_file: Dict[Path, List[DocumentIssue]] = {}
        for issue in needs_review:
            by_file.setdefault(issue.file_path, []).append(issue)

        for file_path, file_issues in sorted(by_file.items()):
            rel_path = file_path.relative_to(PROJECT_ROOT)
            content.append(f"## {rel_path}\n")

            for issue in sorted(file_issues, key=lambda x: x.line_num):
                location = f":{issue.line_num}" if issue.line_num else ""
                content.append(f"- [ ] **{issue.issue_type}**{location} - {issue.description}")

            content.append("")

        content.append("\n---\n")
        content.append("## Instructions\n")
        content.append("1. Review each item above")
        content.append("2. Update the relevant documentation")
        content.append("3. Check off items as you complete them")
        content.append("4. Delete this file when all items are addressed")

        checklist_text = '\n'.join(content)

        if not self.dry_run:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(checklist_text)
            print(f"[OK] Review checklist saved to: {output_file.relative_to(PROJECT_ROOT)}\n")
        else:
            print("REVIEW CHECKLIST (preview):")
            print("="*70)
            print(checklist_text[:1000])
            if len(checklist_text) > 1000:
                print("\n... (truncated)\n")


def main():
    parser = argparse.ArgumentParser(
        description="Intelligent documentation update tool for RegNetAgents",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/update_docs.py --check
  python scripts/update_docs.py --fix
  python scripts/update_docs.py --fix --dates-only
  python scripts/update_docs.py --review-checklist
        """
    )

    parser.add_argument('--check', action='store_true',
                       help='Check for issues without fixing (default mode)')
    parser.add_argument('--fix', action='store_true',
                       help='Automatically fix mechanical issues')
    parser.add_argument('--dates-only', action='store_true',
                       help='Only update dates (use with --fix)')
    parser.add_argument('--review-checklist', action='store_true',
                       help='Generate markdown checklist for manual review')

    args = parser.parse_args()

    # Default to check mode if nothing specified
    if not any([args.check, args.fix, args.review_checklist]):
        args.check = True

    # Create tool instance
    dry_run = not args.fix
    tool = DocUpdateTool(dry_run=dry_run)

    # Scan all docs
    tool.scan_all_docs()

    # Execute requested action
    if args.fix:
        tool.fix_mechanical_issues(dates_only=args.dates_only)

    if args.review_checklist:
        tool.generate_review_checklist()

    # Always show report
    tool.generate_report()


if __name__ == "__main__":
    main()
