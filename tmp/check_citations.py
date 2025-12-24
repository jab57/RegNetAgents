import re

# Read the manuscript
with open('biorxiv/preprint_draft.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Extract all in-text citations
# Modified pattern to avoid matching standalone years (require citation context)
# Citations should be in format (number) or (number,number) or (number-number)
# Avoid matching (2021) which is likely a year
citation_pattern = r'\((\d+(?:[-–,;]\s*\d+)*)\)'
citations = []
positions = []

for match in re.finditer(citation_pattern, content):
    cite_text = match.group(1)

    # Skip if this looks like a year (single 4-digit number >= 1900)
    if cite_text.isdigit() and len(cite_text) == 4 and int(cite_text) >= 1900:
        continue

    # Skip if this is a very large number (likely not a citation)
    if cite_text.isdigit() and int(cite_text) > 100:
        continue

    # Parse individual numbers
    nums = []
    for part in re.split(r'[,;]\s*', cite_text):
        part = part.strip()
        if '-' in part or '–' in part:
            # Range like 32-35
            start, end = re.split(r'[-–]', part)
            nums.extend(range(int(start.strip()), int(end.strip())+1))
        else:
            nums.append(int(part))

    citations.extend(nums)
    positions.append((match.start(), nums))

print("=== CITATION AUDIT RESULTS ===\n")
print(f"Total citation instances: {len(citations)}")
print(f"Unique references cited: {len(set(citations))}")
print(f"Max reference number: {max(citations)}")
print(f"Min reference number: {min(citations)}")
print()

# Check for sequential order on first appearance
seen = set()
first_appearances = []
for pos, nums in positions:
    for num in nums:
        if num not in seen:
            first_appearances.append(num)
            seen.add(num)

print("=== FIRST APPEARANCE ORDER ===")
print(f"Expected: 1, 2, 3, 4, ..., {max(first_appearances)}")
print(f"Actual:   {', '.join(map(str, first_appearances[:20]))}...")
print()

# Check for gaps
all_refs = sorted(set(citations))
expected = list(range(1, max(all_refs) + 1))
missing = set(expected) - set(all_refs)

if missing:
    print("=== GAPS IN REFERENCE NUMBERING ===")
    print(f"Missing reference numbers: {sorted(missing)}")
    print()

# Check if first appearances are in order
out_of_order = []
for i in range(len(first_appearances) - 1):
    if first_appearances[i] > first_appearances[i+1]:
        out_of_order.append((first_appearances[i], first_appearances[i+1]))

if out_of_order:
    print("=== OUT OF ORDER FIRST APPEARANCES ===")
    for prev, curr in out_of_order[:10]:
        print(f"Reference {prev} appears before reference {curr}")
    if len(out_of_order) > 10:
        print(f"... and {len(out_of_order) - 10} more")
    print()

# Now read the References section to verify all cited refs exist
refs_section_match = re.search(r'## References\s*\n(.*?)(?=\n##|\Z)', content, re.DOTALL)
cited_but_not_in_refs = set()
in_refs_but_not_cited = set()
refs_gaps = set()

if refs_section_match:
    refs_text = refs_section_match.group(1)
    # Find all reference numbers in the References section
    ref_entries = re.findall(r'^\s*(\d+)\.\s+', refs_text, re.MULTILINE)
    ref_numbers = [int(n) for n in ref_entries]

    print(f"=== REFERENCES SECTION ===")
    print(f"Total references in References section: {len(ref_numbers)}")
    print(f"Reference numbers present: {min(ref_numbers)} to {max(ref_numbers)}")
    print()

    # Check if all cited refs exist in References
    cited_but_not_in_refs = set(citations) - set(ref_numbers)
    if cited_but_not_in_refs:
        print("=== CITATIONS WITHOUT REFERENCE ENTRIES ===")
        print(f"Cited in text but not in References section: {sorted(cited_but_not_in_refs)}")
        print()

    # Check if all refs are cited
    in_refs_but_not_cited = set(ref_numbers) - set(citations)
    if in_refs_but_not_cited:
        print("=== UNCITED REFERENCES ===")
        print(f"In References section but not cited in text: {sorted(in_refs_but_not_cited)}")
        print()

    # Check for gaps in References numbering
    refs_gaps = set(range(1, max(ref_numbers) + 1)) - set(ref_numbers)
    if refs_gaps:
        print("=== GAPS IN REFERENCES SECTION ===")
        print(f"Missing reference numbers: {sorted(refs_gaps)}")
        print()

print("=== SUMMARY ===")
if not missing and not out_of_order and not cited_but_not_in_refs and not in_refs_but_not_cited:
    print("SUCCESS: All citations are properly ordered and all references match!")
else:
    print("ISSUES FOUND:")
    if missing:
        print(f"  - {len(missing)} gaps in citation numbering")
    if out_of_order:
        print(f"  - {len(out_of_order)} out-of-order first appearances")
    if cited_but_not_in_refs:
        print(f"  - {len(cited_but_not_in_refs)} citations without reference entries")
    if in_refs_but_not_cited:
        print(f"  - {len(in_refs_but_not_cited)} uncited references")
