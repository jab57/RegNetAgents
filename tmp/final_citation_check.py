import re

# Read manuscript
with open('manuscript/preprint_draft.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Split at References section
refs_start = content.find('\n1. Babu MM')
if refs_start == -1:
    refs_start = len(content)

main_text = content[:refs_start]
refs_section = content[refs_start:]

# Extract citations from main text
citation_pattern = r'\((\d+(?:[-–]\d+)?(?:[,;]\s*\d+(?:[-–]\d+)?)*)\)'
all_cited = set()

for match in re.finditer(citation_pattern, main_text):
    cite_text = match.group(1)
    for part in re.split(r'[,;]\s*', cite_text):
        part = part.strip()
        if '-' in part or '–' in part:
            range_parts = re.split(r'[-–]', part)
            if len(range_parts) == 2:
                start_num = int(range_parts[0].strip())
                end_num = int(range_parts[1].strip())
                all_cited.update(range(start_num, end_num + 1))
        else:
            try:
                all_cited.add(int(part))
            except:
                pass

# Extract reference numbers from References section
ref_numbers = set()
for match in re.finditer(r'^\s*(\d+)\.\s+', refs_section, re.MULTILINE):
    ref_numbers.add(int(match.group(1)))

print(f"Citations in main text: {sorted(all_cited)}")
print(f"\nTotal cited: {len(all_cited)}")
print(f"Range: 1-{max(all_cited)}")
print(f"\nReferences in bibliography: {len(ref_numbers)}")
print(f"Range: 1-{max(ref_numbers)}")

# Check perfect match
if all_cited == ref_numbers == set(range(1, 43)):
    print("\n✓ PERFECT! All citations 1-42 match exactly.")
else:
    if all_cited != ref_numbers:
        print(f"\n✗ Mismatch between cited and references")
        print(f"  Cited but not in refs: {sorted(all_cited - ref_numbers)}")
        print(f"  In refs but not cited: {sorted(ref_numbers - all_cited)}")
