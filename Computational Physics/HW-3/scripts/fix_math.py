import re

with open('README.md', 'r', encoding='utf-8') as f:
    text = f.read()

lines = text.split('\n')
new_lines = []
changed = False
for line in lines:
    if line.strip() == '$$' and line != '$$':
        new_lines.append('$$')
        changed = True
    else:
        new_lines.append(line)

new_text = '\n'.join(new_lines)
if text != new_text:
    with open('README.md', 'w', encoding='utf-8') as f:
        f.write(new_text)
    print("Fixed!")
else:
    print("No changes needed in just full line $$")

# now let's also remove any space after $$ for inline or mixed
new_text2 = re.sub(r'\$\$\s+', '$$\n', new_text)
new_text2 = re.sub(r'\s+\$\$', '\n$$', new_text2)

if new_text != new_text2:
    with open('README.md', 'w', encoding='utf-8') as f:
        f.write(new_text2)
    print("Fixed spaces around $$")
