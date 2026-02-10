from Bio import Phylo
from io import StringIO
from itertools import combinations
from Bio.Phylo.NewickIO import Writer
import re



# User-defined: specify the list of taxa to be analyzed and the outgroup taxon here.
taxa = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
outgroup = 'O'


# User-defined: specify input/output file paths here.

input_tree_file = "gamma_test_6.704_0.14916467780429596_M.tre"
quartet_file = "gamma_test_6.704_0.14916467780429596_M_quartet.tre"
adjust_file = "gamma_test_6.704_0.14916467780429596_M_quartet_adjust.tre"
filtered_file = "gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier.tre"
final_output_file = "gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier_modify.tre"
final_weight_file = "gamma_test_6.704_0.14916467780429596_M_quartet_adjust_deleteoutlier_modify_withoutweight_num.tre"

BRANCH_PRECISION = 10
MIN_BRANCH = 0.00001


########################################
# Write the tree (high precision)
########################################

def write_tree_with_precision(tree, handle, precision=10):
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            clade.branch_length = float(f"{clade.branch_length:.{precision}f}")

    writer = Writer([tree])
    writer.branch_length_formatter = lambda x: f"{x:.{precision}f}"
    writer.write(handle)



# Fix the root–O merge issue (if needed)
def adjust_tree_if_needed(tree_str):
    try:
        tree = Phylo.read(StringIO(tree_str.strip()), "newick")
    except Exception:
        return tree_str.strip()

    root = tree.root
    if len(root.clades) != 2:
        return tree_str.strip()

    c1, c2 = root.clades

    if len(c1.get_terminals()) == 1 and c1.get_terminals()[0].name == outgroup:
        o_clade, subtree = c1, c2
    elif len(c2.get_terminals()) == 1 and c2.get_terminals()[0].name == outgroup:
        o_clade, subtree = c2, c1
    else:
        return tree_str.strip()

    adjustment = subtree.branch_length
    if adjustment is None or adjustment <= 0:
        return tree_str.strip()

    o_clade.branch_length = (o_clade.branch_length or 0.0) + adjustment
    subtree.branch_length = 0.0

    out = StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue().strip().replace("\n", "")



# Step 1: Generate quartets
with open(input_tree_file) as f:
    tree_lines = [line.strip() for line in f if line.strip().endswith(";")]

with open(quartet_file, "w") as outfile:
    for triplet in combinations(taxa, 3):
        keep_taxa = set(triplet + (outgroup,))
        for line in tree_lines:
            try:
                tree = Phylo.read(StringIO(line), "newick")
                terminals = tree.get_terminals()

                if len(terminals) != 9:
                    continue

                for t in terminals:
                    if t.name not in keep_taxa:
                        tree.prune(t.name)

                handle = StringIO()
                write_tree_with_precision(tree, handle, BRANCH_PRECISION)
                outfile.write(handle.getvalue().strip() + "\n")
            except Exception as e:
                print("❌ Error:", e)


print("✅ Step 1 finish: quartet_file")



# Step 2: Quartet adjustment / cleanup

with open(quartet_file) as infile, open(adjust_file, "w") as outfile:
    for line in infile:
        line = line.strip()
        if not line:
            continue

        line = line.replace(":0.00000;", ";")

        line = line.replace(":0.00000", f":{MIN_BRANCH}")

        line = adjust_tree_if_needed(line)

        line = line.replace(":0.00000;", ";")
        line = line.replace(":0.00000", "")

        outfile.write(line + "\n")

print("✅ Step 2 finish: quartet_adjust")



# Step 3: Remove quartets with abnormal branch lengths (outliers)


def is_outlier(tree):
    total_length = sum(
        c.branch_length for c in tree.find_clades()
        if c.branch_length is not None
    )

    for c in tree.find_clades():
        if c.branch_length is None:
            continue

        if c.name == outgroup:
            if c.branch_length <= 0.0002:
                return True
        else:
            if c.branch_length <= 0.0002:
                return True
            if c.branch_length > 0.5 * total_length:
                return True

    return False


with open(adjust_file) as fin, open(filtered_file, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue

        tree = Phylo.read(StringIO(line), "newick")
        if not is_outlier(tree):
            fout.write(line + "\n")

print("✅ Step 3 finish: delete outlier")


# Step 4: Normalize Newick format (force O to be on the outside)


with open(filtered_file,'r') as f1, \
     open(final_output_file,'w') as f2:
    lines = f1.readlines()
    for line in lines:
        original_line=line
        for j in line.strip(';\n').split(','):
            if 'O' in j:
                f2.write('(')
                f2.write(j.strip(')'))
                f2.write(',')
                #f2.write(str(line[1:]))
                # ⭐ 这里再删 O
                rest = re.sub(r',O:[0-9\.eE+-]+', '', original_line)
                f2.write(str(rest[1:]))


print(f"✅ Step 4 finish:{final_output_file}")



# Step 5: Convert taxon labels to numbers (force O to be on the outside)


def normalize_tree_line(line: str) -> str:
    # 1. Remove all branch lengths of the form :number
    line = re.sub(r':[0-9.]+', '', line)

    # 2. Remove O (handle multiple possible positions)
    line = re.sub(r',O(?=[,)])', '', line)     # ,O
    line = re.sub(r'(?<=[(,])O,', '', line)    # O,
    line = re.sub(r'\(O\)', '', line)          # (O)

    # 3. Map A,B,C,... to 1,2,3,...
    line = re.sub(r'\bA\b', '1', line)
    line = re.sub(r'\bB\b', '2', line)
    line = re.sub(r'\bC\b', '3', line)
    line = re.sub(r'\bD\b', '4', line)
    line = re.sub(r'\bE\b', '5', line)
    line = re.sub(r'\bF\b', '6', line)
    line = re.sub(r'\bG\b', '7', line)
    line = re.sub(r'\bH\b', '8', line)

    # 4. Clean up potential extra punctuation
    line = re.sub(r',+', ',', line)
    line = re.sub(r'\(,', '(', line)
    line = re.sub(r',\)', ')', line)

    return line.strip()


def convert_tree_file(input_file: str, output_file: str):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            fout.write(normalize_tree_line(line) + '\n')

convert_tree_file(final_output_file, final_weight_file)
print(f"✅ Step 5 finish:{final_weight_file}")