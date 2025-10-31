import sys
import os

def list_files(startpath):
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compare_outputs.py <comparison_input_directory>")
        sys.exit(1)
    comparison_input_dir = sys.argv[1]
    print(f"Comparing outputs in directory: {comparison_input_dir}")
    list_files(comparison_input_dir)
