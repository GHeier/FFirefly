"""
Scans README.md files in each category/calculation/method directory,
extracts the Quick Description, and updates User.md with links.
"""

from pathlib import Path
import re


def get_project_root():
    """Get the project root directory."""
    script_dir = Path(__file__).parent  # src/config/write
    return script_dir.parent.parent.parent  # FFirefly


def extract_quick_description(readme_path):
    """
    Extract the Quick Description section from a README.md file.

    Args:
        readme_path: Path to the README.md file

    Returns:
        The quick description text, or None if not found
    """
    try:
        with open(readme_path, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        return None

    # Find "## Quick Description" section and extract text until next "##"
    pattern = r'## Quick Description\s*\n(.*?)(?=\n##|\Z)'
    match = re.search(pattern, content, re.DOTALL)

    if match:
        description = match.group(1).strip()
        # Skip template placeholder text
        if description.startswith("One sentence description") or \
           description.startswith("One or two sentence description"):
            return None
        return description

    return None


def get_relative_readme_path(readme_path, project_root):
    """Get the relative path from project root to README."""
    return readme_path.relative_to(project_root)


def scan_category_readmes(categories, project_root):
    """
    Scan all category/calculation/method directories for README.md files.

    Args:
        categories: Dict from categories.py
        project_root: Path to project root

    Returns:
        Dict of {category: {calculation: {method: (readme_path, description)}}}
    """
    src_dir = project_root / "src"
    results = {}

    for category, calculations in categories.items():
        results[category] = {}

        for calc_name, methods in calculations.items():
            results[category][calc_name] = {}

            for method_name in methods.keys():
                readme_path = src_dir / category / calc_name / method_name / "README.md"

                if readme_path.exists():
                    description = extract_quick_description(readme_path)
                    rel_path = get_relative_readme_path(readme_path, project_root)
                    results[category][calc_name][method_name] = (rel_path, description)

    return results


def generate_user_md_section(results):
    """
    Generate the markdown section for User.md.

    Args:
        results: Dict from scan_category_readmes

    Returns:
        String containing the markdown content
    """
    lines = []

    for category in sorted(results.keys()):
        calculations = results[category]
        lines.append(f"#### 🔸 `{category}`")
        lines.append("")

        for calc_name in sorted(calculations.keys()):
            methods = calculations[calc_name]
            lines.append(f"- **{calc_name}**")

            for method_name in sorted(methods.keys()):
                rel_path, description = methods[method_name]

                if description:
                    lines.append(f"  - [{method_name}]({rel_path}) - {description}")
                else:
                    lines.append(f"  - [{method_name}]({rel_path})")

        lines.append("")
        lines.append("---")
        lines.append("")

    return "\n".join(lines)


def update_user_md(project_root, new_content):
    """
    Update the User.md file with the new category documentation.

    Args:
        project_root: Path to project root
        new_content: The new markdown content to insert
    """
    user_md_path = project_root / "docs" / "User.md"

    with open(user_md_path, 'r') as f:
        content = f.read()

    # Find the marker: "#### 🔸 `test` *(default)*" section ending with "---"
    # We want to insert after the "---" that follows the test section
    pattern = r'(#### 🔸 `test` \*\(default\)\*.*?---\s*\n)'
    match = re.search(pattern, content, re.DOTALL)

    if not match:
        print("Error: Could not find test section marker in User.md")
        return False

    # Find where the test section ends
    end_pos = match.end()

    # Remove any existing category content after test section
    # (everything from end of test section to end of file, except we keep it clean)
    before_content = content[:end_pos]

    # Build new content
    new_file_content = before_content + "\n" + new_content

    with open(user_md_path, 'w') as f:
        f.write(new_file_content)

    return True


def write(categories=None):
    """
    Main entry point to scan READMEs and update User.md.

    Args:
        categories: Optional categories dict. If None, imports from categories.py
    """
    if categories is None:
        from ..categories import CATEGORIES
        categories = CATEGORIES

    project_root = get_project_root()

    print("Scanning README.md files...")
    results = scan_category_readmes(categories, project_root)

    # Count found descriptions
    total = 0
    with_desc = 0
    for category, calcs in results.items():
        for calc, methods in calcs.items():
            for method, (path, desc) in methods.items():
                total += 1
                if desc:
                    with_desc += 1

    print(f"  Found {total} README files, {with_desc} with Quick Descriptions")

    print("Generating User.md content...")
    new_content = generate_user_md_section(results)

    print("Updating docs/User.md...")
    if update_user_md(project_root, new_content):
        print("✓ Successfully updated User.md")
    else:
        print("✗ Failed to update User.md")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from categories import CATEGORIES
    write(CATEGORIES)
