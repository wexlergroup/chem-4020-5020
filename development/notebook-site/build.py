"""Build the notebook site of the course: the Python cells of the notes, editable in the reader's browser.

Every section of a module that has code cells becomes one notebook. A notebook holds the code cells of its page,
each under the heading and the figure caption that it has in the notes. The notebooks, the figure package
(book/figstyle) and the data files (book/_static) are then built into a JupyterLite site, in which Python runs
in the browser.

Run it after the book is built, because figure numbers, equation numbers and citations are read from the build:

    cd book && myst build --execute --html && cd ..
    python -m venv .lite && .lite/bin/pip install -r development/notebook-site/requirements.txt
    .lite/bin/python development/notebook-site/build.py            # writes book/_build/html/lite

The separate environment keeps the site to JupyterLite and its Python kernel: the build adds every Jupyter
extension that it finds in its own environment.

The site is then served at <address of the book>/lite/. A page links to its notebook as
<address of the book>/lite/notebooks/index.html?path=chapter-01/section-03.ipynb
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SITE = "https://wexlergroup.github.io/chem-4020-5020"
DATA = {".csv", ".tsv", ".txt", ".dat", ".json"}

CELL = re.compile(r"^```\{code-cell\}[^\n]*\n(.*?)^```[ \t]*$", re.S | re.M)
TITLE = re.compile(r"^# +(.+?)\s*$", re.M)
HEADING = re.compile(r"^#{2,4} +(.+?)\s*$", re.M)
FENCE = re.compile(r"^(```|:::).*?^\1[ \t]*$", re.S | re.M)

INTRO = """# Code for {title}

This notebook holds the Python cells of [{title}]({page}) in the lecture notes, each under the heading and \
the figure caption that it has there. It runs in your browser, and nothing is installed on your computer.

- **Run a cell.** Click the cell and press Shift+Enter. The first time, run the cells in order, or choose Run and \
then Run All Cells. Python needs about half a minute to start.
- **Explore.** Change a number in a cell and run the cell again.
- **Start over.** Your changes are saved in this browser only. To get the original notebook back, open the \
[file list]({site}/lite/tree/index.html), right-click this notebook, choose Delete, and open the notebook again \
from the notes."""


def text_of(node):
    """Plain text of a node of the built page."""
    if node.get("type") == "text":
        return node.get("value", "")
    return "".join(text_of(child) for child in node.get("children", []))


def read_build(page):
    """Numbers and citations of a built page: figure label -> number, reference label -> text, citation -> text."""
    figures, references, citations = {}, {}, {}

    def walk(node):
        if isinstance(node, list):
            for item in node:
                walk(item)
        elif isinstance(node, dict):
            kind, label = node.get("type"), node.get("label")
            if kind == "container" and node.get("kind") == "figure" and label and node.get("enumerator"):
                figures[label] = node["enumerator"]
            elif kind == "crossReference" and label:
                references[label] = text_of(node)
            elif kind == "cite" and label:
                citations[label, node.get("kind")] = text_of(node)
            for value in node.values():
                walk(value)

    walk(json.loads(page.read_text())["mdast"])
    return figures, references, citations


def caption(text, source, references, citations):
    """A caption of the notes as Markdown that a notebook can show."""
    def reference(match):
        label = match.group(1)
        if label not in references:
            raise SystemExit(f"{source}: the build has no number for #{label}")
        return references[label]

    def citation(match):
        key = match.group(2)
        kind = "parenthetical" if match.group(1) else "narrative"
        if (key, kind) not in citations:
            raise SystemExit(f"{source}: the build has no citation for @{key}")
        return f"({citations[key, kind]})" if match.group(1) else citations[key, kind]

    def link(match):
        target = (source.parent / match.group(2)).resolve().relative_to(source.parents[1].resolve())
        return f"[{match.group(1)}]({SITE}/{target.with_suffix('').as_posix()})"

    text = re.sub(r"\[\]\(#([^)]+)\)", reference, text)
    text = re.sub(r"(\[)?@([A-Za-z][\w-]*)(?(1)\])", citation, text)
    return re.sub(r"\[([^\]]+)\]\((?!https?:|#)([^)]+\.md)\)", link, text)


def notebook(source, built):
    """The notebook of one page of the notes, or None if the page has no code."""
    text = source.read_text()
    if not CELL.search(text):
        return None
    figures, references, citations = read_build(built)
    title = TITLE.search(text).group(1)
    page = f"{SITE}/{source.relative_to(source.parents[1]).with_suffix('').as_posix()}"
    cells = [("markdown", INTRO.format(title=title, page=page, site=SITE))]
    position, shown = 0, None
    for match in CELL.finditer(text):
        lines = match.group(1).split("\n")
        options = {}
        while lines and lines[0].startswith(":"):
            key, _, value = lines.pop(0)[1:].partition(":")
            options[key] = value.strip()
        between = FENCE.sub("", text[position:match.start()])       # headings only, not "#" lines inside blocks
        headings = HEADING.findall(between)
        notes = []
        if headings and headings[-1] != shown:
            shown = headings[-1]
            notes.append(f"## {shown}")
        if "caption" in options:
            label = options.get("label")
            if label not in figures:
                raise SystemExit(f"{source}: the build has no figure number for {label}")
            notes.append(f"**Figure {figures[label]}.** " + caption(options["caption"], source, references, citations))
        if notes:
            cells.append(("markdown", "\n\n".join(notes)))
        cells.append(("code", "\n".join(lines).strip("\n")))
        position = match.end()
    return {
        "cells": [{"cell_type": "markdown", "metadata": {}, "source": body} if kind == "markdown" else
                  {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": body}
                  for kind, body in cells],
        "metadata": {"kernelspec": {"display_name": "Python (Pyodide)", "language": "python", "name": "python"},
                     "language_info": {"name": "python"}},
        "nbformat": 4, "nbformat_minor": 5,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--book", type=Path, default=REPO / "book", help="the book folder (default: book)")
    parser.add_argument("--out", type=Path, default=None, help="where the site goes (default: book/_build/html/lite)")
    parser.add_argument("--notebooks-only", type=Path, default=None, metavar="FOLDER",
                        help="write the notebooks, figstyle and the data to FOLDER and do not build the site")
    args = parser.parse_args()
    book = args.book.resolve()
    out = (args.out or book / "_build" / "html" / "lite").resolve()
    content = book / "_build" / "site" / "content"
    if not content.is_dir():
        raise SystemExit(f"{content} does not exist: build the book first (myst build --execute --html)")

    with tempfile.TemporaryDirectory() as scratch:
        files = args.notebooks_only.resolve() if args.notebooks_only else Path(scratch) / "files"
        count = 0
        for built in sorted(content.glob("*.json")):
            source = book / json.loads(built.read_text())["location"].lstrip("/")
            if source.suffix != ".md" or not source.parent.name.startswith("chapter-"):
                continue                                   # the sections of the modules; not the appendices
            made = notebook(source, built)
            if made:
                target = files / source.relative_to(book).with_suffix(".ipynb")
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(json.dumps(made, indent=1, ensure_ascii=False) + "\n")
                count += 1
        shutil.copytree(book / "figstyle", files / "figstyle", ignore=shutil.ignore_patterns("__pycache__"),
                        dirs_exist_ok=True)
        for data in sorted((book / "_static").rglob("*")):
            if data.suffix in DATA and data.is_file():
                target = files / data.relative_to(book)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(data, target)
        print(f"{count} notebooks")
        if args.notebooks_only:
            return
        if out.exists():
            shutil.rmtree(out)
        subprocess.run([sys.executable, "-m", "jupyterlite_core", "build", "--contents", str(files),
                        "--output-dir", str(out)], cwd=scratch, check=True)
    print(f"notebook site written to {out}")


if __name__ == "__main__":
    main()
