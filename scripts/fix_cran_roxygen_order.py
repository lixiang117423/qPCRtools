#!/usr/bin/env python3
"""Move CRAN-package globalVariables() blocks ahead of roxygen comments."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R_DIR = ROOT / "CRAN" / "qPCRtools" / "R"

FILES = [
    "CalCurve.R",
    "CalExp2dCt.R",
    "CalExp2ddCt.R",
    "CalExpCurve.R",
    "CalExpRqPCR.R",
    "CalRTable.R",
]


def matching_close(text: str, start: int) -> int:
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    raise ValueError("unbalanced globalVariables() block")


def move_block(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    start = text.find("globalVariables(")
    if start < 0:
        raise ValueError(f"globalVariables() not found: {path}")
    first_function = text.find(" <- function(")
    if first_function < 0 or start > first_function:
        raise ValueError(f"globalVariables() is not before the first function: {path}")

    close = matching_close(text, start + len("globalVariables"))
    end = close + 1
    if text[end : end + 2] == "\r\n":
        end += 2
    elif text[end : end + 1] == "\n":
        end += 1

    block = text[start:end].rstrip()
    remaining = text[:start] + text[end:]
    if not remaining.lstrip().startswith("#'"):
        raise ValueError(f"roxygen block is not adjacent after moving block: {path}")
    path.write_text(block + "\n\n" + remaining, encoding="utf-8")


def main() -> None:
    for name in FILES:
        path = R_DIR / name
        move_block(path)
        print(f"reordered {path}")


if __name__ == "__main__":
    main()
