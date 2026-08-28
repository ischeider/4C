# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Assign the real version number to a pending input file migration fragment.

Run once per merge to `main` (see `.github/workflows/input_version_bump.yml`), after the
merge queue has serialized the actual merge. Looks for one of the three reserved fragment
names (`patch_upgrade.yaml`, `minor_upgrade.yaml`, `major_upgrade.yaml`) under
`input_version/migrations/`; if present, renames it to `<new_version>.yaml` and updates
`input_version/INPUT_VERSION` accordingly. If no fragment is present, this is a no-op (the
merged PR did not require an input file migration).

Prints a short summary and, if `GITHUB_OUTPUT` is set, writes `bumped`, `old_version`,
`new_version`, `fragment` and `new_filename` outputs for the calling workflow to use (e.g. to
craft the commit message).
"""

import argparse
import os
import pathlib
import sys

from four_c_ci.input_version_migration_schema import (
    RESERVED_FRAGMENT_NAMES,
    bump_version,
    format_version,
    parse_version,
    validate_migration_file,
)

MIGRATIONS_DIR = pathlib.Path("input_version/migrations")
VERSION_FILE = pathlib.Path("input_version/INPUT_VERSION")


def _write_github_output(outputs: dict[str, str]) -> None:
    """Append key-value pairs to the `GITHUB_OUTPUT` file, if set.

    Args:
        outputs: Mapping of output name to value
    """
    github_output = os.environ.get("GITHUB_OUTPUT")
    if not github_output:
        return
    with open(github_output, "a") as output_file:
        for key, value in outputs.items():
            output_file.write(f"{key}={value}\n")


def bump(repo_root: pathlib.Path) -> bool:
    """Rename a pending migration fragment and bump `INPUT_VERSION`, if one exists.

    Args:
        repo_root: Path to the root of the 4C repository

    Returns:
        True if a fragment was found and the version was bumped, False otherwise

    Raises:
        SystemExit: If more than one reserved fragment is present, or the fragment fails
            schema validation (both should already be caught by `check-input-version-migration`
            on the pull request, so this indicates the check was bypassed)
    """
    migrations_dir = repo_root / MIGRATIONS_DIR
    version_file = repo_root / VERSION_FILE

    present = [
        migrations_dir / name
        for name in RESERVED_FRAGMENT_NAMES
        if (migrations_dir / name).is_file()
    ]

    if not present:
        print("No pending migration fragment found; nothing to bump.")
        _write_github_output({"bumped": "false"})
        return False

    if len(present) > 1:
        names = ", ".join(str(path) for path in present)
        sys.exit(
            f"error: Found {len(present)} pending migration fragments, expected at most one: "
            f"{names}."
        )

    fragment = present[0]
    errors = validate_migration_file(fragment)
    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        sys.exit(
            f"error: '{fragment}' failed validation, refusing to bump the version."
        )

    current_version = parse_version(version_file.read_text().strip())
    new_version = bump_version(current_version, fragment.name)
    new_version_string = format_version(new_version)
    new_path = migrations_dir / f"{new_version_string}.yaml"

    if new_path.exists():
        sys.exit(f"error: Target migration file '{new_path}' already exists.")

    fragment.rename(new_path)
    version_file.write_text(new_version_string)

    old_version_string = format_version(current_version)
    print(
        f"Bumped input version: {old_version_string} -> {new_version_string} "
        f"(renamed {fragment.name} -> {new_path.name})."
    )
    _write_github_output(
        {
            "bumped": "true",
            "old_version": old_version_string,
            "new_version": new_version_string,
            "fragment": fragment.name,
            "new_filename": new_path.name,
        }
    )
    return True


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=pathlib.Path,
        default=pathlib.Path("."),
        help="Path to the root of the 4C repository (default: current directory)",
    )
    args = parser.parse_args()

    bump(args.repo_root)


if __name__ == "__main__":
    main()
