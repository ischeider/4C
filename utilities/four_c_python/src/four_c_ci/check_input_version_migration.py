# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Validate input file migration fragments under `input_version/migrations/`.

Run as a pre-commit hook on every pull request (see `.pre-commit-config.yaml`). Validates:

- Every migration YAML file's entries against the shared schema.
- That at most one of the three reserved, generically named fragments (`patch_upgrade.yaml`,
  `minor_upgrade.yaml`, `major_upgrade.yaml`) is present at once, since only one pending,
  not-yet-versioned change is allowed at a time (the actual `<version>.yaml` name is
  assigned later, when the fragment is merged, by `bump-input-version`).
- That only the three reserved names (or an already-assigned `<version>.yaml` file) are used;
  a contributor must never hand-pick a concrete version number.
"""

import argparse
import pathlib
import sys

from four_c_ci.input_version_migration_schema import (
    RESERVED_FRAGMENT_NAMES,
    parse_version,
    validate_migration_file,
)


def check_files(paths: list[pathlib.Path]) -> list[str]:
    """Validate a list of migration YAML files.

    Args:
        paths: Paths to migration YAML files (fragments or already-assigned version files)

    Returns:
        A list of human-readable error messages; empty if all files are valid
    """
    errors = []
    reserved_present = []

    for path in paths:
        if path.name in RESERVED_FRAGMENT_NAMES:
            reserved_present.append(path)
        else:
            try:
                parse_version(path.stem)
            except Exception:
                errors.append(
                    f"'{path}' is neither a reserved fragment name ({', '.join(RESERVED_FRAGMENT_NAMES)}) "
                    "nor a valid '<MAJOR>.<MINOR>.<PATCH>.yaml' version file. Contributors must "
                    "use one of the reserved names; the concrete version is assigned automatically."
                )

        errors.extend(validate_migration_file(path))

    if len(reserved_present) > 1:
        names = ", ".join(str(path) for path in reserved_present)
        errors.append(
            "Only one pending migration fragment is allowed at a time, but found "
            f"{len(reserved_present)}: {names}. Combine your changes into a single fragment."
        )

    return errors


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Validate input file migration fragments/version files.",
    )
    parser.add_argument(
        "files", nargs="*", type=pathlib.Path, help="Migration YAML files"
    )
    args = parser.parse_args()

    if not args.files:
        return

    errors = check_files(args.files)

    if errors:
        for error in errors:
            print(f"error: {error}", file=sys.stderr)
        sys.exit(1)

    for path in args.files:
        print(f"OK: {path}")


if __name__ == "__main__":
    main()
