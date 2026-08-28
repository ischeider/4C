# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Shared schema and version helpers for the input file migration workflow.

This mirrors the (intentionally small) migration entry schema used by FourCIPP's
migration tool (see `fourcipp.migration.database` in the FourCIPP repository). It is kept
as a self-contained copy here so that 4C's CI does not need a runtime dependency on
FourCIPP just to validate a contributor's migration fragment. If FourCIPP's schema changes,
this copy should be updated to match.
"""

import pathlib
import re

import yaml

# The only three filenames contributors are allowed to use for a pending migration
# fragment; the actual `<version>.yaml` name is assigned by the `bump-input-version` script
# once the target version is known, i.e. once it is this fragment's turn to be merged.
RESERVED_FRAGMENT_NAMES = (
    "patch_upgrade.yaml",
    "minor_upgrade.yaml",
    "major_upgrade.yaml",
)

# Index into a (major, minor, patch) version tuple that each fragment name bumps.
_BUMP_INDEX = {
    "major_upgrade.yaml": 0,
    "minor_upgrade.yaml": 1,
    "patch_upgrade.yaml": 2,
}

# Fields required on every migration entry, regardless of its type.
UNIVERSAL_FIELDS = {"id", "type", "description"}

# Fields required per migration entry `type`, in addition to `UNIVERSAL_FIELDS`. Keep in
# sync with `fourcipp.migration.database._REQUIRED_FIELDS`.
REQUIRED_FIELDS_BY_TYPE: dict[str, set[str]] = {
    "parameter_removed": {"path"},
    "parameter_renamed": {"path", "new_name"},
    "parameter_default_changed": {"path", "old_default", "new_default"},
    "parameter_added": {"path", "value"},
    "section_renamed": {"path", "new_name"},
    "section_merged": {"old_path", "new_path"},
    "parameter_value_renamed": {"path", "value_map"},
    "parameter_moved": {"old_path", "new_path"},
    "type_renamed": {"path", "new_name"},
    "section_added": {"path", "value"},
}

_VERSION_PATTERN = re.compile(r"(\d+)\.(\d+)\.(\d+)")

Version = tuple[int, int, int]


class MigrationValidationError(Exception):
    """Raised when a migration fragment file fails validation."""


def parse_version(version_string: str) -> Version:
    """Parse a `MAJOR.MINOR.PATCH` version string.

    Args:
        version_string: Version string to parse

    Returns:
        The version as a `(major, minor, patch)` tuple of ints

    Raises:
        MigrationValidationError: If `version_string` is not a valid version
    """
    match = _VERSION_PATTERN.fullmatch(version_string.strip())
    if match is None:
        raise MigrationValidationError(
            f"'{version_string}' is not a valid MAJOR.MINOR.PATCH version string."
        )
    major, minor, patch = match.groups()
    return int(major), int(minor), int(patch)


def format_version(version: Version) -> str:
    """Format a version tuple as a `MAJOR.MINOR.PATCH` string.

    Args:
        version: Version to format

    Returns:
        The formatted version string
    """
    return ".".join(str(part) for part in version)


def bump_version(current: Version, fragment_name: str) -> Version:
    """Bump a version according to which reserved fragment name is used.

    A `major`/`minor` bump resets all less significant components to zero (standard semver
    behavior); a `patch` bump only increments the patch component.

    Args:
        current: Current version
        fragment_name: One of `RESERVED_FRAGMENT_NAMES`

    Returns:
        The bumped version

    Raises:
        MigrationValidationError: If `fragment_name` is not a known reserved fragment name
    """
    if fragment_name not in _BUMP_INDEX:
        raise MigrationValidationError(
            f"'{fragment_name}' is not a reserved fragment name."
        )

    index = _BUMP_INDEX[fragment_name]
    bumped = list(current)
    bumped[index] += 1
    for reset_index in range(index + 1, len(bumped)):
        bumped[reset_index] = 0
    return tuple(bumped)  # type: ignore[return-value]


def validate_migration_entry(entry: dict) -> list[str]:
    """Validate a single migration entry.

    Args:
        entry: Migration entry to validate

    Returns:
        A list of human-readable error messages; empty if the entry is valid
    """
    errors = []

    missing_universal = UNIVERSAL_FIELDS - set(entry)
    if missing_universal:
        errors.append(
            f"Entry {entry} is missing required field(s): {missing_universal}."
        )
        # Without a `type`, further checks are not meaningful.
        return errors

    entry_type = entry["type"]
    if entry_type not in REQUIRED_FIELDS_BY_TYPE:
        errors.append(
            f"Entry '{entry.get('id')}' has unknown type '{entry_type}'. Known types are: "
            f"{sorted(REQUIRED_FIELDS_BY_TYPE)}."
        )
        return errors

    missing_type_fields = REQUIRED_FIELDS_BY_TYPE[entry_type] - set(entry)
    if missing_type_fields:
        errors.append(
            f"Entry '{entry['id']}' of type '{entry_type}' is missing required field(s): "
            f"{missing_type_fields}."
        )

    return errors


def validate_migration_file(path: pathlib.Path) -> list[str]:
    """Load and validate all migration entries from a single fragment/version file.

    Args:
        path: Path to the migration YAML file

    Returns:
        A list of human-readable error messages; empty if the file is valid
    """
    try:
        data = yaml.safe_load(path.read_text())
    except yaml.YAMLError as error:
        return [f"'{path}' is not valid YAML: {error}"]

    entries = data.get("migrations") if isinstance(data, dict) else None
    if not isinstance(entries, list) or not entries:
        return [f"'{path}' must contain a non-empty top-level 'migrations' list."]

    errors = []
    for entry in entries:
        errors.extend(validate_migration_entry(entry))
    return errors
